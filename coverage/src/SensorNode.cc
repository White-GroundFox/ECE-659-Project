#include "SensorNode.h"
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <cmath>

Define_Module(SensorNode);

// ════════════════════════════════════════════════════════════════════════════
//  Tuning constants
// ════════════════════════════════════════════════════════════════════════════

// Pheromone grid (Phase 1 — cascade bias)
static constexpr double PHEROMONE_DEPOSIT      = 1.00;  // mark intensity per activation
static constexpr double PHEROMONE_DECAY        = 0.80;  // per-round grid decay factor
static constexpr double PHEROMONE_SENSE_THRESH = 0.25;  // cell counts as "covered" above this
static constexpr double COVERAGE_BIAS          = 0.25;  // max extra timer penalty (× Ts)

// Cascade edge weights (Steps 3 & 4)
static constexpr double EDGE_REINFORCE    = 0.50;
static constexpr double EDGE_DECAY        = 0.85;
static constexpr double EDGE_PRUNE_THRESH = 0.10;

// Phase 2 perimeter check (post-selection pruning)
static constexpr int    PERIM_SAMPLES = 36;   // points checked around sensing boundary


// ════════════════════════════════ Initialisation ═════════════════════════════

void SensorNode::initialize(int stage)
{
    if (stage == 0) {
        nodeId     = getIndex();
        numNodes   = getParentModule()->par("numSensors");
        rs         = par("sensingRange");
        rc         = par("radioRange");
        areaSize   = par("areaSize");
        initEnergy = par("initialEnergy");
        energy     = initEnergy;
        roundTime  = par("roundTime");
        Td = par("Td");
        Ts = par("Ts");
        Te = par("Te");

        posX = uniform(0, areaSize);
        posY = uniform(0, areaSize);

        state           = UNDECIDED;
        hasPred         = false;
        cascadeParentId = -1;

        // Pheromone grid: 10×10 cells (5m each for a 50m area).
        // Coarse enough to be fast; fine enough to distinguish coverage regions.
        gridRes  = 10;
        cellSize = areaSize / gridRes;
        stigGrid.assign(gridRes, std::vector<double>(gridRes, 0.0));

        sigCoverage = registerSignal("coverage");
        sigActive   = registerSignal("activeNodes");

        updateDisplay();

    } else {
        cModule *parent = getParentModule();
        for (int i = 0; i < numNodes; i++) {
            if (i == nodeId) continue;
            SensorNode *s = check_and_cast<SensorNode*>(
                                parent->getSubmodule("sensor", i));
            if (dist(posX, posY, s->posX, s->posY) <= rc)
                neighbours[i] = {i, s->posX, s->posY};
        }
        EV_INFO << "[Node " << nodeId << "] pos=("
                << posX << "," << posY
                << ")  neighbours=" << neighbours.size() << "\n";

        roundTimer = new cMessage("round", KIND_ROUND);
        scheduleAt(simTime(), roundTimer);
    }
}

// ════════════════════════════════ Message Handler ════════════════════════════

void SensorNode::handleMessage(cMessage *msg)
{
    if (msg->isSelfMessage()) {
        switch (msg->getKind()) {

        // ── New round ──────────────────────────────────────────────────────
        case KIND_ROUND: {
            delete msg; roundTimer = nullptr;
            if (roundCount > 0) {
                double dt = roundTime - Ts;
                energy -= (state == ON) ? dt * 4.0 / 100.0
                                        : dt * 0.01 / 100.0;
                energy = std::max(0.0, energy);
            }
            roundCount++;
            startRound();
            roundTimer = new cMessage("round", KIND_ROUND);
            scheduleAt(simTime() + roundTime, roundTimer);
            break;
        }

        // ── Volunteer as starting node ─────────────────────────────────────
        case KIND_VOLUNTEER:
            delete msg; volunteerTimer = nullptr;
            if (state == UNDECIDED) {
                EV_INFO << "[Node " << nodeId
                        << "] volunteers (energy=" << energy << ")\n";
                hasPred         = false;
                cascadeParentId = -1;
                turnOn();
            }
            break;

        // ── Selected by cascade ────────────────────────────────────────────
        case KIND_SELECT:
            delete msg; selectTimer = nullptr;
            if (state == UNDECIDED) {
                EV_INFO << "[Node " << nodeId << "] selected by cascade\n";
                hasPred         = bestHasPred;
                predX           = bestPredX;
                predY           = bestPredY;
                cascadeParentId = bestPredId;
                turnOn();
            }
            break;

        // ── End of selection phase: remaining UNDECIDED → OFF ─────────────
        case KIND_STEADY:
            delete msg; steadyTimer = nullptr;
            if (state == UNDECIDED) turnOff();
            break;

        // ── Phase 2: post-selection perimeter pruning ──────────────────────
        // Fires at Ts+Te after cascade is fully settled.
        // Each ON node checks whether its entire sensing perimeter is already
        // covered by known active neighbours. If yes → redundant → sleep.
        // This is the core fix: the perimeter check runs when knownActive is
        // fully populated (all CoverageMarkMsgs received), not during the race.
        case KIND_REDUNDANCY_CHECK: {
            delete msg; redundancyTimer = nullptr;
            if (state == ON && isPerimeterRedundant()) {
                EV_INFO << "[Node " << nodeId
                        << "] Phase-2 pruned: perimeter fully covered by "
                        << knownActive.size() << " known active neighbours\n";
                turnOff();
            }
            break;
        }

        // ── Step 4: Pheromone decay + edge pruning ─────────────────────────
        case KIND_PRUNE:
            delete msg; pruneTimer = nullptr;
            decayAndPruneEdges();
            break;

        // ── Measure coverage ───────────────────────────────────────────────
        case KIND_STATS: {
            delete msg; statsTimer = nullptr;
            cModule *par = getParentModule();
            int n = par->par("numSensors");
            int minAlive = -1;
            for (int i = 0; i < n; i++) {
                SensorNode *s = check_and_cast<SensorNode*>(
                                    par->getSubmodule("sensor", i));
                if (s->energy > 0.0) { minAlive = i; break; }
            }
            if (nodeId == minAlive) computeCoverage();
            break;
        }

        default: delete msg;
        }

    } else {
        int kind = msg->getKind();
        if      (kind == 100) handlePowerOn(check_and_cast<PowerOnMsg*>(msg));
        else if (kind == 101) handleCoverageMark(check_and_cast<CoverageMarkMsg*>(msg));
        else                  delete msg;
    }
}

// ════════════════════════════════ Round Management ═══════════════════════════

void SensorNode::startRound()
{
    cancelTimer(volunteerTimer);
    cancelTimer(selectTimer);
    cancelTimer(steadyTimer);
    cancelTimer(statsTimer);
    cancelTimer(pruneTimer);
    cancelTimer(redundancyTimer);

    state           = UNDECIDED;
    hasPred         = false;
    bestTime        = 1e9;
    bestHasPred     = false;
    bestPredId      = -1;
    cascadeParentId = -1;

    // Clear Phase 2 knowledge — will be rebuilt from CoverageMarkMsgs
    knownActive.clear();
    updateDisplay();

    // Volunteer timer: energy-weighted
    double lostFrac = 1.0 - energy / initEnergy;
    lostFrac = std::max(0.0, std::min(lostFrac, 0.9));
    double vt = exponential(0.05 + lostFrac * 0.40);
    vt = std::min(vt, Ts * 0.9);
    volunteerTimer = new cMessage("vol", KIND_VOLUNTEER);
    scheduleAt(simTime() + vt, volunteerTimer);

    // Steady timer: end of Phase 1 (cascade selection)
    steadyTimer = new cMessage("steady", KIND_STEADY);
    scheduleAt(simTime() + Ts, steadyTimer);

    // Phase 4: pheromone decay fires at mid-round
    pruneTimer = new cMessage("prune", KIND_PRUNE);
    scheduleAt(simTime() + roundTime * 0.5, pruneTimer);

    // Stats: near end of steady phase
    statsTimer = new cMessage("stats", KIND_STATS);
    scheduleAt(simTime() + roundTime - 1.0, statsTimer);

    // Phase 2 redundancy check: fires at Ts + Te + jitter for ALL nodes.
    // The random jitter (uniform across Te) staggers the checks so nodes
    // decide sequentially rather than simultaneously.  Early-checking nodes
    // turn off; later-checking nodes query the live network state and
    // correctly see fewer active neighbours → stay ON if needed.
    redundancyTimer = new cMessage("redcheck", KIND_REDUNDANCY_CHECK);
    scheduleAt(simTime() + Ts + Te + uniform(0, Te), redundancyTimer);
}

// ════════════════════════════════ OGDC + Stigmergy Phase 1 ═══════════════════

void SensorNode::handlePowerOn(PowerOnMsg *m)
{
    if (state != UNDECIDED) { delete m; return; }

    // ── Step 1: Update pheromone grid from environmental marks ────────────
    markCoverage(m->senderX, m->senderY);
    if (m->hasPred && m->predX >= 0)
        markCoverage(m->predX, m->predY);

    // Cancel volunteer timer — a neighbour is already seeding this area
    cancelTimer(volunteerTimer);

    double selectTime;
    double myPredX   = m->senderX;
    double myPredY   = m->senderY;
    bool   myHasPred = true;

    if (!m->hasPred) {
        // Theorem 1: ideal B is at distance √3·rs from A
        double dA  = dist(posX, posY, m->senderX, m->senderY);
        if (dA < 1e-6) { delete m; return; }
        double dev = std::fabs(dA - std::sqrt(3.0) * rs);
        selectTime = Td + (dev / rs) * Ts;
    } else {
        // Theorem 2: next node closest to target positions T1, T2
        double dT1 = (m->t1x >= 0) ? dist(posX, posY, m->t1x, m->t1y) : 1e9;
        double dT2 = (m->t2x >= 0) ? dist(posX, posY, m->t2x, m->t2y) : 1e9;
        double best = std::min(dT1, dT2);
        if (best > rs) { delete m; return; }
        selectTime = Td + (best / rs) * Ts;
    }

    // ── Step 2: Bias timer AWAY from pheromone-saturated areas ───────────
    // Read the pheromone level at this node's position.
    // This reflects activations from PREVIOUS rounds (current round marks are
    // too sparse to be useful during the cascade phase).
    // Nodes in freshly uncovered territory → short timer → win the race.
    // Nodes near previously active positions → longer timer → yield to fresher nodes.
    double cov = localCoverage(posX, posY);
    selectTime += cov * COVERAGE_BIAS * Ts;

    EV_INFO << "[Node " << nodeId << "] pheromone="
            << cov * 100.0 << "% selectTime=" << selectTime << "s\n";

    if (selectTime < bestTime) {
        bestTime    = selectTime;
        bestHasPred = myHasPred;
        bestPredX   = myPredX;
        bestPredY   = myPredY;
        bestPredId  = m->senderId;
        cancelTimer(selectTimer);
        selectTimer = new cMessage("sel", KIND_SELECT);
        scheduleAt(simTime() + selectTime, selectTimer);
    }
    delete m;
}

void SensorNode::turnOn()
{
    state = ON;
    cancelTimer(volunteerTimer);
    cancelTimer(selectTimer);
    energy -= 20.0 * Td / 100.0;
    energy  = std::max(0.0, energy);
    updateDisplay();

    // ── Step 1: Mark this position in local pheromone grid ────────────────
    markCoverage(posX, posY);

    // ── Step 3: Reinforce the cascade edge that led to this activation ────
    if (cascadeParentId >= 0) {
        edgeWeight[cascadeParentId] += EDGE_REINFORCE;
        EV_INFO << "[Node " << nodeId << "] reinforced edge from node "
                << cascadeParentId
                << " weight=" << edgeWeight[cascadeParentId] << "\n";
    }

    // Broadcast OGDC cascade message + stigmergic coverage mark
    broadcastPowerOn();
    broadcastCoverageMark();
}

void SensorNode::turnOff()
{
    state = OFF;
    updateDisplay();
}

void SensorNode::broadcastPowerOn()
{
    double t1x=-1, t1y=-1, t2x=-1, t2y=-1;
    if (hasPred)
        computeTargets(predX, predY, posX, posY, t1x, t1y, t2x, t2y);

    PowerOnMsg tmpl;
    tmpl.senderId = nodeId;
    tmpl.senderX  = posX;  tmpl.senderY = posY;
    tmpl.hasPred  = hasPred;
    tmpl.predX    = predX; tmpl.predY   = predY;
    tmpl.predId   = (hasPred && cascadeParentId >= 0) ? cascadeParentId : -1;
    tmpl.t1x = t1x; tmpl.t1y = t1y;
    tmpl.t2x = t2x; tmpl.t2y = t2y;

    for (auto &kv : neighbours) {
        cModule *mod = getParentModule()->getSubmodule("sensor", kv.first);
        sendDirect(tmpl.dup(), mod, "radioIn");
    }
}

void SensorNode::broadcastCoverageMark()
{
    CoverageMarkMsg tmpl;
    tmpl.senderId     = nodeId;
    tmpl.senderX      = posX;
    tmpl.senderY      = posY;
    tmpl.sensingRange = rs;

    for (auto &kv : neighbours) {
        cModule *mod = getParentModule()->getSubmodule("sensor", kv.first);
        sendDirect(tmpl.dup(), mod, "radioIn");
    }
}

// Receives a CoverageMarkMsg from an active neighbour.
// Serves both phases:
//   Phase 1: updates local pheromone grid (future round bias)
//   Phase 2: populates knownActive (used by perimeter pruning at Ts+Te)
void SensorNode::handleCoverageMark(CoverageMarkMsg *m)
{
    // Phase 1: update pheromone
    markCoverage(m->senderX, m->senderY, PHEROMONE_DEPOSIT);

    // Phase 2: record active neighbour position for perimeter check
    bool found = false;
    for (auto &k : knownActive) {
        if (std::fabs(k.first - m->senderX) < 1e-6 &&
            std::fabs(k.second - m->senderY) < 1e-6) { found = true; break; }
    }
    if (!found)
        knownActive.push_back({m->senderX, m->senderY});

    delete m;
}

// ════════════════════════════════ Geometry ═══════════════════════════════════

void SensorNode::computeTargets(double ax, double ay,
                                 double bx, double by,
                                 double &t1x, double &t1y,
                                 double &t2x, double &t2y)
{
    t1x = t1y = t2x = t2y = -1;
    double dAB = dist(ax, ay, bx, by);
    if (dAB < 1e-6 || dAB > 2.0*rs) return;

    double mx = (ax+bx)*0.5, my = (ay+by)*0.5;
    double ux = (bx-ax)/dAB, uy = (by-ay)/dAB;
    double px = -uy,          py =  ux;
    double h  = std::sqrt(rs*rs - (dAB*0.5)*(dAB*0.5));

    double o1x = mx + h*px, o1y = my + h*py;
    double o2x = mx - h*px, o2y = my - h*py;

    auto makeTarget = [&](double ox, double oy, double &tx, double &ty) {
        double dx = ox - mx, dy = oy - my;
        double len = std::sqrt(dx*dx + dy*dy);
        if (len < 1e-9) { tx = ty = -1; return; }
        tx = ox + (dx/len)*rs;
        ty = oy + (dy/len)*rs;
    };
    makeTarget(o1x, o1y, t1x, t1y);
    makeTarget(o2x, o2y, t2x, t2y);
}

// ════════════════════════════════ Stigmergy Phase 1 ══════════════════════════

// Deposit pheromone on all grid cells within rs of (cx, cy)
void SensorNode::markCoverage(double cx, double cy, double deposit)
{
    for (int gx = 0; gx < gridRes; gx++) {
        for (int gy = 0; gy < gridRes; gy++) {
            double ccx = (gx + 0.5) * cellSize;
            double ccy = (gy + 0.5) * cellSize;
            if (dist(ccx, ccy, cx, cy) <= rs)
                stigGrid[gx][gy] = std::min(1.0, stigGrid[gx][gy] + deposit);
        }
    }
}

// Fraction of cells within rs of (cx, cy) that exceed the sensing threshold.
// Range: 0.0 (no marks nearby) → 1.0 (fully saturated from previous rounds).
double SensorNode::localCoverage(double cx, double cy) const
{
    int total = 0, marked = 0;
    for (int gx = 0; gx < gridRes; gx++) {
        for (int gy = 0; gy < gridRes; gy++) {
            double ccx = (gx + 0.5) * cellSize;
            double ccy = (gy + 0.5) * cellSize;
            if (dist(ccx, ccy, cx, cy) <= rs) {
                total++;
                if (stigGrid[gx][gy] >= PHEROMONE_SENSE_THRESH) marked++;
            }
        }
    }
    return (total > 0) ? (double)marked / total : 0.0;
}

// Step 4: decay pheromone grid + prune weak cascade edges
void SensorNode::decayAndPruneEdges()
{
    // Evaporate pheromone — areas not recently activated fade over rounds
    for (int gx = 0; gx < gridRes; gx++)
        for (int gy = 0; gy < gridRes; gy++)
            stigGrid[gx][gy] *= PHEROMONE_DECAY;

    // Decay and prune cascade edges
    std::vector<int> toPrune;
    for (auto &kv : edgeWeight) {
        kv.second *= EDGE_DECAY;
        if (kv.second < EDGE_PRUNE_THRESH)
            toPrune.push_back(kv.first);
    }
    for (int id : toPrune) {
        EV_INFO << "[Node " << nodeId << "] pruned edge to node " << id << "\n";
        edgeWeight.erase(id);
    }
}

// ════════════════════════════════ Stigmergy Phase 2 ══════════════════════════

// 36-point perimeter check — runs AFTER the cascade has fully settled.
// By the time this fires (Ts + Te), knownActive contains positions of ALL
// active neighbours (populated from CoverageMarkMsgs received during 0..Ts+Te).
//
// For each sample point on the sensing boundary of this node:
//   - If no known active neighbour covers it → this node is NOT redundant
// If ALL 36 points are covered by known neighbours → this node is redundant
// and can safely turn OFF without losing any coverage.
bool SensorNode::isPerimeterRedundant() const
{
    // Query CURRENT live ON nodes from the network instead of using the
    // stale knownActive cache.  This is essential: nodes that already turned
    // off during Phase 2 (because they fired their staggered timer earlier)
    // must NOT be counted as coverage — the cache would still list them.
    cModule *parent = getParentModule();
    int n = parent->par("numSensors");

    std::vector<std::pair<double,double>> liveActive;
    for (int i = 0; i < n; i++) {
        if (i == nodeId) continue;   // don't count ourselves
        SensorNode *s = check_and_cast<SensorNode*>(
                            parent->getSubmodule("sensor", i));
        if (s->state == ON)
            liveActive.push_back({s->posX, s->posY});
    }
    if (liveActive.empty()) return false;

    const double TWO_PI = 2.0 * M_PI;
    for (int i = 0; i < PERIM_SAMPLES; i++) {
        double angle = TWO_PI * i / PERIM_SAMPLES;
        double px = posX + rs * std::cos(angle);
        double py = posY + rs * std::sin(angle);

        // Clamp to deployment area
        px = std::max(0.0, std::min(areaSize, px));
        py = std::max(0.0, std::min(areaSize, py));

        bool covered = false;
        for (auto &a : liveActive) {
            if (dist(px, py, a.first, a.second) <= rs) {
                covered = true;
                break;
            }
        }
        if (!covered) return false;   // uncovered perimeter point → stay ON
    }
    return true;   // all 36 perimeter points covered by live nodes → sleep
}

// ════════════════════════════════ Coverage & Display ═════════════════════════

void SensorNode::computeCoverage()
{
    cModule *parent = getParentModule();
    int n = parent->par("numSensors");

    std::vector<std::pair<double,double>> active;
    int onCount = 0, offCount = 0, deadCount = 0;

    for (int i = 0; i < n; i++) {
        SensorNode *s = check_and_cast<SensorNode*>(
                            parent->getSubmodule("sensor", i));
        if (s->energy <= 0.05 * s->initEnergy) {
            deadCount++;                         // count as dead, not ON or OFF
        } else if (s->state == ON) {
            active.push_back({s->posX, s->posY});
            onCount++;
        } else if (s->state == OFF) {
            offCount++;
        }
    }

    // Fine grid: 5 cells per metre → 250×250 for a 50m area (0.04m² per cell)
    const int    cellsPerMeter = 5;
    const int    G             = (int)(areaSize * cellsPerMeter);
    const double cs            = 1.0 / cellsPerMeter;
    int covered = 0;
    for (int gx = 0; gx < G; gx++) {
        for (int gy = 0; gy < G; gy++) {
            double cx = (gx + 0.5) * cs, cy = (gy + 0.5) * cs;
            for (auto &p : active) {
                if (dist(p.first, p.second, cx, cy) <= rs) { covered++; break; }
            }
        }
    }
    double cov = (double)covered / ((long long)G * G);

    EV_INFO << "\n========================================\n"
            << " Round        : " << roundCount          << "\n"
            << " Coordinator  : Node " << nodeId         << "\n"
            << " Active (ON)  : " << onCount             << "\n"
            << " Sleeping(OFF): " << offCount            << "\n"
            << " Dead/Low bat : " << deadCount << "\n"
            << " Grid         : " << G << "x" << G
                                  << " (" << cs << "m)\n"
            << " Coverage     : " << cov * 100.0         << "%\n"
            << " Sim time     : " << simTime()           << "s\n"
            << "========================================\n\n";

    std::cout << "\n========================================"  << std::endl;
    std::cout << " Round        : " << roundCount             << std::endl;
    std::cout << " Coordinator  : Node " << nodeId            << std::endl;
    std::cout << " Active (ON)  : " << onCount                << std::endl;
    std::cout << " Sleeping(OFF): " << offCount               << std::endl;
    std::cout << " Dead/Low bat : " << deadCount << std::endl;
    std::cout << " Grid         : " << G << "x" << G
              << " (" << cs << "m)"                          << std::endl;
    std::cout << " Coverage     : " << cov * 100.0 << "%"    << std::endl;
    std::cout << " Sim time     : " << simTime() << "s"      << std::endl;
    std::cout << "========================================"   << std::endl;

    emit(sigCoverage, cov);
    emit(sigActive,   (double)onCount);
}

void SensorNode::updateDisplay()
{
    if (!hasGUI()) return;
    const double SCALE = 10.0;
    char buf[256];
    if (energy <= 0.05 * initEnergy) {
        // Dead or low battery — solid black dot
        std::snprintf(buf, sizeof(buf),
            "p=%.0f,%.0f;b=10,10,oval,#000000,#000000,1",
            posX * SCALE, posY * SCALE);
    } else if (state == ON) {
        double ringDiam = 2.0 * rs * SCALE;
        std::snprintf(buf, sizeof(buf),
                "p=%.0f,%.0f;b=10,10,oval,#00cc00,#006600,1;r=%.0f,#00cc00",
                posX * SCALE, posY * SCALE,
                rs * SCALE);
    } else if (state == OFF) {
        std::snprintf(buf, sizeof(buf),
            "p=%.0f,%.0f;b=8,8,oval,#cc0000,black,1",
            posX * SCALE, posY * SCALE);
    } else {
        std::snprintf(buf, sizeof(buf),
            "p=%.0f,%.0f;b=8,8,oval,#aaaaaa,black,1",
            posX * SCALE, posY * SCALE);
    }
    getDisplayString().parse(buf);
}

void SensorNode::finish()
{
    if (nodeId == 0) {
        std::cout << "\n>>> SIMULATION FINISHED after "
                  << roundCount << " rounds <<<" << std::endl;
        computeCoverage();
    }
}
