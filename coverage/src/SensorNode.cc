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
static constexpr double PHEROMONE_DEPOSIT      = 1.00;
static constexpr double PHEROMONE_DECAY        = 0.80;
static constexpr double PHEROMONE_SENSE_THRESH = 0.25;
static constexpr double COVERAGE_BIAS          = 0.25;

// Cascade edge weights (Steps 3 & 4)
static constexpr double EDGE_REINFORCE    = 0.50;
static constexpr double EDGE_DECAY        = 0.85;
static constexpr double EDGE_PRUNE_THRESH = 0.10;

// Phase 2 perimeter check
static constexpr int    PERIM_SAMPLES = 36;

// Greedy-MSC
static constexpr double WAKEUP_BUFFER      = 30.0;  // wake up this many seconds before estimated death
static constexpr double LOW_BATTERY_THRESH = 0.05;   // 5% energy triggers warning
static constexpr double GROUP_COLLECT_TIME = 0.5;    // seconds to collect score broadcasts


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
        useGreedyMSC = par("useGreedyMSC");

        posX = uniform(0, areaSize);
        posY = uniform(0, areaSize);

        state           = UNDECIDED;
        hasPred         = false;
        cascadeParentId = -1;

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
                energy = getCurrentEnergy();
                phaseStartTime   = simTime();
                phaseStartEnergy = energy;
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

        // ── Periodic coverage report ───────────────────────────────────────
        case KIND_STATS: {
            delete msg; statsTimer = nullptr;
            cModule *parent = getParentModule();
            int n = parent->par("numSensors");
            int minAlive = -1;
            for (int i = 0; i < n; i++) {
                SensorNode *s = check_and_cast<SensorNode*>(
                                    parent->getSubmodule("sensor", i));
                if (s->getCurrentEnergy() > 0.0) { minAlive = i; break; }
            }
            if (nodeId == minAlive)
                computeCoverage();
            lastReportedRound = roundCount;

            // Reschedule
            double statsInterval = par("statsInterval").doubleValue();
            simtime_t nextFire = simTime() + statsInterval;
            simtime_t roundEnd = roundStartTime + roundTime;
            if (nextFire < roundEnd) {
                statsTimer = new cMessage("stats", KIND_STATS);
                scheduleAt(nextFire, statsTimer);
            }
            break;
        }

        // ── Greedy-MSC Phase 3 step 1: broadcast score ────────────────────
        // Every sleep node with known active neighbours broadcasts its
        // contribution score and preferred active node (most critical).
        case KIND_GROUP_BROADCAST: {
            delete msg; groupBroadcastTimer = nullptr;
            if (!useGreedyMSC || state != OFF) break;
            if (knownActiveInfo.empty()) break;

            // Pick the most critical active neighbour (lowest energy)
            int    bestActiveId = -1;
            double lowestEnergy = 1e9;
            for (auto &info : knownActiveInfo) {
                if (info.energy < lowestEnergy) {
                    lowestEnergy = info.energy;
                    bestActiveId = info.id;
                }
            }
            if (bestActiveId < 0) break;

            assignedActiveNodeId = bestActiveId;

            // Compute contribution score: high pheromone (covers many
            // active areas) × high energy (will last long as replacement)
            double phero = localCoverage(posX, posY);
            double eFrac = std::max(0.01, energy / initEnergy);
            myGroupScore = phero * eFrac;

            // Add ourselves to the collected scores
            collectedScores.push_back({nodeId, myGroupScore, assignedActiveNodeId});

            // Broadcast to neighbours
            broadcastGroupScore();

            EV_INFO << "[Node " << nodeId << "] broadcast group score="
                    << myGroupScore << " for active " << bestActiveId << "\n";
            break;
        }

        // ── Greedy-MSC Phase 3 step 2: compute group assignment ───────────
        // After collecting all neighbours' scores, each node ranks itself
        // among all sleep nodes that target the same active node.
        // Rank 1 → Group 1 (first replacement), Rank 2 → Group 2, etc.
        case KIND_GROUPING: {
            delete msg; groupingTimer = nullptr;
            if (!useGreedyMSC || state != OFF) break;
            if (assignedActiveNodeId < 0) break;

            computeGroupAssignment();
            break;
        }

        // ── Active node low-battery warning ────────────────────────────────
        case KIND_LOW_BATTERY: {
            delete msg; lowBatteryTimer = nullptr;
            if (state != ON) break;

            double drainRate   = 4.0 / 100.0;
            double currentE    = getCurrentEnergy();
            double timeToZero  = currentE / drainRate;
            double deathTime   = simTime().dbl() + timeToZero;

            int coveredId = (assignedActiveNodeId >= 0)
                            ? assignedActiveNodeId : nodeId;

            EV_INFO << "[Node " << nodeId << "] LOW BATTERY warning"
                    << " (covering active " << coveredId
                    << ", death≈" << deathTime << "s)\n";
            broadcastLowBatteryWarning(coveredId, deathTime);
            break;
        }

        // ── Replacement wake-up ────────────────────────────────────────────
        case KIND_WAKEUP: {
            delete msg; wakeupTimer = nullptr;
            if (state != OFF) break;

            EV_INFO << "[Node " << nodeId << "] waking up as Group "
                    << myGroupNumber << " replacement for active node "
                    << assignedActiveNodeId << "\n";
            turnOnAsReplacement();
            break;
        }

        case KIND_DISPLAY_REFRESH:
            delete msg; refreshTimer = nullptr;
            updateDisplay();
            if (simTime() < (roundCount * roundTime) - 1.0) {
                refreshTimer = new cMessage("refresh", KIND_DISPLAY_REFRESH);
                scheduleAt(simTime() + 50.0, refreshTimer);
            }
            break;

        default: delete msg;
        }

    } else {
        int kind = msg->getKind();
        if      (kind == 100) handlePowerOn(check_and_cast<PowerOnMsg*>(msg));
        else if (kind == 101) handleCoverageMark(check_and_cast<CoverageMarkMsg*>(msg));
        else if (kind == 102) handleGroupScore(check_and_cast<GroupScoreMsg*>(msg));
        else if (kind == 103) handleLowBatteryWarning(check_and_cast<LowBatteryWarningMsg*>(msg));
        else                  delete msg;
    }
}

// ════════════════════════════════ Round Management ═══════════════════════════

void SensorNode::startRound()
{
    roundStartTime = simTime();

    cancelTimer(volunteerTimer);
    cancelTimer(selectTimer);
    cancelTimer(steadyTimer);
    cancelTimer(statsTimer);
    cancelTimer(pruneTimer);
    cancelTimer(redundancyTimer);
    cancelTimer(refreshTimer);
    cancelTimer(groupBroadcastTimer);
    cancelTimer(groupingTimer);
    cancelTimer(lowBatteryTimer);
    cancelTimer(wakeupTimer);

    state           = UNDECIDED;
    hasPred         = false;
    bestTime        = 1e9;
    bestHasPred     = false;
    bestPredId      = -1;
    cascadeParentId = -1;

    knownActive.clear();

    // Reset Greedy-MSC state
    if (useGreedyMSC) {
        myGroupNumber        = -1;
        assignedActiveNodeId = -1;
        myGroupScore         = 0;
        collectedScores.clear();
        warningsReceived.clear();
        knownActiveInfo.clear();
    }

    updateDisplay();

    // Volunteer timer: energy-weighted
    double lostFrac = 1.0 - energy / initEnergy;
    lostFrac = std::max(0.0, std::min(lostFrac, 0.9));
    double vt = exponential(0.05 + lostFrac * 0.40);
    vt = std::min(vt, Ts * 0.9);
    volunteerTimer = new cMessage("vol", KIND_VOLUNTEER);
    scheduleAt(simTime() + vt, volunteerTimer);

    // Steady timer: end of Phase 1
    steadyTimer = new cMessage("steady", KIND_STEADY);
    scheduleAt(simTime() + Ts, steadyTimer);

    // Pheromone decay at mid-round
    pruneTimer = new cMessage("prune", KIND_PRUNE);
    scheduleAt(simTime() + roundTime * 0.5, pruneTimer);

    // Phase 2 redundancy check
    redundancyTimer = new cMessage("redcheck", KIND_REDUNDANCY_CHECK);
    scheduleAt(simTime() + Ts + Te + uniform(0, Te), redundancyTimer);

    // Greedy-MSC Phase 3: two-step grouping
    // Step 1 at Ts+3Te: broadcast scores
    // Step 2 at Ts+3Te+0.5s: compute group assignment
    // Stats report after grouping at Ts+3Te+1s
    if (useGreedyMSC) {
        simtime_t groupPhaseStart = simTime() + Ts + 3.0 * Te;
        groupBroadcastTimer = new cMessage("grpbcast", KIND_GROUP_BROADCAST);
        scheduleAt(groupPhaseStart, groupBroadcastTimer);

        groupingTimer = new cMessage("grouping", KIND_GROUPING);
        scheduleAt(groupPhaseStart + GROUP_COLLECT_TIME, groupingTimer);

        statsTimer = new cMessage("stats", KIND_STATS);
        scheduleAt(groupPhaseStart + GROUP_COLLECT_TIME + 0.1, statsTimer);
    } else {
        statsTimer = new cMessage("stats", KIND_STATS);
        scheduleAt(simTime() + Ts + 3.0 * Te, statsTimer);
    }

    // Refresh display
    refreshTimer = new cMessage("refresh", KIND_DISPLAY_REFRESH);
    scheduleAt(simTime() + Ts + 50.0, refreshTimer);
}


// ════════════════════════════════ OGDC + Stigmergy Phase 1 ═══════════════════

void SensorNode::handlePowerOn(PowerOnMsg *m)
{
    if (state != UNDECIDED) { delete m; return; }

    markCoverage(m->senderX, m->senderY);
    if (m->hasPred && m->predX >= 0)
        markCoverage(m->predX, m->predY);

    cancelTimer(volunteerTimer);

    double selectTime;
    double myPredX   = m->senderX;
    double myPredY   = m->senderY;
    bool   myHasPred = true;

    if (!m->hasPred) {
        double dA  = dist(posX, posY, m->senderX, m->senderY);
        if (dA < 1e-6) { delete m; return; }
        double dev = std::fabs(dA - std::sqrt(3.0) * rs);
        selectTime = Td + (dev / rs) * Ts;
    } else {
        double dT1 = (m->t1x >= 0) ? dist(posX, posY, m->t1x, m->t1y) : 1e9;
        double dT2 = (m->t2x >= 0) ? dist(posX, posY, m->t2x, m->t2y) : 1e9;
        double best = std::min(dT1, dT2);
        if (best > rs) { delete m; return; }
        selectTime = Td + (best / rs) * Ts;
    }

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
    phaseStartTime   = simTime();
    phaseStartEnergy = energy;
    cancelTimer(volunteerTimer);
    cancelTimer(selectTimer);
    energy -= 20.0 * Td / 100.0;
    energy  = std::max(0.0, energy);
    updateDisplay();

    markCoverage(posX, posY);

    if (cascadeParentId >= 0) {
        edgeWeight[cascadeParentId] += EDGE_REINFORCE;
        EV_INFO << "[Node " << nodeId << "] reinforced edge from node "
                << cascadeParentId
                << " weight=" << edgeWeight[cascadeParentId] << "\n";
    }

    broadcastPowerOn();
    broadcastCoverageMark();

    // Schedule low-battery warning
    if (useGreedyMSC) {
        double fivePercent       = LOW_BATTERY_THRESH * initEnergy;
        double drainRate         = 4.0 / 100.0;
        double timeToWarning     = (phaseStartEnergy - fivePercent) / drainRate;
        if (timeToWarning > 0) {
            cancelTimer(lowBatteryTimer);
            lowBatteryTimer = new cMessage("lowbat", KIND_LOW_BATTERY);
            scheduleAt(simTime() + timeToWarning, lowBatteryTimer);
        }
    }
}

void SensorNode::turnOff()
{
    state = OFF;
    phaseStartTime   = simTime();
    phaseStartEnergy = energy;
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
    tmpl.senderId        = nodeId;
    tmpl.senderX         = posX;
    tmpl.senderY         = posY;
    tmpl.sensingRange    = rs;
    tmpl.remainingEnergy = getCurrentEnergy();

    for (auto &kv : neighbours) {
        cModule *mod = getParentModule()->getSubmodule("sensor", kv.first);
        sendDirect(tmpl.dup(), mod, "radioIn");
    }
}

void SensorNode::handleCoverageMark(CoverageMarkMsg *m)
{
    markCoverage(m->senderX, m->senderY, PHEROMONE_DEPOSIT);

    bool found = false;
    for (auto &k : knownActive) {
        if (std::fabs(k.first - m->senderX) < 1e-6 &&
            std::fabs(k.second - m->senderY) < 1e-6) { found = true; break; }
    }
    if (!found)
        knownActive.push_back({m->senderX, m->senderY});

    if (useGreedyMSC) {
        bool infoFound = false;
        for (auto &info : knownActiveInfo) {
            if (info.id == m->senderId) {
                info.energy = m->remainingEnergy;
                infoFound = true;
                break;
            }
        }
        if (!infoFound)
            knownActiveInfo.push_back(
                {m->senderId, m->senderX, m->senderY, m->remainingEnergy});
    }

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

void SensorNode::decayAndPruneEdges()
{
    for (int gx = 0; gx < gridRes; gx++)
        for (int gy = 0; gy < gridRes; gy++)
            stigGrid[gx][gy] *= PHEROMONE_DECAY;

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

bool SensorNode::isPerimeterRedundant() const
{
    cModule *parent = getParentModule();
    int n = parent->par("numSensors");

    std::vector<std::pair<double,double>> liveActive;
    for (int i = 0; i < n; i++) {
        if (i == nodeId) continue;
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

        px = std::max(0.0, std::min(areaSize, px));
        py = std::max(0.0, std::min(areaSize, py));

        bool covered = false;
        for (auto &a : liveActive) {
            if (dist(px, py, a.first, a.second) <= rs) {
                covered = true;
                break;
            }
        }
        if (!covered) return false;
    }
    return true;
}

// ════════════════════════════════ Greedy-MSC Phase 3 ═════════════════════════

// ── Step 1: broadcast score ───────────────────────────────────────────────
void SensorNode::broadcastGroupScore()
{
    GroupScoreMsg tmpl;
    tmpl.senderId         = nodeId;
    tmpl.score            = myGroupScore;
    tmpl.assignedActiveId = assignedActiveNodeId;

    for (auto &kv : neighbours) {
        cModule *mod = getParentModule()->getSubmodule("sensor", kv.first);
        sendDirect(tmpl.dup(), mod, "radioIn");
    }
}

// ── Handle incoming score from another sleep node ─────────────────────────
void SensorNode::handleGroupScore(GroupScoreMsg *m)
{
    if (useGreedyMSC)
        collectedScores.push_back({m->senderId, m->score, m->assignedActiveId});
    delete m;
}

// ── Step 2: deterministic group assignment ────────────────────────────────
// Among all sleep nodes that target the same active node, rank by score
// (descending).  Rank 1 = Group 1, rank 2 = Group 2, etc.
// Ties broken by node ID (lower ID wins higher rank).
void SensorNode::computeGroupAssignment()
{
    if (assignedActiveNodeId < 0) return;

    // Collect all competitors for the same active node (including self)
    std::vector<std::pair<double, int>> competitors;  // {score, nodeId}
    for (auto &sc : collectedScores) {
        if (sc.assignedActiveId == assignedActiveNodeId)
            competitors.push_back({sc.score, sc.nodeId});
    }

    // Sort: highest score first; tie-break by lowest nodeId
    std::sort(competitors.begin(), competitors.end(),
        [](const std::pair<double,int> &a, const std::pair<double,int> &b) {
            if (std::fabs(a.first - b.first) > 1e-9)
                return a.first > b.first;   // higher score = better rank
            return a.second < b.second;     // lower ID wins tie
        });

    // Find my rank
    int rank = 1;
    for (auto &c : competitors) {
        if (c.second == nodeId) break;
        rank++;
    }

    myGroupNumber = rank;

    EV_INFO << "[Node " << nodeId << "] assigned Group " << myGroupNumber
            << " for active node " << assignedActiveNodeId
            << " (score=" << myGroupScore
            << ", competitors=" << competitors.size() << ")\n";
}

// ── Low-battery warning broadcast ─────────────────────────────────────────
void SensorNode::broadcastLowBatteryWarning(int coveredId, double deathTime)
{
    LowBatteryWarningMsg tmpl;
    tmpl.senderId           = nodeId;
    tmpl.senderX            = posX;
    tmpl.senderY            = posY;
    tmpl.coveredActiveId    = coveredId;
    tmpl.estimatedDeathTime = deathTime;

    for (auto &kv : neighbours) {
        cModule *mod = getParentModule()->getSubmodule("sensor", kv.first);
        sendDirect(tmpl.dup(), mod, "radioIn");
    }
}

// ── Handle incoming low-battery warning ───────────────────────────────────
// Count warnings for each original active node.  When the count matches
// this node's group number, it is next in the chain and schedules wake-up.
void SensorNode::handleLowBatteryWarning(LowBatteryWarningMsg *m)
{
    if (!useGreedyMSC || state != OFF || myGroupNumber < 1) { delete m; return; }

    warningsReceived[m->coveredActiveId]++;
    int warnCount = warningsReceived[m->coveredActiveId];

    if (assignedActiveNodeId == m->coveredActiveId &&
        myGroupNumber == warnCount)
    {
        double wakeupTime = m->estimatedDeathTime - WAKEUP_BUFFER;
        if (wakeupTime <= simTime().dbl())
            wakeupTime = simTime().dbl() + 0.1;

        cancelTimer(wakeupTimer);
        wakeupTimer = new cMessage("wakeup", KIND_WAKEUP);
        scheduleAt(wakeupTime, wakeupTimer);

        EV_INFO << "[Node " << nodeId << "] Group " << myGroupNumber
                << " scheduling wake-up at t=" << wakeupTime
                << " for active node " << assignedActiveNodeId << "\n";
    }

    delete m;
}

// ── Turn ON as a scheduled replacement ────────────────────────────────────
void SensorNode::turnOnAsReplacement()
{
    state = ON;
    phaseStartTime   = simTime();
    phaseStartEnergy = energy;
    energy -= 20.0 * Td / 100.0;
    energy  = std::max(0.0, energy);
    updateDisplay();

    markCoverage(posX, posY);
    broadcastCoverageMark();

    // Schedule our own low-battery warning for the next group
    if (useGreedyMSC) {
        double fivePercent       = LOW_BATTERY_THRESH * initEnergy;
        double drainRate         = 4.0 / 100.0;
        double timeToWarning     = (phaseStartEnergy - fivePercent) / drainRate;
        if (timeToWarning > 0) {
            cancelTimer(lowBatteryTimer);
            lowBatteryTimer = new cMessage("lowbat", KIND_LOW_BATTERY);
            scheduleAt(simTime() + timeToWarning, lowBatteryTimer);
        }
    }
}

// ════════════════════════════════ Coverage & Display ═════════════════════════

double SensorNode::getCurrentEnergy() const
{
    double rate = (state == ON) ? 4.0 / 100.0 : 0.01 / 100.0;
    double elapsed = (simTime() - phaseStartTime).dbl();
    return std::max(0.0, phaseStartEnergy - elapsed * rate);
}

void SensorNode::computeCoverage()
{
    cModule *parent = getParentModule();
    int n = parent->par("numSensors");

    std::vector<std::pair<double,double>> active;
    int onCount = 0, offCount = 0, deadCount = 0;
    double sumEnergy = 0.0, maxEnergy = 0.0, minEnergy = 1e9;
    int    maxNode = -1, minNode = -1, avgNode = -1;
    double avgDiff = 1e9;

    std::map<int,int> groupCounts;
    int ungroupedSleep = 0;

    for (int i = 0; i < n; i++) {
        SensorNode *s = check_and_cast<SensorNode*>(
                            parent->getSubmodule("sensor", i));
        double liveE = s->getCurrentEnergy();
        sumEnergy += liveE;
        if (liveE > maxEnergy) { maxEnergy = liveE; maxNode = i; }
        if (liveE < minEnergy) { minEnergy = liveE; minNode = i; }
        if (s->getCurrentEnergy() <= LOW_BATTERY_THRESH * s->initEnergy) {
            deadCount++;
        } else if (s->state == ON) {
            active.push_back({s->posX, s->posY});
            onCount++;
        } else if (s->state == OFF) {
            offCount++;
        }

        if (useGreedyMSC && s->state == OFF &&
            s->getCurrentEnergy() > LOW_BATTERY_THRESH * s->initEnergy)
        {
            if (s->myGroupNumber > 0)
                groupCounts[s->myGroupNumber]++;
            else
                ungroupedSleep++;
        }
    }

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

    double avgEnergy = sumEnergy / n;
    for (int i = 0; i < n; i++) {
        SensorNode *s = check_and_cast<SensorNode*>(
                            parent->getSubmodule("sensor", i));
        double diff = std::fabs(s->getCurrentEnergy() - avgEnergy);
        if (diff < avgDiff) { avgDiff = diff; avgNode = i; }
    }

    double cov = (double)covered / ((long long)G * G);

    EV_INFO << "\n+---------- ROUND " << roundCount << " REPORT ----------+\n"
            << " Coordinator  : Node " << nodeId         << "\n"
            << " Active (ON)  : " << onCount             << "\n"
            << " Sleeping(OFF): " << offCount            << "\n"
            << " Dead/Low bat : " << deadCount           << "\n"
            << " Grid         : " << G << "x" << G
                                  << " (" << cs << "m)\n"
            << " Avg energy   : " << avgEnergy << " (Node " << avgNode << ")\n"
            << " Max energy   : " << maxEnergy << " (Node " << maxNode << ")\n"
            << " Min energy   : " << minEnergy << " (Node " << minNode << ")\n"
            << " Coverage     : " << cov * 100.0         << "%\n"
            << " Sim time     : " << simTime()           << "s\n";
    if (useGreedyMSC) {
        EV_INFO << " Greedy-MSC   : ON\n";
        for (auto &gc : groupCounts)
            EV_INFO << "   Group " << gc.first << "    : " << gc.second << " nodes\n";
        EV_INFO << "   Ungrouped  : " << ungroupedSleep << " nodes\n";
    }
    EV_INFO << "+--------------------------------------+\n\n";

    std::cout << "\n+---------- ROUND " << roundCount << " REPORT ----------+\n" << std::endl;
    std::cout << " Round        : " << roundCount             << std::endl;
    std::cout << " Coordinator  : Node " << nodeId            << std::endl;
    std::cout << " Active (ON)  : " << onCount                << std::endl;
    std::cout << " Sleeping(OFF): " << offCount               << std::endl;
    std::cout << " Dead/Low bat : " << deadCount              << std::endl;
    std::cout << " Grid         : " << G << "x" << G
              << " (" << cs << "m)"                           << std::endl;
    std::cout << " Avg energy   : " << avgEnergy << " (Node " << avgNode << ")" << std::endl;
    std::cout << " Max energy   : " << maxEnergy << " (Node " << maxNode << ")" << std::endl;
    std::cout << " Min energy   : " << minEnergy << " (Node " << minNode << ")" << std::endl;
    std::cout << " Coverage     : " << cov * 100.0 << "%"     << std::endl;
    std::cout << " Sim time     : " << simTime() << "s"       << std::endl;
    if (useGreedyMSC) {
        std::cout << " Greedy-MSC   : ON"                     << std::endl;
        for (auto &gc : groupCounts)
            std::cout << "   Group " << gc.first << "    : " << gc.second << " nodes" << std::endl;
        std::cout << "   Ungrouped  : " << ungroupedSleep << " nodes" << std::endl;
    }
    std::cout << "+--------------------------------------+\n\n" << std::endl;
    std::cout.flush();

    emit(sigCoverage, cov);
    emit(sigActive,   (double)onCount);
}

void SensorNode::updateDisplay()
{
    if (!hasGUI()) return;
    const double SCALE = 10.0;
    char buf[256];
    if (getCurrentEnergy() <= LOW_BATTERY_THRESH * initEnergy) {
        std::snprintf(buf, sizeof(buf),
            "p=%.0f,%.0f;b=10,10,oval,#000000,#000000,1",
            posX * SCALE, posY * SCALE);
    } else if (state == ON) {
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
    cModule *par = getParentModule();
    int n = par->par("numSensors");
    int minAlive = -1;
    for (int i = 0; i < n; i++) {
        SensorNode *s = check_and_cast<SensorNode*>(
                            par->getSubmodule("sensor", i));
        if (s->getCurrentEnergy() > 0.0) { minAlive = i; break; }
    }
    if (nodeId == minAlive) {
        std::cout << "\n>>> SIMULATION FINISHED after "
                  << roundCount << " rounds <<<" << std::endl;
        computeCoverage();
    }
}
