#include "SensorNode.h"
#include <algorithm>
#include <cstdio>
#include <iostream>

Define_Module(SensorNode);

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

        state   = UNDECIDED;
        hasPred = false;

        sigCoverage = registerSignal("coverage");
        sigActive   = registerSignal("activeNodes");

        updateDisplay();

    } else {
        cModule *parent = getParentModule();
        for (int i = 0; i < numNodes; i++) {
            if (i == nodeId) continue;
            SensorNode *s = check_and_cast<SensorNode*>(
                                parent->getSubmodule("sensor", i));
            if (dist(posX,posY, s->posX,s->posY) <= rc)
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

            // Deduct energy from the PREVIOUS round's steady-state phase
            if (roundCount > 0) {
                // Steady phase duration = roundTime - Ts
                double dt = roundTime - Ts;
                // Paper energy ratio: ON(idle/receive)=4, OFF(sleep)=0.01
                // Scale down by 100 so nodes last many rounds
                if (state == ON)
                    energy -= dt * 4.0 / 100.0;
                else
                    energy -= dt * 0.01 / 100.0;
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
                hasPred = false;
                turnOn();
            }
            break;

        // ── Selected by cascade ────────────────────────────────────────────
        case KIND_SELECT:
            delete msg; selectTimer = nullptr;
            if (state == UNDECIDED) {
                EV_INFO << "[Node " << nodeId << "] selected by cascade\n";
                hasPred  = bestHasPred;
                predX    = bestPredX;
                predY    = bestPredY;
                turnOn();
            }
            break;

        // ── End of selection phase: remaining UNDECIDED → OFF ─────────────
        case KIND_STEADY:
            delete msg; steadyTimer = nullptr;
            if (state == UNDECIDED) turnOff();
            break;

        // ── Measure coverage AFTER selection phase settles ─────────────────
        case KIND_STATS:
            delete msg; statsTimer = nullptr;
            if (nodeId == 0) computeCoverage();
            break;

        default:
            delete msg;
        }

    } else {
        if (msg->getKind() == 100)
            handlePowerOn(check_and_cast<PowerOnMsg*>(msg));
        else
            delete msg;
    }
}

// ════════════════════════════════ Round Management ═══════════════════════════

void SensorNode::startRound()
{
    cancelTimer(volunteerTimer);
    cancelTimer(selectTimer);
    cancelTimer(steadyTimer);
    cancelTimer(statsTimer);

    state       = UNDECIDED;
    hasPred     = false;
    bestTime    = 1e9;
    bestHasPred = false;
    knownActive.clear();   // reset neighbourhood knowledge for the new round
    updateDisplay();

    // -- Volunteer timer ----------------------------------------------------
    // Paper: probabilistic, related to remaining energy
    // Use exponential distribution so statistically only 1-2 nodes
    // volunteer before the first POWER_ON message suppresses the rest
    // Mean = 0.1s means ~1-2 nodes fire before 0.3s
    // Higher energy ? lower mean ? volunteers sooner
    double lostFrac = 1.0 - energy / initEnergy;
    lostFrac = std::max(0.0, std::min(lostFrac, 0.9));

    // Mean volunteer time scales from 0.05s (full energy) to 0.45s (low energy)
    // All well below Ts=1s so at least one node always volunteers
    double meanVt = 0.05 + lostFrac * 0.40;
    double vt     = exponential(meanVt);
    vt = std::min(vt, Ts * 0.9);   // hard cap: always before steady timer

    volunteerTimer = new cMessage("vol", KIND_VOLUNTEER);
    scheduleAt(simTime() + vt, volunteerTimer);

    // Steady timer: ends selection phase, remaining UNDECIDED ? OFF
    steadyTimer = new cMessage("steady", KIND_STEADY);
    scheduleAt(simTime() + Ts, steadyTimer);

    // Stats timer: measure coverage after cascade settles
    if (nodeId == 0) {
        statsTimer = new cMessage("stats", KIND_STATS);
        scheduleAt(simTime() + Ts + 0.5, statsTimer);
    }
}

// ════════════════════════════════ OGDC Core ══════════════════════════════════

void SensorNode::handlePowerOn(PowerOnMsg *m)
{
    if (state != UNDECIDED) { delete m; return; }

    // ── Update our knowledge of which nodes are active ─────────────────────
    // Record the sender (and its predecessor if known) so isRedundant() has
    // an accurate picture of already-covered area around us.
    auto addKnown = [&](double x, double y) {
        for (auto &k : knownActive)
            if (std::fabs(k.first-x) < 1e-6 && std::fabs(k.second-y) < 1e-6) return;
        knownActive.push_back({x, y});
    };
    addKnown(m->senderX, m->senderY);
    if (m->hasPred && m->predX >= 0) addKnown(m->predX, m->predY);

    // ── Redundancy check ───────────────────────────────────────────────────
    // If every point on my sensing perimeter is already covered by a known
    // active node, activating me adds zero new area → skip.
    if (isRedundant()) { delete m; return; }

    // This node is now a CASCADE CANDIDATE, not a starting node
    cancelTimer(volunteerTimer);

    double selectTime;
    double myPredX   = m->senderX;
    double myPredY   = m->senderY;
    bool   myHasPred = true;

    if (!m->hasPred) {
        // ── Selecting node B (first neighbour of starting node A) ──────────
        // Theorem 1: optimal B is at distance √3·rs from A
        // Timer is shorter for nodes closer to the ideal distance
        double dA  = dist(posX, posY, m->senderX, m->senderY);
        if (dA < 1e-6) { delete m; return; }
        double dev = std::fabs(dA - std::sqrt(3.0) * rs);
        selectTime = Td + (dev / rs) * Ts;
    } else {
        // ── Selecting C, D, … (cascade) ───────────────────────────────────
        // Theorem 2: next node closest to one of the two target positions
        double dT1 = (m->t1x >= 0) ? dist(posX,posY, m->t1x, m->t1y) : 1e9;
        double dT2 = (m->t2x >= 0) ? dist(posX,posY, m->t2x, m->t2y) : 1e9;
        double best = std::min(dT1, dT2);
        if (best > rs) { delete m; return; }   // too far from any target
        selectTime = Td + (best / rs) * Ts;
    }

    // Accept only if this gives a shorter timer than previous messages
    if (selectTime < bestTime) {
        bestTime    = selectTime;
        bestHasPred = myHasPred;
        bestPredX   = myPredX;
        bestPredY   = myPredY;

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
    // Transmit cost: paper ratio 20, for one Td period
    energy -= 20.0 * Td / 100.0;
    energy  = std::max(0.0, energy);
    updateDisplay();
    broadcastPowerOn();
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
    tmpl.senderX  = posX;    tmpl.senderY  = posY;
    tmpl.hasPred  = hasPred;
    tmpl.predX    = predX;   tmpl.predY    = predY;
    tmpl.t1x = t1x; tmpl.t1y = t1y;
    tmpl.t2x = t2x; tmpl.t2y = t2y;

    for (auto &kv : neighbours) {
        cModule *mod = getParentModule()->getSubmodule("sensor", kv.first);
        sendDirect(tmpl.dup(), mod, "radioIn");
    }
}

// ════════════════════════════════ Geometry ═══════════════════════════════════
//
//  Given two active nodes A and B, compute the two optimal next
//  positions T1, T2 for the next node in the cascade.
//  Derived from Theorems 1 & 2 of the paper.
//
void SensorNode::computeTargets(double ax, double ay,
                                 double bx, double by,
                                 double &t1x, double &t1y,
                                 double &t2x, double &t2y)
{
    t1x = t1y = t2x = t2y = -1;

    double dAB = dist(ax,ay, bx,by);
    if (dAB < 1e-6 || dAB > 2.0*rs) return;   // disks don't intersect

    double mx = (ax+bx)*0.5,  my = (ay+by)*0.5;  // midpoint M
    double ux = (bx-ax)/dAB,  uy = (by-ay)/dAB;  // unit vector A→B
    double px = -uy,           py = ux;            // perpendicular

    // Distance from M to each crossing point
    double h = std::sqrt(rs*rs - (dAB*0.5)*(dAB*0.5));

    double o1x = mx + h*px,  o1y = my + h*py;    // crossing point 1
    double o2x = mx - h*px,  o2y = my - h*py;    // crossing point 2

    // Target = crossing_point pushed outward by rs
    auto makeTarget = [&](double ox, double oy,
                          double &tx, double &ty) {
        double dx = ox - mx, dy = oy - my;
        double len = std::sqrt(dx*dx + dy*dy);
        if (len < 1e-9) { tx = ty = -1; return; }
        tx = ox + (dx/len) * rs;
        ty = oy + (dy/len) * rs;
    };
    makeTarget(o1x, o1y, t1x, t1y);
    makeTarget(o2x, o2y, t2x, t2y);
}

// ════════════════════════════════ Coverage & Display ═════════════════════════

// Perimeter coverage check (based on Lemma 1 / CCP concept from the paper).
//
// We sample SAMPLES equally-spaced points on the boundary of this node's
// sensing disk.  If EVERY sample point is within rs of at least one known
// active neighbour, then this node's entire sensing area is already covered
// by neighbours — activating it would add no new coverage (redundant).
//
// 36 samples (one every 10°) gives <1% geometric error and is fast to compute.
bool SensorNode::isRedundant() const
{
    if (knownActive.empty()) return false;   // no active neighbours known → not redundant

    const int    SAMPLES = 36;
    const double TWO_PI  = 2.0 * M_PI;

    for (int i = 0; i < SAMPLES; i++) {
        double angle = TWO_PI * i / SAMPLES;
        double px = posX + rs * std::cos(angle);
        double py = posY + rs * std::sin(angle);

        // Clamp to deployment area so boundary nodes aren't penalised
        px = std::max(0.0, std::min(areaSize, px));
        py = std::max(0.0, std::min(areaSize, py));

        bool covered = false;
        for (auto &a : knownActive) {
            if (dist(px, py, a.first, a.second) <= rs) {
                covered = true;
                break;
            }
        }
        if (!covered) return false;   // found an uncovered perimeter point
    }
    return true;   // all perimeter points covered → this node is redundant
}

void SensorNode::computeCoverage()
{
    cModule *parent = getParentModule();
    int n = parent->par("numSensors");

    std::vector<std::pair<double,double>> active;
    int onCount  = 0;
    int offCount = 0;

    for (int i = 0; i < n; i++) {
        SensorNode *s = check_and_cast<SensorNode*>(
                            parent->getSubmodule("sensor", i));
        if (s->state == ON) {
            active.push_back({s->posX, s->posY});
            onCount++;
        } else if (s->state == OFF) {
            offCount++;
        }
    }

    // Paper: 50×50 grid, coverage = fraction of 1m² cells covered
    const int G = (int)areaSize;
    int covered = 0;
    for (int gx = 0; gx < G; gx++) {
        for (int gy = 0; gy < G; gy++) {
            double cx = gx + 0.5, cy = gy + 0.5;
            for (auto &p : active) {
                if (dist(p.first, p.second, cx, cy) <= rs) {
                    covered++;
                    break;
                }
            }
        }
    }

    double cov = (double)covered / (G * G);

    // EV_INFO → shows in OMNeT++ console as INFO
    EV_INFO << "\n========================================\n"
            << " Round        : " << roundCount        << "\n"
            << " Active (ON)  : " << onCount           << "\n"
            << " Sleeping(OFF): " << offCount          << "\n"
            << " Coverage     : " << cov * 100.0       << "%\n"
            << " Sim time     : " << simTime()         << "s\n"
            << "========================================\n\n";

    // std::cout → always visible regardless of log level
    std::cout << "\n========================================"  << std::endl;
    std::cout << " Round        : " << roundCount             << std::endl;
    std::cout << " Active (ON)  : " << onCount                << std::endl;
    std::cout << " Sleeping(OFF): " << offCount               << std::endl;
    std::cout << " Coverage     : " << cov * 100.0 << "%"     << std::endl;
    std::cout << " Sim time     : " << simTime() << "s"       << std::endl;
    std::cout << "========================================"    << std::endl;

    emit(sigCoverage, cov);
    emit(sigActive,   (double)onCount);
}

void SensorNode::updateDisplay()
{
    if (!hasGUI()) return;

    // Scale: world coords 0..areaSize -> display pixels 0..areaSize*10
    const double SCALE = 10.0;

    char buf[256];
    if (state == ON) {
        // Green ring sized to the sensing range — makes coverage immediately
        // visible in the GUI.  Diameter in pixels = 2 * rs * SCALE.
        double ringDiam = 2.0 * rs * SCALE;
        std::snprintf(buf, sizeof(buf),
            "p=%.0f,%.0f;b=%.0f,%.0f,oval,,#00cc00,2",
            posX * SCALE, posY * SCALE,
            ringDiam, ringDiam);
    } else if (state == OFF) {
        // Small red dot, no range ring
        std::snprintf(buf, sizeof(buf),
            "p=%.0f,%.0f;b=8,8,oval,#cc0000,black,1",
            posX * SCALE, posY * SCALE);
    } else {
        // UNDECIDED: grey dot
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
