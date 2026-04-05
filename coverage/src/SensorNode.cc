#include "SensorNode.h"
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <cmath>

Define_Module(SensorNode);

// ════════════════════════════════════════════════════════════════════════════
//  Tuning constants
// ════════════════════════════════════════════════════════════════════════════

static constexpr double PHEROMONE_DEPOSIT      = 1.00;
static constexpr double PHEROMONE_DECAY        = 0.80;
static constexpr double PHEROMONE_SENSE_THRESH = 0.25;
static constexpr double COVERAGE_BIAS          = 0.25;

static constexpr double EDGE_REINFORCE    = 0.50;
static constexpr double EDGE_DECAY        = 0.85;
static constexpr double EDGE_PRUNE_THRESH = 0.10;

static constexpr int    PERIM_SAMPLES = 36;

static constexpr double WAKEUP_BUFFER      = 30.0;
static constexpr double LOW_BATTERY_THRESH = 0.05;
static constexpr double RESELECT_COOLDOWN_BASE = 5.0;


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
        useGreedyMSC         = par("useGreedyMSC");
        clusterDistThreshold = par("clusterDistThreshold");
        greedyW              = par("greedyW");
        warmupRounds         = par("warmupRounds");
        warmupRoundTime      = par("warmupRoundTime");
        ewClusterThreshold   = par("ewClusterThreshold");
        useReselection       = par("useReselection");
        coverageThreshold    = par("coverageThreshold");
        maxReselections      = par("maxReselections");

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
            double rt = (useGreedyMSC && roundCount <= warmupRounds)
                        ? warmupRoundTime : roundTime;
            scheduleAt(simTime() + rt, roundTimer);
            break;
        }

        case KIND_VOLUNTEER:
            delete msg; volunteerTimer = nullptr;
            if (state == UNDECIDED) {
                hasPred = false; cascadeParentId = -1;
                turnOn();
            }
            break;

        case KIND_SELECT:
            delete msg; selectTimer = nullptr;
            if (state == UNDECIDED) {
                hasPred         = bestHasPred;
                predX           = bestPredX;
                predY           = bestPredY;
                cascadeParentId = bestPredId;
                turnOn();
            }
            break;

        case KIND_STEADY:
            delete msg; steadyTimer = nullptr;
            if (state == UNDECIDED) turnOff();
            break;

        case KIND_REDUNDANCY_CHECK: {
            delete msg; redundancyTimer = nullptr;
            if (state == ON && isPerimeterRedundant()) {
                EV_INFO << "[Node " << nodeId << "] Phase-2 pruned\n";
                turnOff();
            }
            break;
        }

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
            bool firstFire = (lastReportedRound != roundCount);
            if (nodeId == minAlive)
                computeCoverage(firstFire);
            if (firstFire) lastReportedRound = roundCount;

            double statsInterval = par("statsInterval").doubleValue();
            simtime_t nextFire = simTime() + statsInterval;
            simtime_t roundEnd = roundStartTime + roundTime;
            if (nextFire < roundEnd) {
                statsTimer = new cMessage("stats", KIND_STATS);
                scheduleAt(nextFire, statsTimer);
            }
            break;
        }

        // ════════════════ Greedy-MSC Phase 1: active node clustering ════════

        case KIND_CLUSTER_SEED: {
            delete msg; clusterSeedTimer = nullptr;
            if (!useGreedyMSC || state != ON) break;
            if (myClusterId >= 0) break;   // already joined a cluster via direct invite

            // ── Compute hub score from receivedChildEdgeWeights ──────────────
            // Hub score = sum of edge weights children hold pointing to this node.
            // High hub score → structural backbone → seed a cluster.
            double hubScore = 0.0;
            for (auto &kv : receivedChildEdgeWeights) hubScore += kv.second;

            if (receivedChildEdgeWeights.empty()) {
                // No children reported — not a backbone node this round.
                // Fall back to avg-distance seeding.
                myClusterId = nodeId;
                myAvgActiveDist = 0; int ac = 0;
                for (auto &info : knownActiveInfo) {
                    double d = dist(posX, posY, info.x, info.y);
                    if (d <= rs) { myAvgActiveDist += d; ac++; }
                }
                if (ac > 0) myAvgActiveDist /= ac;
                broadcastClusterInvite();
                EV_INFO << "[Node " << nodeId << "] fallback cluster (no children)"
                        << " avgDist=" << myAvgActiveDist << "\n";
                break;
            }

            // ── EW-based: stagger seed timer by hub score ────────────────────
            // High hub-score nodes fire earlier → their invites arrive before
            // lower-score competitors, so conflict resolution is implicit.
            static constexpr double SEED_WINDOW = 0.4;
            double seedDelay = SEED_WINDOW / (1.0 + hubScore);

            // Re-schedule via a lambda: cancel existing and fire after delay.
            // Since we're already IN the timer handler, just do the work after
            // the delay by rescheduling a new timer.  Use clusterDoneTimer slot
            // as a proxy is wrong — just schedule a new self-message with a
            // small unique kind. Instead, do the invite work inline after delay
            // by posting a new KIND_CLUSTER_SEED — but that risks re-entry.
            // Cleanest: do the work here, having arrived via the staggered timer
            // (the stagger is applied at scheduleAt time in startRound — see below).
            // HERE we just execute: create cluster and send direct invites.

            myClusterId = nodeId;

            // ── Find qualifying children ─────────────────────────────────────
            // Best child EW for this backbone node
            double bestChildEW = 0.0;
            for (auto &kv : receivedChildEdgeWeights)
                if (kv.second > bestChildEW) bestChildEW = kv.second;

            double ewThreshold = ewClusterThreshold * bestChildEW;

            // Collect qualifying children, sorted descending by EW
            std::vector<std::pair<double,int>> qualifying; // (ew, childId)
            for (auto &kv : receivedChildEdgeWeights)
                if (kv.second >= ewThreshold)
                    qualifying.push_back({kv.second, kv.first});
            std::sort(qualifying.begin(), qualifying.end(),
                      [](const std::pair<double,int> &a, const std::pair<double,int> &b){
                          return a.first > b.first;
                      });

            // Send direct invites to qualifying children
            for (auto &qc : qualifying) {
                int childId = qc.second;
                double ew   = qc.first;
                // Only invite if child is currently sleeping (OFF state) and
                // within radio range
                if (neighbours.find(childId) == neighbours.end()) continue;

                ClusterInviteMsg *inv = new ClusterInviteMsg("CLUSTER_INV");
                inv->senderId          = nodeId;
                inv->clusterHeadId     = myClusterId;
                inv->avgDist           = 0;       // not used for direct invites
                inv->directInvite      = true;
                inv->offeredEdgeWeight = ew;
                sendDirect(inv, getParentModule()->getSubmodule("sensor", childId), "radioIn");
            }

            EV_INFO << "[Node " << nodeId << "] EW-cluster " << myClusterId
                    << " hubScore=" << hubScore
                    << " bestChildEW=" << bestChildEW
                    << " threshold=" << ewThreshold
                    << " invited=" << qualifying.size() << "\n";
            break;
        }

        // Phase 1 ends: active nodes that haven't joined any cluster become
        // self-clusters.  All active nodes then broadcast their cluster ID.
        case KIND_CLUSTER_DONE: {
            delete msg; clusterDoneTimer = nullptr;
            if (!useGreedyMSC) break;

            if (state == ON) {
                if (myClusterId < 0) {
                    myClusterId = nodeId;   // self-cluster
                    EV_INFO << "[Node " << nodeId
                            << "] self-cluster (no group joined)\n";
                }
                broadcastClusterAnnounce();
            }
            break;
        }

        // ════════ Greedy-MSC Phase 2+3: sleep assignment + group ranking ════

        case KIND_GROUPING: {
            delete msg; groupingTimer = nullptr;
            if (!useGreedyMSC || state != OFF) break;
            if (activeNodeClusterMap.empty()) break;

            // ── Phase 2: determine which cluster this sleep node joins ──
            // Count how many active neighbours belong to each cluster
            std::map<int, int> clusterNeighCount;
            for (auto &info : knownActiveInfo) {
                auto it = activeNodeClusterMap.find(info.id);
                if (it != activeNodeClusterMap.end())
                    clusterNeighCount[it->second]++;
            }

            if (clusterNeighCount.empty()) break;

            if (clusterNeighCount.size() == 1) {
                // All active neighbours in one cluster → join it
                myClusterId = clusterNeighCount.begin()->first;
            } else {
                // Multiple clusters → join the one with FEWER active
                // neighbours nearby (that cluster needs more backup)
                int bestCluster = -1;
                int fewestActive = 1e9;
                for (auto &kv : clusterNeighCount) {
                    if (kv.second < fewestActive) {
                        fewestActive = kv.second;
                        bestCluster  = kv.first;
                    }
                }
                myClusterId = bestCluster;
            }

            // ── Phase 3: Greedy-MSC group ranking within cluster ────────
            // Find active nodes in my cluster
            std::vector<ActiveNeighbourInfo> clusterTargets;
            for (auto &info : knownActiveInfo) {
                auto it = activeNodeClusterMap.find(info.id);
                if (it != activeNodeClusterMap.end() && it->second == myClusterId)
                    clusterTargets.push_back(info);
            }
            if (clusterTargets.empty()) break;

            // Pick the nearest active node in my cluster as primary target
            int    targetId  = -1;
            double closestD  = 1e9;
            double targetX   = 0, targetY = 0;
            for (auto &t : clusterTargets) {
                double d = dist(posX, posY, t.x, t.y);
                if (d < closestD) {
                    closestD = d;
                    targetId = t.id;
                    targetX  = t.x;
                    targetY  = t.y;
                }
            }
            assignedActiveNodeId = targetId;

            // Rank among sleep neighbours in the same cluster targeting
            // the same active node.  Closer distance = better rank = Group 1.
            double myDist = closestD;
            int rank = 1;

            for (auto &kv : neighbours) {
                int nid  = kv.first;
                double nx = kv.second.x;
                double ny = kv.second.y;

                // Skip active nodes
                bool isActive = false;
                for (auto &info : knownActiveInfo)
                    if (info.id == nid) { isActive = true; break; }
                if (isActive) continue;

                // Would this sleep neighbour also be in our cluster?
                // Check: is its nearest cluster-active-node also our target?
                // (Approximate: check if our target is the nearest active
                // node in our cluster to this neighbour too)
                double nDistToTarget = dist(nx, ny, targetX, targetY);

                // Check if this neighbour would pick the same target
                bool sameTarget = true;
                for (auto &t : clusterTargets) {
                    if (dist(nx, ny, t.x, t.y) < nDistToTarget - 1e-9) {
                        sameTarget = false;
                        break;
                    }
                }
                if (!sameTarget) continue;

                // Check if also in same cluster (neighbour must also see
                // mainly this cluster's active nodes — approximate check)
                // For simplicity, count it if it's targeting the same node
                if (nDistToTarget < myDist - 1e-9 ||
                    (std::fabs(nDistToTarget - myDist) < 1e-9 && nid < nodeId))
                {
                    rank++;
                }
            }

            myGroupInCluster = rank;

            EV_INFO << "[Node " << nodeId << "] cluster=" << myClusterId
                    << " Group " << myGroupInCluster
                    << " for active " << assignedActiveNodeId
                    << " (dist=" << myDist << ")\n";
            break;
        }

        // ════════════════ Greedy-MSC Phase 4: wake-up chain ═════════════════

        case KIND_LOW_BATTERY: {
            delete msg; lowBatteryTimer = nullptr;
            if (state != ON) break;

            double drainRate  = 4.0 / 100.0;
            double currentE   = getCurrentEnergy();
            double timeToZero = currentE / drainRate;
            double deathTime  = simTime().dbl() + timeToZero;

            int coveredId = (assignedActiveNodeId >= 0)
                            ? assignedActiveNodeId : nodeId;

            EV_INFO << "[Node " << nodeId << "] LOW BATTERY (covering active "
                    << coveredId << " cluster=" << myClusterId
                    << " death≈" << deathTime << "s)\n";
            broadcastLowBatteryWarning(coveredId, myClusterId, deathTime);
            break;
        }

        case KIND_WAKEUP: {
            delete msg; wakeupTimer = nullptr;
            if (state != OFF) break;

            EV_INFO << "[Node " << nodeId << "] waking as Group "
                    << myGroupInCluster << " cluster=" << myClusterId
                    << " for active " << assignedActiveNodeId << "\n";
            turnOnAsReplacement();
            break;
        }

        // ════════════════ Reselection ═══════════════════════════════════════

        case KIND_RESELECT_STEADY:
            delete msg; reselectSteadyTimer = nullptr;
            if (state == UNDECIDED) turnOff();
            break;

        case KIND_RESELECT_REDUND: {
            delete msg; reselectRedundTimer = nullptr;
            if (state == ON && isPerimeterRedundant()) {
                EV_INFO << "[Node " << nodeId << "] Reselection pruned\n";
                turnOff();
            }
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
        // ── Network messages ───────────────────────────────────────────────
        int kind = msg->getKind();
        if      (kind == 100) handlePowerOn(check_and_cast<PowerOnMsg*>(msg));
        else if (kind == 101) handleCoverageMark(check_and_cast<CoverageMarkMsg*>(msg));
        else if (kind == 103) handleLowBatteryWarning(check_and_cast<LowBatteryWarningMsg*>(msg));
        else if (kind == 104) handleReselectTrigger(check_and_cast<ReselectTriggerMsg*>(msg));
        else if (kind == 105) handleClusterInvite(check_and_cast<ClusterInviteMsg*>(msg));
        else if (kind == 106) handleClusterAnnounce(check_and_cast<ClusterAnnounceMsg*>(msg));
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
    cancelTimer(clusterSeedTimer);
    cancelTimer(clusterDoneTimer);
    cancelTimer(groupingTimer);
    cancelTimer(lowBatteryTimer);
    cancelTimer(wakeupTimer);
    cancelTimer(reselectSteadyTimer);
    cancelTimer(reselectRedundTimer);

    state           = UNDECIDED;
    hasPred         = false;
    bestTime        = 1e9;
    bestHasPred     = false;
    bestPredId      = -1;
    cascadeParentId = -1;

    knownActive.clear();

    if (useGreedyMSC) {
        myClusterId         = -1;
        myAvgActiveDist     = 0;
        activeNodeClusterMap.clear();
        myGroupInCluster    = -1;
        assignedActiveNodeId = -1;
        warningsReceived.clear();
        knownActiveInfo.clear();
        receivedChildEdgeWeights.clear();
    }

    if (useReselection) {
        lastReselectTime = -1;
        reselectCount    = 0;
        covAfterReselect = -1;
        reselectGaveUp   = false;
    }

    updateDisplay();

    // Volunteer timer: energy-weighted
    double lostFrac = 1.0 - energy / initEnergy;
    lostFrac = std::max(0.0, std::min(lostFrac, 0.9));
    double vt = exponential(0.05 + lostFrac * 0.40);
    vt = std::min(vt, Ts * 0.9);
    volunteerTimer = new cMessage("vol", KIND_VOLUNTEER);
    scheduleAt(simTime() + vt, volunteerTimer);

    // Steady timer
    steadyTimer = new cMessage("steady", KIND_STEADY);
    scheduleAt(simTime() + Ts, steadyTimer);

    // Pheromone + edge weight decay at mid-round.
    // Runs during warmup rounds too — this is what creates differentiated
    // edge weights before Greedy-MSC clustering activates.
    double currentRoundTime = (useGreedyMSC && roundCount <= warmupRounds)
                              ? warmupRoundTime : roundTime;
    pruneTimer = new cMessage("prune", KIND_PRUNE);
    scheduleAt(simTime() + currentRoundTime * 0.5, pruneTimer);

    // Phase 2 redundancy check
    redundancyTimer = new cMessage("redcheck", KIND_REDUNDANCY_CHECK);
    scheduleAt(simTime() + Ts + Te + uniform(0, Te), redundancyTimer);

    // Greedy-MSC 4-phase scheduling — only after warmup rounds complete.
    if (useGreedyMSC && roundCount > warmupRounds) {
        simtime_t p1Start = simTime() + Ts + 3.0 * Te;

        // Stagger the seed timer by hub score so high-hub backbone nodes
        // fire first. Hub score is sum of child EWs pointing to this node,
        // accumulated in receivedChildEdgeWeights during the selection phase.
        // At startRound time the map is already cleared, so we use the
        // edgeWeight map inversely: approximate hub score as max own EW
        // (a node with high outgoing weight is also likely a backbone recipient).
        // The full hub score is recomputed in KIND_CLUSTER_SEED from
        // receivedChildEdgeWeights populated during the new round's cascade.
        static constexpr double SEED_WINDOW = 0.4;
        double approxHub = 0.0;
        for (auto &kv : edgeWeight) approxHub += kv.second;
        double seedDelay = SEED_WINDOW / (1.0 + approxHub);

        clusterSeedTimer = new cMessage("cluster_seed", KIND_CLUSTER_SEED);
        scheduleAt(p1Start + seedDelay, clusterSeedTimer);

        clusterDoneTimer = new cMessage("cluster_done", KIND_CLUSTER_DONE);
        scheduleAt(p1Start + Ts + 0.4, clusterDoneTimer);  // +0.4 covers full stagger window

        groupingTimer = new cMessage("grouping", KIND_GROUPING);
        scheduleAt(p1Start + Ts + 0.9, groupingTimer);

        statsTimer = new cMessage("stats", KIND_STATS);
        scheduleAt(p1Start + Ts + 1.1, statsTimer);
    } else {
        statsTimer = new cMessage("stats", KIND_STATS);
        scheduleAt(simTime() + Ts + 3.0 * Te, statsTimer);
    }

    refreshTimer = new cMessage("refresh", KIND_DISPLAY_REFRESH);
    scheduleAt(simTime() + Ts + 50.0, refreshTimer);
}


// ════════════════════════════════ OGDC Phase 1 ═══════════════════════════════

void SensorNode::handlePowerOn(PowerOnMsg *m)
{
    if (state != UNDECIDED) { delete m; return; }

    markCoverage(m->senderX, m->senderY);
    if (m->hasPred && m->predX >= 0)
        markCoverage(m->predX, m->predY);

    cancelTimer(volunteerTimer);

    double selectTime;
    double myPredX = m->senderX, myPredY = m->senderY;
    bool myHasPred = true;

    if (!m->hasPred) {
        double dA = dist(posX, posY, m->senderX, m->senderY);
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

    if (cascadeParentId >= 0)
        edgeWeight[cascadeParentId] += EDGE_REINFORCE;

    broadcastPowerOn();
    broadcastCoverageMark();

    if (useGreedyMSC) scheduleLowBatteryWarning();
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
    if (hasPred) computeTargets(predX, predY, posX, posY, t1x, t1y, t2x, t2y);

    PowerOnMsg tmpl;
    tmpl.senderId = nodeId;
    tmpl.senderX = posX; tmpl.senderY = posY;
    tmpl.hasPred = hasPred;
    tmpl.predX = predX; tmpl.predY = predY;
    tmpl.predId = (hasPred && cascadeParentId >= 0) ? cascadeParentId : -1;
    tmpl.t1x = t1x; tmpl.t1y = t1y; tmpl.t2x = t2x; tmpl.t2y = t2y;

    for (auto &kv : neighbours)
        sendDirect(tmpl.dup(), getParentModule()->getSubmodule("sensor", kv.first), "radioIn");
}

void SensorNode::broadcastCoverageMark()
{
    // Include the edge weight this node holds for its cascade parent so the
    // parent can accumulate its hub score from received CoverageMarkMsgs.
    double ewToParent = 0.0;
    if (cascadeParentId >= 0) {
        auto it = edgeWeight.find(cascadeParentId);
        if (it != edgeWeight.end()) ewToParent = it->second;
    }

    CoverageMarkMsg tmpl;
    tmpl.senderId                 = nodeId;
    tmpl.senderX                  = posX;
    tmpl.senderY                  = posY;
    tmpl.sensingRange             = rs;
    tmpl.remainingEnergy          = getCurrentEnergy();
    tmpl.senderCascadeParentId    = cascadeParentId;
    tmpl.senderEdgeWeightToParent = ewToParent;

    for (auto &kv : neighbours)
        sendDirect(tmpl.dup(), getParentModule()->getSubmodule("sensor", kv.first), "radioIn");
}

void SensorNode::handleCoverageMark(CoverageMarkMsg *m)
{
    markCoverage(m->senderX, m->senderY, PHEROMONE_DEPOSIT);

    bool found = false;
    for (auto &k : knownActive)
        if (std::fabs(k.first - m->senderX) < 1e-6 &&
            std::fabs(k.second - m->senderY) < 1e-6) { found = true; break; }
    if (!found) knownActive.push_back({m->senderX, m->senderY});

    if (useGreedyMSC) {
        bool infoFound = false;
        for (auto &info : knownActiveInfo)
            if (info.id == m->senderId) {
                info.energy = m->remainingEnergy;
                infoFound = true;
                break;
            }
        if (!infoFound)
            knownActiveInfo.push_back({m->senderId, m->senderX, m->senderY, m->remainingEnergy});

        // Hub score accumulation: if this sender was activated by ME (nodeId),
        // record how strongly it holds that edge — this is my hub score contribution.
        if (m->senderCascadeParentId == nodeId && m->senderEdgeWeightToParent > 1e-9)
            receivedChildEdgeWeights[m->senderId] = m->senderEdgeWeightToParent;
    }
    delete m;
}

// ════════════════════════════════ Geometry ═══════════════════════════════════

void SensorNode::computeTargets(double ax, double ay, double bx, double by,
                                 double &t1x, double &t1y, double &t2x, double &t2y)
{
    t1x = t1y = t2x = t2y = -1;
    double dAB = dist(ax, ay, bx, by);
    if (dAB < 1e-6 || dAB > 2.0*rs) return;
    double mx = (ax+bx)*0.5, my = (ay+by)*0.5;
    double ux = (bx-ax)/dAB, uy = (by-ay)/dAB;
    double px = -uy, py = ux;
    double h = std::sqrt(rs*rs - (dAB*0.5)*(dAB*0.5));
    double o1x = mx+h*px, o1y = my+h*py, o2x = mx-h*px, o2y = my-h*py;
    auto makeTarget = [&](double ox, double oy, double &tx, double &ty) {
        double dx = ox-mx, dy = oy-my, len = std::sqrt(dx*dx+dy*dy);
        if (len < 1e-9) { tx = ty = -1; return; }
        tx = ox + (dx/len)*rs; ty = oy + (dy/len)*rs;
    };
    makeTarget(o1x, o1y, t1x, t1y);
    makeTarget(o2x, o2y, t2x, t2y);
}

// ════════════════════════════════ Stigmergy Phase 1 ══════════════════════════

void SensorNode::markCoverage(double cx, double cy, double deposit)
{
    for (int gx = 0; gx < gridRes; gx++)
        for (int gy = 0; gy < gridRes; gy++) {
            double ccx = (gx+0.5)*cellSize, ccy = (gy+0.5)*cellSize;
            if (dist(ccx, ccy, cx, cy) <= rs)
                stigGrid[gx][gy] = std::min(1.0, stigGrid[gx][gy] + deposit);
        }
}

double SensorNode::localCoverage(double cx, double cy) const
{
    int total = 0, marked = 0;
    for (int gx = 0; gx < gridRes; gx++)
        for (int gy = 0; gy < gridRes; gy++) {
            double ccx = (gx+0.5)*cellSize, ccy = (gy+0.5)*cellSize;
            if (dist(ccx, ccy, cx, cy) <= rs) {
                total++;
                if (stigGrid[gx][gy] >= PHEROMONE_SENSE_THRESH) marked++;
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
        if (kv.second < EDGE_PRUNE_THRESH) toPrune.push_back(kv.first);
    }
    for (int id : toPrune) edgeWeight.erase(id);
}

// ════════════════════════════════ Stigmergy Phase 2 ══════════════════════════

bool SensorNode::isPerimeterRedundant() const
{
    cModule *parent = getParentModule();
    int n = parent->par("numSensors");
    std::vector<std::pair<double,double>> liveActive;
    for (int i = 0; i < n; i++) {
        if (i == nodeId) continue;
        SensorNode *s = check_and_cast<SensorNode*>(parent->getSubmodule("sensor", i));
        if (s->state == ON && s->getCurrentEnergy() > LOW_BATTERY_THRESH * s->initEnergy)
            liveActive.push_back({s->posX, s->posY});
    }
    if (liveActive.empty()) return false;
    const double TWO_PI = 2.0 * M_PI;
    for (int i = 0; i < PERIM_SAMPLES; i++) {
        double angle = TWO_PI * i / PERIM_SAMPLES;
        double px = std::max(0.0, std::min(areaSize, posX + rs*std::cos(angle)));
        double py = std::max(0.0, std::min(areaSize, posY + rs*std::sin(angle)));
        bool covered = false;
        for (auto &a : liveActive)
            if (dist(px, py, a.first, a.second) <= rs) { covered = true; break; }
        if (!covered) return false;
    }
    return true;
}

// ════════════════════════════════ Greedy-MSC Phase 1: Clustering ═════════════

void SensorNode::broadcastClusterInvite()
{
    ClusterInviteMsg tmpl;
    tmpl.senderId      = nodeId;
    tmpl.clusterHeadId = myClusterId;
    tmpl.avgDist       = myAvgActiveDist;

    // Send only to active neighbours (within rc)
    for (auto &kv : neighbours) {
        // Check if this neighbour is active
        bool isActive = false;
        for (auto &info : knownActiveInfo)
            if (info.id == kv.first) { isActive = true; break; }
        if (isActive)
            sendDirect(tmpl.dup(), getParentModule()->getSubmodule("sensor", kv.first), "radioIn");
    }
}

void SensorNode::handleClusterInvite(ClusterInviteMsg *m)
{
    if (!useGreedyMSC) { delete m; return; }
    if (myClusterId >= 0) { delete m; return; }   // already claimed

    if (m->directInvite) {
        // EW-based direct invite from a backbone node — accept immediately.
        // Because backbone nodes are staggered by hub score, the first invite
        // received is from the highest hub-score parent: the strongest claim.
        // Do NOT cascade — cluster membership is defined only by the backbone's
        // direct children, not their neighbours.
        myClusterId = m->clusterHeadId;
        cancelTimer(clusterSeedTimer);
        EV_INFO << "[Node " << nodeId << "] joined EW-cluster " << myClusterId
                << " (direct invite from " << m->senderId
                << " ew=" << m->offeredEdgeWeight << ")\n";
    } else {
        // Fallback avg-distance cascade invite (from non-backbone nodes)
        if (state != ON) { delete m; return; }
        double myAvg = 0; int cnt = 0;
        for (auto &info : knownActiveInfo) {
            double d = dist(posX, posY, info.x, info.y);
            if (d <= rs) { myAvg += d; cnt++; }
        }
        if (cnt > 0) myAvg /= cnt;
        myAvgActiveDist = myAvg;
        double maxD = std::max(myAvg, m->avgDist);
        double diff = (maxD > 1e-9) ? std::fabs(myAvg - m->avgDist) / maxD : 0;
        if (diff <= clusterDistThreshold) {
            myClusterId = m->clusterHeadId;
            cancelTimer(clusterSeedTimer);
            EV_INFO << "[Node " << nodeId << "] joined fallback cluster "
                    << myClusterId << " diff=" << diff*100.0 << "%\n";
            broadcastClusterInvite();
        }
    }
    delete m;
}

void SensorNode::broadcastClusterAnnounce()
{
    ClusterAnnounceMsg tmpl;
    tmpl.senderId  = nodeId;
    tmpl.clusterId = myClusterId;

    // Broadcast to ALL neighbours (active + sleep)
    for (auto &kv : neighbours)
        sendDirect(tmpl.dup(), getParentModule()->getSubmodule("sensor", kv.first), "radioIn");
}

void SensorNode::handleClusterAnnounce(ClusterAnnounceMsg *m)
{
    if (useGreedyMSC)
        activeNodeClusterMap[m->senderId] = m->clusterId;
    delete m;
}

// ════════════════════════════════ Greedy-MSC Phase 4: Wake-up ════════════════

void SensorNode::scheduleLowBatteryWarning()
{
    double drainRate = 4.0 / 100.0;
    double warnEnergy;

    if (greedyW >= 1.0) {
        // w=1: warn at 5% of initial energy (run until nearly dead)
        warnEnergy = LOW_BATTERY_THRESH * initEnergy;
    } else {
        // w<1: warn when (1-w) fraction of current energy remains
        // e.g. w=0.5 → warn at 50% of start energy → hand over at half-life
        warnEnergy = phaseStartEnergy * (1.0 - greedyW);
    }

    double timeToWarning = (phaseStartEnergy - warnEnergy) / drainRate;
    if (timeToWarning > 0) {
        cancelTimer(lowBatteryTimer);
        lowBatteryTimer = new cMessage("lowbat", KIND_LOW_BATTERY);
        scheduleAt(simTime() + timeToWarning, lowBatteryTimer);
    }
}

void SensorNode::broadcastLowBatteryWarning(int coveredId, int clustId, double deathTime)
{
    LowBatteryWarningMsg tmpl;
    tmpl.senderId = nodeId;
    tmpl.senderX = posX; tmpl.senderY = posY;
    tmpl.coveredActiveId = coveredId;
    tmpl.clusterId = clustId;
    tmpl.estimatedDeathTime = deathTime;

    for (auto &kv : neighbours)
        sendDirect(tmpl.dup(), getParentModule()->getSubmodule("sensor", kv.first), "radioIn");
}

void SensorNode::handleLowBatteryWarning(LowBatteryWarningMsg *m)
{
    if (!useGreedyMSC || state != OFF || myGroupInCluster < 1) { delete m; return; }

    // Only respond to warnings from my own cluster
    if (m->clusterId != myClusterId) { delete m; return; }

    warningsReceived[m->coveredActiveId]++;
    int warnCount = warningsReceived[m->coveredActiveId];

    // Fix: use <= instead of == so that groups with multiple nodes (G1=2, G2=2 etc.)
    // don't cause the warning counter to race past higher group numbers.
    // The wakeupTimer == nullptr guard prevents double-scheduling if the counter
    // jumps multiple steps in rapid succession.
    if (assignedActiveNodeId == m->coveredActiveId &&
        myGroupInCluster <= warnCount &&
        wakeupTimer == nullptr)
    {
        double wakeupTime = m->estimatedDeathTime - WAKEUP_BUFFER;
        if (wakeupTime <= simTime().dbl()) wakeupTime = simTime().dbl() + 0.1;

        wakeupTimer = new cMessage("wakeup", KIND_WAKEUP);
        scheduleAt(wakeupTime, wakeupTimer);

        EV_INFO << "[Node " << nodeId << "] Group " << myGroupInCluster
                << " cluster=" << myClusterId
                << " wake at t=" << wakeupTime
                << " for active " << assignedActiveNodeId
                << " (warnCount=" << warnCount << ")\n";
    }
    delete m;
}

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

    if (useGreedyMSC) scheduleLowBatteryWarning();
}

// ════════════════════════════════ Reselection ════════════════════════════════

void SensorNode::triggerReselection()
{
    reselectCount++;
    lastReselectTime = simTime();

    std::cout << "\n>>> RESELECTION #" << reselectCount
              << " triggered at t=" << simTime() << "s <<<\n" << std::endl;

    cModule *parent = getParentModule();
    int n = parent->par("numSensors");
    for (int i = 0; i < n; i++) {
        ReselectTriggerMsg *trigger = new ReselectTriggerMsg("RESELECT");
        trigger->senderId = nodeId;
        sendDirect(trigger, parent->getSubmodule("sensor", i), "radioIn");
    }
}

void SensorNode::handleReselectTrigger(ReselectTriggerMsg *m)
{
    delete m;
    if (!useReselection) return;
    if (state != OFF) return;
    if (getCurrentEnergy() <= LOW_BATTERY_THRESH * initEnergy) return;

    state = UNDECIDED;
    hasPred = false;
    bestTime = 1e9; bestHasPred = false; bestPredId = -1;
    cascadeParentId = -1;

    stigGrid.assign(gridRes, std::vector<double>(gridRes, 0.0));

    double lostFrac = 1.0 - getCurrentEnergy() / initEnergy;
    lostFrac = std::max(0.0, std::min(lostFrac, 0.9));
    double vt = exponential(0.05 + lostFrac * 0.40);
    vt = std::min(vt, Ts * 0.9);
    cancelTimer(volunteerTimer);
    volunteerTimer = new cMessage("vol", KIND_VOLUNTEER);
    scheduleAt(simTime() + vt, volunteerTimer);

    cancelTimer(reselectSteadyTimer);
    reselectSteadyTimer = new cMessage("resel_steady", KIND_RESELECT_STEADY);
    scheduleAt(simTime() + Ts, reselectSteadyTimer);

    cancelTimer(reselectRedundTimer);
    reselectRedundTimer = new cMessage("resel_redund", KIND_RESELECT_REDUND);
    scheduleAt(simTime() + Ts + Te + uniform(0, Te), reselectRedundTimer);

    updateDisplay();
}

// ════════════════════════════════ Coverage & Display ═════════════════════════

double SensorNode::getCurrentEnergy() const
{
    double rate = (state == ON) ? 4.0 / 100.0 : 0.01 / 100.0;
    double elapsed = (simTime() - phaseStartTime).dbl();
    return std::max(0.0, phaseStartEnergy - elapsed * rate);
}

void SensorNode::computeCoverage(bool printReport)
{
    cModule *parent = getParentModule();
    int n = parent->par("numSensors");

    std::vector<std::pair<double,double>> active;
    int onCount = 0, offCount = 0, deadCount = 0;
    double sumEnergy = 0.0, maxEnergy = 0.0, minEnergy = 1e9;
    int maxNode = -1, minNode = -1, avgNode = -1;
    double avgDiff = 1e9;

    // Cluster + group stats
    std::map<int, int> clusterSizes;
    std::map<int, std::map<int,int>> clusterGroupCounts;

    for (int i = 0; i < n; i++) {
        SensorNode *s = check_and_cast<SensorNode*>(parent->getSubmodule("sensor", i));
        double liveE = s->getCurrentEnergy();
        sumEnergy += liveE;
        if (liveE > maxEnergy) { maxEnergy = liveE; maxNode = i; }
        if (liveE < minEnergy) { minEnergy = liveE; minNode = i; }
        if (liveE <= LOW_BATTERY_THRESH * s->initEnergy) {
            deadCount++;
        } else if (s->state == ON) {
            active.push_back({s->posX, s->posY}); onCount++;
        } else if (s->state == OFF) {
            offCount++;
        }

        if (useGreedyMSC && s->myClusterId >= 0) {
            clusterSizes[s->myClusterId]++;
            if (s->state == OFF && s->myGroupInCluster > 0 &&
                liveE > LOW_BATTERY_THRESH * s->initEnergy)
                clusterGroupCounts[s->myClusterId][s->myGroupInCluster]++;
        }
    }

    const int G = (int)(areaSize * 5);
    const double cs = 0.2;
    int covered = 0;
    for (int gx = 0; gx < G; gx++)
        for (int gy = 0; gy < G; gy++) {
            double cx = (gx+0.5)*cs, cy = (gy+0.5)*cs;
            for (auto &p : active)
                if (dist(p.first, p.second, cx, cy) <= rs) { covered++; break; }
        }

    double avgEnergy = sumEnergy / n;
    for (int i = 0; i < n; i++) {
        SensorNode *s = check_and_cast<SensorNode*>(parent->getSubmodule("sensor", i));
        double diff = std::fabs(s->getCurrentEnergy() - avgEnergy);
        if (diff < avgDiff) { avgDiff = diff; avgNode = i; }
    }

    double cov = (double)covered / ((long long)G * G);

    if (printReport) {
        EV_INFO << "\n+---------- ROUND " << roundCount << " REPORT ----------+\n"
                << " Active (ON)  : " << onCount   << "\n"
                << " Sleeping(OFF): " << offCount  << "\n"
                << " Dead/Low bat : " << deadCount << "\n"
                << " Avg energy   : " << avgEnergy << "\n"
                << " Coverage     : " << cov*100.0 << "%\n"
                << " Sim time     : " << simTime() << "s\n";
        if (useReselection) {
            EV_INFO << " Reselections : " << reselectCount;
            if (reselectGaveUp) EV_INFO << " (gave up)";
            EV_INFO << "\n";
        }
        if (useGreedyMSC) {
            EV_INFO << " Greedy-MSC   : ON  w=" << greedyW << "\n"
                    << " Clusters     : " << clusterSizes.size() << "\n";
            for (auto &cs : clusterSizes) {
                EV_INFO << "   Cluster " << cs.first << " (" << cs.second << " nodes)";
                auto git = clusterGroupCounts.find(cs.first);
                if (git != clusterGroupCounts.end())
                    for (auto &g : git->second)
                        EV_INFO << "  G" << g.first << "=" << g.second;
                EV_INFO << "\n";
            }
        }
        EV_INFO << "+--------------------------------------+\n\n";

        std::cout << "\n+---------- ROUND " << roundCount << " REPORT ----------+\n" << std::endl;
        std::cout << " Active (ON)  : " << onCount   << std::endl;
        std::cout << " Sleeping(OFF): " << offCount  << std::endl;
        std::cout << " Dead/Low bat : " << deadCount << std::endl;
        std::cout << " Avg energy   : " << avgEnergy << std::endl;
        std::cout << " Coverage     : " << cov*100.0 << "%" << std::endl;
        std::cout << " Sim time     : " << simTime() << "s" << std::endl;
        if (useReselection) {
            std::cout << " Reselections : " << reselectCount;
            if (reselectGaveUp) std::cout << " (gave up)";
            std::cout << std::endl;
        }
        if (useGreedyMSC) {
            std::cout << " Greedy-MSC   : ON  w=" << greedyW << std::endl;
            std::cout << " Clusters     : " << clusterSizes.size() << std::endl;
            for (auto &cs : clusterSizes) {
                std::cout << "   Cluster " << cs.first << " (" << cs.second << " nodes)";
                auto git = clusterGroupCounts.find(cs.first);
                if (git != clusterGroupCounts.end())
                    for (auto &g : git->second)
                        std::cout << "  G" << g.first << "=" << g.second;
                std::cout << std::endl;
            }
        }
        std::cout << "+--------------------------------------+\n" << std::endl;
        std::cout.flush();
    }

    emit(sigCoverage, cov);
    emit(sigActive, (double)onCount);

    // Reselection trigger
    if (useReselection && !reselectGaveUp && cov < coverageThreshold) {
        double statsInterval = par("statsInterval").doubleValue();
        double cooldown = std::max(RESELECT_COOLDOWN_BASE, Ts + 3.0*Te + 2.0*statsInterval);
        bool cooldownPassed = (lastReselectTime < 0 ||
                               (simTime() - lastReselectTime).dbl() > cooldown);

        if (cooldownPassed) {
            // Give up only when there are genuinely no sleeping nodes left
            // to recruit. The per-node reselectCount counter is unreliable
            // as a give-up signal because the coordinator node changes when
            // the lowest-alive-index node dies, resetting the counter.
            // offCount == 0 is the true termination condition: if no OFF
            // nodes exist, no reselection cascade can help.
            if (offCount == 0) {
                reselectGaveUp = true;
                std::cout << "\n>>> RESELECTION gave up: no sleeping nodes "
                          << "remain to recruit <<<\n" << std::endl;
                return;
            }
            if (reselectCount >= maxReselections) {
                reselectGaveUp = true;
                std::cout << "\n>>> RESELECTION limit (" << maxReselections
                          << ") reached <<<\n" << std::endl;
                return;
            }
            covAfterReselect = cov;
            triggerReselection();
        }
    } else if (useReselection && cov >= coverageThreshold && reselectCount > 0) {
        covAfterReselect = cov;
    }
}

void SensorNode::updateDisplay()
{
    if (!hasGUI()) return;
    const double S = 10.0;
    char buf[256];
    if (getCurrentEnergy() <= LOW_BATTERY_THRESH * initEnergy)
        std::snprintf(buf, sizeof(buf), "p=%.0f,%.0f;b=10,10,oval,#000000,#000000,1", posX*S, posY*S);
    else if (state == ON)
        std::snprintf(buf, sizeof(buf), "p=%.0f,%.0f;b=10,10,oval,#00cc00,#006600,1;r=%.0f,#00cc00",
                      posX*S, posY*S, rs*S);
    else if (state == OFF)
        std::snprintf(buf, sizeof(buf), "p=%.0f,%.0f;b=8,8,oval,#cc0000,black,1", posX*S, posY*S);
    else
        std::snprintf(buf, sizeof(buf), "p=%.0f,%.0f;b=8,8,oval,#aaaaaa,black,1", posX*S, posY*S);
    getDisplayString().parse(buf);
}

void SensorNode::finish()
{
    cModule *par = getParentModule();
    int n = par->par("numSensors");
    int minAlive = -1;
    for (int i = 0; i < n; i++) {
        SensorNode *s = check_and_cast<SensorNode*>(par->getSubmodule("sensor", i));
        if (s->getCurrentEnergy() > 0.0) { minAlive = i; break; }
    }
    if (nodeId == minAlive) {
        std::cout << "\n>>> SIMULATION FINISHED after "
                  << roundCount << " rounds <<<" << std::endl;
        computeCoverage();
    }
}
