#ifndef SENSORNODE_H_
#define SENSORNODE_H_

#pragma once
#include <omnetpp.h>
#include <map>
#include <set>
#include <vector>
#include <cmath>
using namespace omnetpp;

// ── POWER_ON network message ─────────────────────────────────────────────────
class PowerOnMsg : public cMessage {
  public:
    int    senderId = -1;
    double senderX  = 0,  senderY  = 0;
    bool   hasPred  = false;
    double predX    = -1, predY    = -1;
    int    predId   = -1;
    double t1x = -1, t1y = -1;
    double t2x = -1, t2y = -1;

    PowerOnMsg(const char *name = "POWER_ON") : cMessage(name, 100) {}
    virtual PowerOnMsg *dup() const override { return new PowerOnMsg(*this); }
};

// ── COVERAGE_MARK stigmergy message ──────────────────────────────────────────
class CoverageMarkMsg : public cMessage {
  public:
    int    senderId        = -1;
    double senderX         = 0, senderY = 0;
    double sensingRange    = 0;
    double remainingEnergy = 0;

    CoverageMarkMsg(const char *name = "COV_MARK") : cMessage(name, 101) {}
    virtual CoverageMarkMsg *dup() const override { return new CoverageMarkMsg(*this); }
};

// ── CLUSTER_INVITE message (Phase 1) ─────────────────────────────────────────
// Broadcast by an active node to invite active neighbours into its cluster.
// Cascades: a node that joins re-broadcasts with its own avgDist.
class ClusterInviteMsg : public cMessage {
  public:
    int    senderId      = -1;
    int    clusterHeadId = -1;   // original seed node of this cluster
    double avgDist       = 0;    // sender's average distance to its active neighbours

    ClusterInviteMsg(const char *name = "CLUSTER_INV") : cMessage(name, 105) {}
    virtual ClusterInviteMsg *dup() const override { return new ClusterInviteMsg(*this); }
};

// ── CLUSTER_ANNOUNCE message (Phase 1 → Phase 2 bridge) ─────────────────────
// Broadcast by each active node after clustering settles to tell all
// neighbours (active + sleep) which cluster it belongs to.
class ClusterAnnounceMsg : public cMessage {
  public:
    int senderId  = -1;
    int clusterId = -1;

    ClusterAnnounceMsg(const char *name = "CLUSTER_ANN") : cMessage(name, 106) {}
    virtual ClusterAnnounceMsg *dup() const override { return new ClusterAnnounceMsg(*this); }
};

// ── LOW_BATTERY_WARNING message (Phase 4) ────────────────────────────────────
class LowBatteryWarningMsg : public cMessage {
  public:
    int    senderId           = -1;
    double senderX            = 0, senderY = 0;
    int    coveredActiveId    = -1;
    int    clusterId          = -1;
    double estimatedDeathTime = 0;

    LowBatteryWarningMsg(const char *name = "LOW_BAT") : cMessage(name, 103) {}
    virtual LowBatteryWarningMsg *dup() const override { return new LowBatteryWarningMsg(*this); }
};

// ── RESELECT_TRIGGER message ─────────────────────────────────────────────────
class ReselectTriggerMsg : public cMessage {
  public:
    int senderId = -1;

    ReselectTriggerMsg(const char *name = "RESELECT") : cMessage(name, 104) {}
    virtual ReselectTriggerMsg *dup() const override { return new ReselectTriggerMsg(*this); }
};

// ── Self-message kinds ────────────────────────────────────────────────────────
enum SelfKind {
    KIND_ROUND            = 1,
    KIND_VOLUNTEER        = 2,
    KIND_SELECT           = 3,
    KIND_STEADY           = 4,
    KIND_STATS            = 5,
    KIND_PRUNE            = 6,
    KIND_REDUNDANCY_CHECK = 7,
    KIND_DISPLAY_REFRESH  = 8,
    KIND_GROUPING         = 10,  // Phase 2+3: sleep node cluster assign + Greedy-MSC
    KIND_LOW_BATTERY      = 11,
    KIND_WAKEUP           = 12,
    KIND_RESELECT_STEADY  = 13,
    KIND_RESELECT_REDUND  = 14,
    KIND_CLUSTER_SEED     = 15,  // Phase 1: active node starts/joins clustering
    KIND_CLUSTER_DONE     = 16   // Phase 1 end: finalize + announce cluster IDs
};

// ── Neighbour record ──────────────────────────────────────────────────────────
struct Neighbour { int id; double x, y; };

// ── Active neighbour info ─────────────────────────────────────────────────────
struct ActiveNeighbourInfo { int id; double x, y, energy; };

// ── SensorNode ────────────────────────────────────────────────────────────────
class SensorNode : public cSimpleModule
{
  public:
    enum State { UNDECIDED, ON, OFF } state = UNDECIDED;
    double posX = 0, posY = 0, rs = 10;
    double getCurrentEnergy() const;
    double energy = 0, initEnergy = 1;

  private:
    int    nodeId, numNodes;
    double rc, areaSize;
    double roundTime, Td, Ts, Te;
    int    roundCount = 0;
    int    lastReportedRound = 0;
    simtime_t roundStartTime   = 0;
    simtime_t phaseStartTime   = 0;
    double    phaseStartEnergy = 0;

    std::map<int,Neighbour> neighbours;

    // ── Stigmergy Phase 1 — pheromone grid ───────────────────────────────────
    std::vector<std::vector<double>> stigGrid;
    int    gridRes  = 10;
    double cellSize = 5.0;

    std::map<int,double> edgeWeight;
    int cascadeParentId = -1;

    // ── Stigmergy Phase 2 — perimeter pruning ────────────────────────────────
    std::vector<std::pair<double,double>> knownActive;

    // ── Greedy-MSC — shared state ────────────────────────────────────────────
    bool   useGreedyMSC          = false;
    double clusterDistThreshold  = 0.25;
    double greedyW               = 1.0;
    std::vector<ActiveNeighbourInfo> knownActiveInfo;

    // ── Phase 1: active node clustering ──────────────────────────────────────
    int    myClusterId       = -1;
    double myAvgActiveDist   = 0;
    std::map<int, int> activeNodeClusterMap;   // activeNodeId → clusterId (from announce)

    // ── Phase 2+3: sleep node assignment + Greedy-MSC ────────────────────────
    int  myGroupInCluster     = -1;
    int  assignedActiveNodeId = -1;

    // ── Phase 4: wake-up chain ───────────────────────────────────────────────
    std::map<int, int> warningsReceived;

    // ── Reselection — reactive coverage restoration ──────────────────────────
    bool   useReselection      = false;
    double coverageThreshold   = 0.90;
    int    maxReselections     = 50;
    simtime_t lastReselectTime = -1;
    int    reselectCount       = 0;
    double covAfterReselect    = -1;
    bool   reselectGaveUp      = false;

    // ── Cascade candidate tracking ────────────────────────────────────────────
    double bestTime    = 1e9;
    double bestPredX   = -1, bestPredY = -1;
    bool   bestHasPred = false;
    int    bestPredId  = -1;

    double predX = -1, predY = -1;
    bool   hasPred = false;

    // ── Timers ────────────────────────────────────────────────────────────────
    cMessage *roundTimer          = nullptr;
    cMessage *volunteerTimer      = nullptr;
    cMessage *selectTimer         = nullptr;
    cMessage *steadyTimer         = nullptr;
    cMessage *statsTimer          = nullptr;
    cMessage *pruneTimer          = nullptr;
    cMessage *redundancyTimer     = nullptr;
    cMessage *refreshTimer        = nullptr;
    cMessage *clusterSeedTimer    = nullptr;
    cMessage *clusterDoneTimer    = nullptr;
    cMessage *groupingTimer       = nullptr;
    cMessage *lowBatteryTimer     = nullptr;
    cMessage *wakeupTimer         = nullptr;
    cMessage *reselectSteadyTimer = nullptr;
    cMessage *reselectRedundTimer = nullptr;

    simsignal_t sigCoverage, sigActive;

    // ── OGDC helpers ──────────────────────────────────────────────────────────
    void startRound();
    void handlePowerOn(PowerOnMsg *m);
    void turnOn();
    void turnOff();
    void broadcastPowerOn();
    void computeTargets(double ax, double ay, double bx, double by,
                        double &t1x, double &t1y,
                        double &t2x, double &t2y);

    // ── Stigmergy helpers ─────────────────────────────────────────────────────
    void   markCoverage(double cx, double cy, double deposit = 1.0);
    double localCoverage(double cx, double cy) const;
    void   decayAndPruneEdges();
    void   broadcastCoverageMark();
    void   handleCoverageMark(CoverageMarkMsg *m);
    bool   isPerimeterRedundant() const;

    // ── Greedy-MSC Phase 1: clustering ────────────────────────────────────────
    void handleClusterInvite(ClusterInviteMsg *m);
    void handleClusterAnnounce(ClusterAnnounceMsg *m);
    void broadcastClusterInvite();
    void broadcastClusterAnnounce();

    // ── Greedy-MSC Phase 4: wake-up chain ─────────────────────────────────────
    void handleLowBatteryWarning(LowBatteryWarningMsg *m);
    void broadcastLowBatteryWarning(int coveredId, int clustId, double deathTime);
    void turnOnAsReplacement();
    void scheduleLowBatteryWarning();

    // ── Reselection helpers ───────────────────────────────────────────────────
    void triggerReselection();
    void handleReselectTrigger(ReselectTriggerMsg *m);

    // ── Display & stats ───────────────────────────────────────────────────────
    void computeCoverage(bool printReport = true);
    void updateDisplay();

    void cancelTimer(cMessage *&t) { cancelAndDelete(t); t = nullptr; }

    double dist(double x1, double y1, double x2, double y2) const {
        double dx = x2-x1, dy = y2-y1;
        return std::sqrt(dx*dx + dy*dy);
    }

  protected:
    virtual int  numInitStages() const override { return 2; }
    virtual void initialize(int stage) override;
    virtual void handleMessage(cMessage *msg) override;
    virtual void finish() override;
};

#endif /* SENSORNODE_H_ */
