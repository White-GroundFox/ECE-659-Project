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

// ── LOW_BATTERY_WARNING message (Greedy-MSC only) ────────────────────────────
class LowBatteryWarningMsg : public cMessage {
  public:
    int    senderId           = -1;
    double senderX            = 0, senderY = 0;
    int    coveredActiveId    = -1;
    double estimatedDeathTime = 0;

    LowBatteryWarningMsg(const char *name = "LOW_BAT") : cMessage(name, 103) {}
    virtual LowBatteryWarningMsg *dup() const override { return new LowBatteryWarningMsg(*this); }
};

// ── RESELECT_TRIGGER message ─────────────────────────────────────────────────
// Broadcast by the coordinator when coverage drops below the threshold.
// OFF nodes with sufficient energy re-enter UNDECIDED and run a new cascade.
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
    KIND_GROUPING         = 10,  // Greedy-MSC: local group assignment
    KIND_LOW_BATTERY      = 11,
    KIND_WAKEUP           = 12,
    KIND_RESELECT_STEADY  = 13,  // Reselection: end of cascade window
    KIND_RESELECT_REDUND  = 14   // Reselection: post-cascade pruning
};

// ── Neighbour record ──────────────────────────────────────────────────────────
struct Neighbour { int id; double x, y; };

// ── Active neighbour info (for Greedy-MSC grouping) ───────────────────────────
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

    // ── Greedy-MSC Phase 3 — replacement grouping ────────────────────────────
    bool useGreedyMSC = false;
    std::vector<ActiveNeighbourInfo> knownActiveInfo;
    int  myGroupNumber        = -1;
    int  assignedActiveNodeId = -1;
    double myGroupScore       = 0;
    std::map<int, int> warningsReceived;

    // ── Reselection — reactive coverage restoration ──────────────────────────
    bool   useReselection      = false;
    double coverageThreshold   = 0.90;
    int    maxReselections     = 50;      // safety limit
    simtime_t lastReselectTime = -1;      // cooldown: prevent rapid retriggering
    int    reselectCount       = 0;       // how many reselections so far
    double covAfterReselect    = -1;      // coverage measured after last reselection settled
    bool   reselectGaveUp      = false;   // true = stop trying, not enough nodes

    // ── Cascade candidate tracking ────────────────────────────────────────────
    double bestTime    = 1e9;
    double bestPredX   = -1, bestPredY = -1;
    bool   bestHasPred = false;
    int    bestPredId  = -1;

    double predX = -1, predY = -1;
    bool   hasPred = false;

    // ── Timers ────────────────────────────────────────────────────────────────
    cMessage *roundTimer         = nullptr;
    cMessage *volunteerTimer     = nullptr;
    cMessage *selectTimer        = nullptr;
    cMessage *steadyTimer        = nullptr;
    cMessage *statsTimer         = nullptr;
    cMessage *pruneTimer         = nullptr;
    cMessage *redundancyTimer    = nullptr;
    cMessage *refreshTimer       = nullptr;
    cMessage *groupingTimer      = nullptr;
    cMessage *lowBatteryTimer    = nullptr;
    cMessage *wakeupTimer        = nullptr;
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

    // ── Stigmergy Phase 1 helpers ─────────────────────────────────────────────
    void   markCoverage(double cx, double cy, double deposit = 1.0);
    double localCoverage(double cx, double cy) const;
    void   decayAndPruneEdges();
    void   broadcastCoverageMark();
    void   handleCoverageMark(CoverageMarkMsg *m);

    // ── Stigmergy Phase 2 helpers ─────────────────────────────────────────────
    bool isPerimeterRedundant() const;

    // ── Greedy-MSC Phase 3 helpers ────────────────────────────────────────────
    void handleLowBatteryWarning(LowBatteryWarningMsg *m);
    void broadcastLowBatteryWarning(int coveredId, double deathTime);
    void turnOnAsReplacement();

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
