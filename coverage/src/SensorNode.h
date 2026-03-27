#ifndef SENSORNODE_H_
#define SENSORNODE_H_

#pragma once
#include <omnetpp.h>
#include <map>
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
// Broadcast when a node turns ON. Serves dual purpose:
//   Phase 1 (during cascade 0..Ts): updates pheromone grid to bias timers
//   Phase 2 (Ts..Ts+Te):            populates knownActive for perimeter pruning
class CoverageMarkMsg : public cMessage {
  public:
    int    senderId     = -1;
    double senderX      = 0, senderY = 0;
    double sensingRange = 0;

    CoverageMarkMsg(const char *name = "COV_MARK") : cMessage(name, 101) {}
    virtual CoverageMarkMsg *dup() const override { return new CoverageMarkMsg(*this); }
};

// ── Self-message kinds ────────────────────────────────────────────────────────
enum SelfKind {
    KIND_ROUND            = 1,
    KIND_VOLUNTEER        = 2,
    KIND_SELECT           = 3,
    KIND_STEADY           = 4,
    KIND_STATS            = 5,
    KIND_PRUNE            = 6,   // mid-round: pheromone decay + edge pruning
    KIND_REDUNDANCY_CHECK = 7,    // Phase 2: post-selection perimeter self-pruning
    KIND_DISPLAY_REFRESH  = 8

};

// ── Neighbour record ──────────────────────────────────────────────────────────
struct Neighbour { int id; double x, y; };

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
    simtime_t phaseStartTime   = 0;   // when current ON/OFF phase began
    double    phaseStartEnergy = 0;   // energy at the start of current phase

    std::map<int,Neighbour> neighbours;

    // ── Stigmergy Phase 1 — pheromone grid ───────────────────────────────────
    // Local grid tracking "how saturated" each area is from previous/current
    // round activations. Used to add a penalty to the cascade select timer so
    // nodes near already-covered areas lose the timer race to nodes in fresh areas.
    std::vector<std::vector<double>> stigGrid;   // [gridRes][gridRes]
    int    gridRes  = 10;
    double cellSize = 5.0;

    std::map<int,double> edgeWeight;   // cascade edge weights (Steps 3 & 4)
    int cascadeParentId = -1;

    // ── Stigmergy Phase 2 — perimeter pruning ────────────────────────────────
    // After the cascade completes (at Ts), each ON node waits Te seconds while
    // collecting positions of all other active nodes via CoverageMarkMsg.
    // At Ts+Te it runs the 36-point perimeter check: if every boundary point of
    // its sensing disk is covered by a known active neighbour, it turns OFF.
    std::vector<std::pair<double,double>> knownActive;

    // ── Cascade candidate tracking ────────────────────────────────────────────
    double bestTime    = 1e9;
    double bestPredX   = -1, bestPredY = -1;
    bool   bestHasPred = false;
    int    bestPredId  = -1;

    double predX = -1, predY = -1;
    bool   hasPred = false;

    // ── Timers ────────────────────────────────────────────────────────────────
    cMessage *roundTimer      = nullptr;
    cMessage *volunteerTimer  = nullptr;
    cMessage *selectTimer     = nullptr;
    cMessage *steadyTimer     = nullptr;
    cMessage *statsTimer      = nullptr;   // per-round coverage report
    cMessage *pruneTimer      = nullptr;
    cMessage *redundancyTimer = nullptr;   // Phase 2 self-pruning
    cMessage *refreshTimer    = nullptr;

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

    // ── Display & stats ───────────────────────────────────────────────────────
    void computeCoverage();
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
