#ifndef SENSORNODE_H_
#define SENSORNODE_H_

#pragma once
#include <omnetpp.h>
#include <map>
#include <vector>
#include <cmath>
using namespace omnetpp;

// -- POWER_ON network message -------------------------------------------------
class PowerOnMsg : public cMessage {
  public:
    int    senderId = -1;
    double senderX  = 0,  senderY  = 0;
    bool   hasPred  = false;
    double predX    = -1, predY    = -1;
    double t1x = -1, t1y = -1;   // cascade target 1  (-1 = none)
    double t2x = -1, t2y = -1;   // cascade target 2  (-1 = none)

    PowerOnMsg(const char *name = "POWER_ON") : cMessage(name, 100) {}
    virtual PowerOnMsg *dup() const override { return new PowerOnMsg(*this); }
};

// -- Self-message kinds -------------------------------------------------------
enum SelfKind { KIND_ROUND=1, KIND_VOLUNTEER=2, KIND_SELECT=3, KIND_STEADY=4, KIND_STATS=5 };

// -- Neighbour record ---------------------------------------------------------
struct Neighbour { int id; double x, y; };


// -- SensorNode ---------------------------------------------------------------
class SensorNode : public cSimpleModule
{
  public:
    // Public so node-0 can read them during coverage computation
    enum State { UNDECIDED, ON, OFF } state = UNDECIDED;
    double posX = 0, posY = 0, rs = 10;

  private:
    int    nodeId, numNodes;
    double rc, energy, initEnergy, areaSize;
    double roundTime, Td, Ts, Te;
    int    roundCount = 0;

    std::map<int,Neighbour> neighbours;

    // Positions of active nodes this node has learned about this round
    // (populated from received POWER_ON messages; used for redundancy check)
    std::vector<std::pair<double,double>> knownActive;

    // Best pending cascade candidate (updated on each received POWER_ON)
    double bestTime    = 1e9;
    double bestPredX   = -1, bestPredY = -1;
    bool   bestHasPred = false;

    // Predecessor confirmed when this node turns ON
    double predX = -1, predY = -1;
    bool   hasPred = false;

    // Timers
    cMessage *roundTimer    = nullptr;
    cMessage *volunteerTimer= nullptr;
    cMessage *selectTimer   = nullptr;
    cMessage *steadyTimer   = nullptr;
    cMessage *statsTimer    = nullptr;

    simsignal_t sigCoverage, sigActive;

    // -- Helpers
    void startRound();
    void handlePowerOn(PowerOnMsg *m);
    void turnOn();
    void turnOff();
    void broadcastPowerOn();
    void computeTargets(double ax,double ay,double bx,double by,
                        double &t1x,double &t1y,
                        double &t2x,double &t2y);
    void computeCoverage();
    // Returns true if this node's sensing area is already fully covered
    // by the active nodes it has learned about — i.e., activating would
    // add zero new coverage (redundant node).
    bool isRedundant() const;
    void updateDisplay();
    void cancelTimer(cMessage *&t){ cancelAndDelete(t); t=nullptr; }

    double dist(double x1,double y1,double x2,double y2) const {
        double dx=x2-x1, dy=y2-y1;
        return std::sqrt(dx*dx+dy*dy);
    }

  protected:
    virtual int  numInitStages() const override { return 2; }
    virtual void initialize(int stage) override;
    virtual void handleMessage(cMessage *msg) override;
    virtual void finish() override;
};

#endif /* SENSORNODE_H_ */
