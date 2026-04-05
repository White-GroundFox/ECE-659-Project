# Bio-Inspired WSN Coverage Protocol
### *Physarum*-Inspired Stigmergic Coverage with Greedy-MSC Lifetime Extension

A wireless sensor network (WSN) coverage simulation built on **OMNeT++ 6.3.0**,
implementing the OGDC density control algorithm extended with:
- Stigmergic pheromone-based coverage bias
- Edge-weight reinforcement (Physarum vein topology)
- Decentralised Greedy-MSC proactive replacement clustering
- Reactive coverage reselection

---

## Requirements

| Software | Version |
|---|---|
| OMNeT++ | 6.3.0 |
| C++ compiler | C++17 or later |
| Python (optional) | 3.8+ (for result merging script) |
| matplotlib (optional) | Any recent version (for plot generation) |

---

## Project Structure

```
coverage/
├── simulations/
│   ├── OGDCNetwork.ned       Network topology definition
│   ├── SensorNode.ned        Node parameters, signals, gates
│   ├── omnetpp.ini           Simulation configurations
│   ├── merge_vectors.py      Post-processing: merge multi-coordinator result vectors
│   └── results/              OMNeT++ output files (.vec, .sca, .vci)
└── src/
    ├── SensorNode.cc         Main simulation logic (~1200 lines)
    ├── SensorNode.h          Class and message definitions
    └── package.ned           Package declaration
```

---

## Building

Open the project in the OMNeT++ IDE:

1. File → Import → General → Existing Projects into Workspace
2. Select the `coverage/` root folder
3. Right-click project → Build Project

Or from the command line inside the `coverage/` folder:

```bash
opp_makemake -f --deep
make -j4
```

---

## Running

### From the OMNeT++ IDE

1. Open `simulations/omnetpp.ini`
2. Select a configuration from the Run menu:
   - **General** — 100 nodes, default settings
   - **N100** — 100 sensors (quick test)
   - **N300** — 300 sensors
   - **N500** — 500 sensors
   - **N1000** — 1000 sensors (paper environment)
3. Click Run or Debug


---

## Key Configuration Parameters

Edit `simulations/omnetpp.ini` to change behaviour.

### Core OGDC Parameters

| Parameter | Default | Description |
|---|---|---|
| `numSensors` | 100 | Number of deployed sensor nodes |
| `sensingRange` | 10 m | Sensing radius $r_s$ |
| `radioRange` | 20 m | Radio range $r_c = 2r_s$ |
| `areaSize` | 50 m | Field size ($A \times A$ m²) |
| `initialEnergy` | 100 | Energy units per node |
| `roundTime` | 1200 s | Duration of one operational round |
| `Ts` | 1 s | OGDC selection window |
| `Td` | 10 ms | Cascade propagation delay |
| `Te` | 0.2 s | Pruning jitter window |

### Lifetime Extension: Greedy-MSC

Set `useGreedyMSC = true` to enable. **Do not combine with `useReselection = true`.**

| Parameter | Default | Description |
|---|---|---|
| `useGreedyMSC` | false | Enable 4-phase Greedy-MSC clustering |
| `warmupRounds` | 10 | Short OGDC-only rounds before clustering activates |
| `warmupRoundTime` | 3 s | Duration of each warmup round |
| `ewClusterThreshold` | 0.5 | Child qualifies if its edge weight ≥ threshold × best child EW |
| `clusterDistThreshold` | 0.25 | Fallback: join cluster if avg-distance differs by ≤ 25% |
| `greedyW` | 1.0 | Active time weight (1.0 = disjoint groups, 0.5 = overlapping) |

> **Important:** `roundTime` must need to be larger than `sim-time-limit` to test baseline, reselection mode, and greedy MSC mode.
> !!!! A ratio of 5–8× is recommended (e.g. `roundTime = 19000s`, `sim-time-limit = 18000s`) to reproduce the result in the report!!!!

### Lifetime Extension: Reactive Reselection

Set `useReselection = true` to enable. **Do not combine with `useGreedyMSC = true`.**

| Parameter | Default | Description |
|---|---|---|
| `useReselection` | false | Enable reactive coverage reselection |
| `coverageThreshold` | 0.90 | Trigger reselection when coverage drops below this |
| `maxReselections` | 50 | Safety limit on reselection attempts |
| `statsInterval` | 10 s | Coverage measurement interval |

---

## Protocol Modes

Three configurations are available by combining the above flags:

| Mode | `useGreedyMSC` | `useReselection` | Description |
|---|---|---|---|
| **Baseline** | false | false | Stigmergic OGDC only, no lifetime extension |
| **Greedy-MSC** | true | false | Proactive replacement via edge-weight clustering |
| **Reselection** | false | true | Reactive cascade re-election on coverage drop |

---

## Reading the Console Output

The coordinator node (lowest-alive index) prints a round report to stdout:

```
+---------- ROUND 2 REPORT ----------+
 Active (ON)  : 18
 Sleeping(OFF): 82
 Dead/Low bat : 0
 Avg energy   : 91.25
 Coverage     : 99.07%
 Sim time     : 1201.6s
 Greedy-MSC   : ON  w=1.0
 Clusters     : 5
   Cluster 6  (39 nodes)  G1=2  G2=2  G3=6 ...
+--------------------------------------+
```

**Coverage** is measured on a 250×250 grid (cell size 0.2 m²).
**Clusters** are only shown when `useGreedyMSC = true`.
**G1=2** means 2 sleeping nodes are assigned to Group 1 of that cluster.

Edge weights per node are also printed before each selection phase:
```
 Edge weights per node (parent → child : weight):
   Node  91 -> Node   4 : 0.7321
   Node   9 -> Node  14 : 0.5681
   ...
```
High-weight edges form the Physarum vein backbone used by Greedy-MSC clustering.

---

## GUI Visualisation (Greedy-MSC mode)

When running with the OMNeT++ Qtenv GUI and `useGreedyMSC = true`:

| Visual | Meaning |
|---|---|
| Large coloured circle + sensing ring | Active (ON) node |
| Small coloured dot (same hue) | Sleeping backup node in same cluster |
| Thicker dark border | Higher accumulated edge weight (vein node) |
| Black dot | Dead node (energy ≤ 5% of initial) |
| Grey circle | UNDECIDED node (during selection phase) |

Each cluster is assigned a unique colour from a 15-colour palette. The
border thickness (1–5 px) encodes the node's maximum edge weight,
revealing the structural backbone of the cascade topology.

---

## Post-Processing: Merging Result Vectors

Coverage and active-node statistics are recorded from whichever node
is currently the lowest-alive coordinator. As nodes die, the coordinator
changes, creating multiple separate vector records in the `.vec` file.
Use the provided script to merge them into a single continuous time series:

```bash
cd simulations/results

# Basic merge to CSV
python3 merge_vectors.py N100-#0.vec merged.csv

# Merge and generate publication-quality PNG
python3 merge_vectors.py N100-#0.vec merged.csv --plot --title "Greedy-MSC N=100"
```

The script outputs:
- A CSV with columns: `time_s`, `coverage_ratio`, `active_nodes`
- A 300 DPI PNG with coverage and active-node subplots (if `--plot` is given)

---

## Tuning Constants (SensorNode.cc)

The following constants at the top of `SensorNode.cc` control the
stigmergy and edge-weight behaviour and can be adjusted without
recompiling the NED files:

| Constant | Value | Effect |
|---|---|---|
| `PHEROMONE_DEPOSIT` | 1.00 | Amount deposited per CoverageMark received |
| `PHEROMONE_DECAY` | 0.80 | Multiplier applied at mid-round |
| `PHEROMONE_SENSE_THRESH` | 0.25 | Minimum grid value to count as "covered" |
| `COVERAGE_BIAS` | 0.25 | Maximum stigmergic delay as fraction of $T_s$ |
| `EDGE_REINFORCE` | 0.50 | Added to edge weight on each cascade activation |
| `EDGE_DECAY` | 0.85 | Multiplier applied to all edges at mid-round |
| `EDGE_PRUNE_THRESH` | 0.10 | Edges below this are discarded |
| `WAKEUP_BUFFER` | 30.0 s | Wake-up lead time before estimated node death |
| `LOW_BATTERY_THRESH` | 0.05 | Fraction of initial energy below which node is "dead" |

---

## Authors

- Apisan S Kaneshan — University of Waterloo
- Guanfeng Wu — University of Waterloo
- Luna Nga Do — University of Waterloo

*ECE Department, University of Waterloo, 2026*
