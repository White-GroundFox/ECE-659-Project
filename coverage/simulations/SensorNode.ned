



simple SensorNode {
    parameters:
        double sensingRange @unit(m) = default(10m);
        double radioRange   @unit(m) = default(20m);
        double areaSize     @unit(m) = default(50m);
        double initialEnergy          = default(5000);
        double roundTime    @unit(s) = default(1000s);
        double Td           @unit(s) = default(0.01s);
        double Ts           @unit(s) = default(1s);
        double Te           @unit(s) = default(0.2s);
        bool   useGreedyMSC           = default(false);
        double clusterDistThreshold   = default(0.25);  // Phase 1 fallback clustering similarity
        double greedyW                = default(1.0);   // active time weight w
        int    warmupRounds           = default(10);    // rounds before EW clustering activates
        double warmupRoundTime @unit(s) = default(3s); // short round duration during warmup
        double ewClusterThreshold     = default(0.5);  // child qualifies if EW >= threshold * best
        bool   useReselection         = default(false);
        double coverageThreshold      = default(0.90);
        int    maxReselections        = default(50);
        double statsInterval @unit(s) = default(100s);
		
		@signal[coverage](type=double);
		@signal[activeNodes](type=double);
		
        @statistic[coverage](title="Coverage ratio";
            record=vector,last; interpolationmode=none);
        @statistic[activeNodes](title="Active nodes";
            record=vector,last; interpolationmode=none);

        @display("i=misc/sensor2");
    gates:
        input radioIn @directIn;
}