



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