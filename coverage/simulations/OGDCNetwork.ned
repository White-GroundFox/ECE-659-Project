
network OGDCNetwork {
    parameters:
        int numSensors = default(100);
        @display("bgb=500,500");   // 500px = 50m at 10px/m
    submodules:
        sensor[numSensors]: SensorNode { }
}