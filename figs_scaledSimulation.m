% rng(12072018);
% N = 30;

generator = @() generateScaledSimulation(N);
includeLinear = false;

ResultsSummary_scaled = simulationFigures(generator,includeLinear);