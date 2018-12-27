% rng(12072018);
% N = 30;

generator = @() generatePlainSimulation(N);
includeLinear = true;

ResultsSummary_plain = simulationFigures(generator,includeLinear);