% rng(12092018);
% N = 30;

generator = @() generateRotatedSimulation(N);
includeLinear = false;

ResultsSummary_rotated = simulationFigures(generator,includeLinear);