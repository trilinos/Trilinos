function soln = soln_sine(cubPointsPhysCoord)

% Define the frequencies
freq=pi*2;

% Retrieve dimensions.
dim = size(cubPointsPhysCoord, 1);
numCubPts = size(cubPointsPhysCoord,2);
numCells = size(cubPointsPhysCoord,3);

soln = ones(numCubPts,numCells);
for i=1:dim
  soln = soln .* squeeze(sin(freq*cubPointsPhysCoord(i,:,:)));
end
