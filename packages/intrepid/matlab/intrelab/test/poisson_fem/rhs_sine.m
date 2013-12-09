function f = rhs_sine(usr_par)

%%%%%%%%%%% get computational mesh cubature points physical frame
cubPointsPhysCoord=usr_par.cubPointsPhysCoord;

% Define the frequencies
freq=pi*2;

% Retrieve spatial dimension.
dim = usr_par.spaceDim;

%%%%%%%%%%% build right hand side
numCubPts = size(cubPointsPhysCoord,2);
numCells = size(cubPointsPhysCoord,3);

ff = dim*(freq^2)*ones(numCubPts,numCells);
for i=1:dim
  ff = ff .* squeeze(sin(freq*cubPointsPhysCoord(i,:,:)));
end

%%%%%%%%%%% integrate right hand side
f = zeros(usr_par.numFields, usr_par.numCells);
intrepid_integrate(f, ff, ...
    usr_par.weighted_transformed_val_at_cub_points, 'COMP_BLAS');
