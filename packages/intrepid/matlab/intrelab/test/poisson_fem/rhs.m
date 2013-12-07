function f = rhs(usr_par)

%%%%%%%%%%% get computational mesh cubature points physical frame
cubPointsPhysCoord=usr_par.cubPointsPhysCoord;

% Define the frequencies
freq=pi*2;

%%%%%%%%%%% build right hand side
nVert = size(cubPointsPhysCoord,2);
numCells = size(cubPointsPhysCoord,3);
ff = zeros(nVert,numCells);
ff(:,:) = sin(freq*cubPointsPhysCoord(1,:,:)) .* ...
    sin(freq*cubPointsPhysCoord(2,:,:));

%%%%%%%%%%% integrate right hand side
f = zeros(usr_par.numFields, usr_par.numCells);
intrepid_integrate(f, ff, ...
    usr_par.weighted_transformed_val_at_cub_points, 'COMP_BLAS');
