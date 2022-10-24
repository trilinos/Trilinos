function plotDarcy(numprocs)

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates
N = length(nodes);

%% PLOT PRESSURE
data_obj = importdata('state.txt', ' ', 2);  %% we need to skip the first two lines
state = data_obj.data;
data_obj = importdata('map_state.txt', ' ', 9);  %% we need to skip the first 9 lines
map_state = data_obj.data;
map_state = map_state(1:2:end)+1;
[tmp, state_permute] = sort(map_state);
state = state(state_permute);  %% we need to permute the state according to parallel maps

state = state(1:N); %% nodal values only

vtkwrite('state', ...
          adj, ...
          nodes, ...
          state, ...
          'Tetrahedron', ...
          'state');

if (numprocs>1)
  for i=1:numprocs
    cid_files{i} = "cellid_"       + num2str(i-1) + "_" + num2str(numprocs) + ".txt";
    per_files{i} = "permeability_" + num2str(i-1) + "_" + num2str(numprocs) + ".txt";
    vel_files{i} = "velocity_"     + num2str(i-1) + "_" + num2str(numprocs) + ".txt";
  end
else
  cid_files{1} = "cellid.txt";
  per_files{1} = "permeability.txt";
  vel_files{1} = "velocity.txt";
end

cids = [];
pers = [];
vels = [];

for i=1:length(cid_files)
  data_obj = importdata(cid_files{i}, ' ');
  cids = [cids; data_obj];
  data_obj = importdata(per_files{i}, ' ');
  pers = [pers; data_obj];
  data_obj = importdata(vel_files{i}, ' ');
  vels = [vels; data_obj];
end

cids = cids + 1;
vtk_pers = pers(cids);
vtk_vels = vels(cids, :);

vtkwritecell('permeability', ...
             adj, ...
             nodes, ...
             vtk_pers, ...
             'Tetrahedron', ...
             'permeability');

vtkwritecell('velocity', ...
             adj, ...
             nodes, ...
             vtk_vels, ...
             'Tetrahedron', ...
             'velocity');
