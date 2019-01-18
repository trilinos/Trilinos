
function plotControl()

adj   = load('cell_to_node_quad.txt') + 1;       %% load node adjacency, increment by 1 for 1-based indexing
nodes = load('nodes.txt');                       %% load node coordinates
N     = 3*length(nodes);                         %% determine number of nodes
nt    = numel(dir('control.*.txt'));             %% determine number of time steps

%% Read in state
Ux = cell(nt,1);
Uy = cell(nt,1);
for i = 0:nt-1
  data_obj  = importdata(['control.',int2str(i),'.txt'], ' ', 2);      %% we need to skip the first two lines
  state     = data_obj.data;
  data_obj  = importdata(['map_control.',int2str(i),'.txt'], ' ', 9);  %% we need to skip the first 9 lines
  map_state = data_obj.data;
  map_state = map_state(1:2:end)+1;
  [tmp, state_permute] = sort(map_state);
  state = state(state_permute);  %% we need to permute the state according to parallel maps

  Ux{i+1} = state(1:3:N);  %% Extract x-velocity
  Uy{i+1} = state(2:3:N);  %% Extract y-velocity
end

axsize = 400;
figure('Position', [100 100 4*axsize 4*axsize]);
for i=1:nt-1
  quiver(nodes(:,1), nodes(:,2), Ux{i+1}, Uy{i+1});
  axis('equal','tight');
  xlim([-5,-3]);
  xlabel('x');
  ylabel('y');
  title(['Time/T = ',num2str(i/(nt-1))],'fontsize',16);
  set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');% tightfig;

  drawnow
end
