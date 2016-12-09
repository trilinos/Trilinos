
adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

data_obj = importdata('state.txt', ' ', 2);  %% we need to skip the first two lines
state = data_obj.data;

stateR = state(1:2:end);
stateI = state(2:2:end);

figure(1)
trisurf(adj, nodes(:,1), nodes(:,2), stateR);
shading interp;
view(0,90)
title('Real State','fontsize',24)
colorbar
axis equal
axis tight

figure(2)
trisurf(adj, nodes(:,1), nodes(:,2), stateI);
shading interp;
view(0,90)
title('Imaginary State','fontsize',24)
colorbar
axis equal
axis tight

figure(3)
trisurf(adj, nodes(:,1), nodes(:,2), sqrt(stateR.^2 + stateI.^2));
shading interp;
view(0,90)
title('State Magnitude','fontsize',24)
colorbar
axis equal
axis tight

data_obj = importdata('control.txt', ' ', 2);  %% we need to skip the first two lines
control = data_obj.data;

controlR = control(1:2:end);
controlI = control(2:2:end);

figure(4)
trisurf(adj, nodes(:,1), nodes(:,2), controlR);
shading interp;
view(0,90)
title('Real Control','fontsize',24)
colorbar
axis equal
axis tight

figure(5)
trisurf(adj, nodes(:,1), nodes(:,2), controlI);
shading interp;
view(0,90)
title('Imaginary Control','fontsize',24)
colorbar
axis equal
axis tight

figure(6)
trisurf(adj, nodes(:,1), nodes(:,2), sqrt(controlR.^2 + controlI.^2));
shading interp;
view(0,90)
title('Control Magnitude','fontsize',24)
colorbar
axis equal
axis tight

data_obj = importdata('state_uncontrolled.txt', ' ', 2);  %% we need to skip the first two lines
state = data_obj.data;

stateR = state(1:2:end);
stateI = state(2:2:end);

figure(7)
trisurf(adj, nodes(:,1), nodes(:,2), stateR);
shading interp;
view(0,90)
title('Real Uncontrolled State','fontsize',24)
colorbar
axis equal
axis tight

figure(8)
trisurf(adj, nodes(:,1), nodes(:,2), stateI);
shading interp;
view(0,90)
title('Imaginary Uncontrolled State','fontsize',24)
colorbar
axis equal
axis tight

figure(9)
trisurf(adj, nodes(:,1), nodes(:,2), sqrt(stateR.^2 + stateI.^2));
shading interp;
view(0,90)
title('Uncontrolled State Magnitude','fontsize',24)
colorbar
axis equal
axis tight
