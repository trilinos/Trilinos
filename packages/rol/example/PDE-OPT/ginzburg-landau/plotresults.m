
function plotresults(Nx,Ny,x0,x1,y0,y1)
Nx = Nx+1;
Ny = Ny+1;

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

figure(4)
contour(reshape(nodes(:,1),Nx,Ny), reshape(nodes(:,2),Nx,Ny), reshape(sqrt(stateR.^2 + stateI.^2),Nx,Ny));
title('State Magnitude','fontsize',24)
axis equal
axis tight

figure(5)
contour(reshape(nodes(:,1),Nx,Ny), reshape(nodes(:,2),Nx,Ny), reshape(stateR,Nx,Ny));
title('Real State','fontsize',24)
axis equal
axis tight

figure(6)
contour(reshape(nodes(:,1),Nx,Ny), reshape(nodes(:,2),Nx,Ny), reshape(stateI,Nx,Ny));
title('Imaginary State','fontsize',24)
axis equal
axis tight

data_obj = importdata('control.txt', ' ', 2);  %% we need to skip the first two lines
control = data_obj.data;

controlR = control(1:2:end);
controlI = control(2:2:end);

figure(7)
ind = find(nodes(:,2)==y0);
plot(nodes(ind,1),controlR(ind),'b','linewidth',3), hold on
plot(nodes(ind,1),controlI(ind),'r','linewidth',3), hold off
title(['Control, y=',int2str(y0)],'fontsize',24)
axis square

figure(8)
ind = find(nodes(:,2)==y1);
plot(nodes(ind,1),controlR(ind),'b','linewidth',3), hold on
plot(nodes(ind,1),controlI(ind),'r','linewidth',3), hold off
title(['Control, y=',int2str(y1)],'fontsize',24)
axis square

figure(9)
ind = find(nodes(:,1)==x0);
plot(nodes(ind,2),controlR(ind),'b','linewidth',3), hold on
plot(nodes(ind,2),controlI(ind),'r','linewidth',3), hold off
title(['Control, x=',int2str(x0)],'fontsize',24)
axis square

figure(10)
ind = find(nodes(:,1)==x1);
plot(nodes(ind,2),controlR(ind),'b','linewidth',3), hold on
plot(nodes(ind,2),controlI(ind),'r','linewidth',3), hold off
title(['Control, x=',int2str(x1)],'fontsize',24)
axis square

data_obj = importdata('state_uncontrolled.txt', ' ', 2);  %% we need to skip the first two lines
state = data_obj.data;

stateR = state(1:2:end);
stateI = state(2:2:end);

figure(11)
trisurf(adj, nodes(:,1), nodes(:,2), stateR);
shading interp;
view(0,90)
title('Real Uncontrolled State','fontsize',24)
colorbar
axis equal
axis tight

figure(12)
trisurf(adj, nodes(:,1), nodes(:,2), stateI);
shading interp;
view(0,90)
title('Imaginary Uncontrolled State','fontsize',24)
colorbar
axis equal
axis tight

figure(13)
trisurf(adj, nodes(:,1), nodes(:,2), sqrt(stateR.^2 + stateI.^2));
shading interp;
view(0,90)
title('Uncontrolled State Magnitude','fontsize',24)
colorbar
axis equal
axis tight

figure(14)
contour(reshape(nodes(:,1),Nx,Ny), reshape(nodes(:,2),Nx,Ny), reshape(sqrt(stateR.^2 + stateI.^2),Nx,Ny));
title('Uncontrolled State Magnitude','fontsize',24)
axis equal
axis tight

figure(15)
contour(reshape(nodes(:,1),Nx,Ny), reshape(nodes(:,2),Nx,Ny), reshape(stateR,Nx,Ny));
title('Real Uncontrolled State','fontsize',24)
axis equal
axis tight

figure(16)
contour(reshape(nodes(:,1),Nx,Ny), reshape(nodes(:,2),Nx,Ny), reshape(stateI,Nx,Ny));
title('Imaginary Uncontrolled State','fontsize',24)
axis equal
axis tight
