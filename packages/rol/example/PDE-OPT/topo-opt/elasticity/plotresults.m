function plotresults(nx,ny)

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

data_obj = importdata('density.txt', ' ', 2);  %% we need to skip the first two lines
dens = data_obj.data(1:2:end,:);
figure
patch('Faces',adj,'Vertices',nodes,'FaceVertexCData',dens,'FaceColor','interp','EdgeColor','none');
colormap(flipud(gray))
axis equal;
axis tight;
box on;
xlabel('x');,
ylabel('y');
title('Interpolated nodal density');

figure
xnodes = nodes(:,1);
ynodes = nodes(:,2);
celldens = dens(adj);
avgcelldens = sum(celldens,2)/4;
%maxcelldens = max(celldens,[],2);
patch(xnodes(adj)',ynodes(adj)',ones(size(xnodes(adj)')),'CData',avgcelldens,'FaceColor','flat','EdgeColor','none');
colormap(flipud(gray))
axis equal;
axis tight;
box on;
xlabel('x');
ylabel('y');
title('Averaged nodal density');

figure
scatter(xnodes,ynodes,40*dens,'o','filled','k')
axis equal;
axis tight;
box on;
xlabel('x');
ylabel('y');
title('Bubble plot nodal density');

figure
v1 = 0.1;
v2 = 1.0;
V  = [v1,v2];
contourf(reshape(nodes(:,1),nx+1,ny+1), reshape(nodes(:,2),nx+1,ny+1), reshape(dens,nx+1,ny+1), V);
colormap([0 0 0]);
axis equal;
axis tight;
box on
xlabel('x');
ylabel('y');
title('Contoured nodal density');

data_obj = importdata('state.txt', ' ', 2);  %% we need to skip the first two lines
state = data_obj.data;

N = 2*length(nodes);
% Extract x-displacement
Ux = state(1:2:N);
% Extract y-displacement
Uy = state(2:2:N);

% Mask, depending on density.
mask = dens > min(dens)*1e2;
UxPlot = Ux.*mask;
UyPlot = Uy.*mask;

figure,
trisurf(adj, nodes(:,1), nodes(:,2), UxPlot);
shading interp;
view(0,90)
axis equal
axis tight
box on
caxis(gca, max(abs(caxis(gca))) * [-1 1]);
newmap = bwr;
colormap(newmap)
colorbar

figure,
trisurf(adj, nodes(:,1), nodes(:,2), UyPlot);
shading interp;
view(0,90)
axis equal
axis tight
box on
caxis(gca, max(abs(caxis(gca))) * [-1 1]);
newmap = bwr;
colormap(newmap)
colorbar

figure
xnodes = nodes(:,1)+Ux;
ynodes = nodes(:,2)+Uy;
celldens = dens(adj);
avgcelldens = sum(celldens,2)/4;
%maxcelldens = max(celldens,[],2);
patch(xnodes(adj)',ynodes(adj)',ones(size(xnodes(adj)')),'CData',avgcelldens,'FaceColor','flat','EdgeColor','none');
colormap(flipud(gray))
axis equal;
axis tight;
box on;
xlabel('x');
ylabel('y');
title('Deformed density');
