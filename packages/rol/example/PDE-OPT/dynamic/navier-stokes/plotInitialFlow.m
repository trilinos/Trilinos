adj   = load('cell_to_node_quad.txt') + 1;       %% load node adjacency, increment by 1 for 1-based indexing
nodes = load('nodes.txt');                       %% load node coordinates
n     = length(nodes);
N     = 3*length(nodes);                         %% determine number of nodes

%% Read in state
state  = importdata('initial_condition_Re200.txt',' ',2);
state  = state.data;
map    = importdata('map_initial_condition_Re200.txt',' ',9);
map    = map.data;
map    = map(1:2:end)+1;
[tmp, perm] = sort(map);
state  = state(perm);  %% we need to permute the state according to parallel maps

Ux     = state(1:3:N);
Uy     = state(2:3:N);
P      = state(3:3:N);
Um     = sqrt(Ux.^2 + Uy.^2);

minUx  = min(Ux); maxUx = max(Ux);
minUy  = min(Uy); maxUy = max(Uy);
minUm  = min(Um); maxUm = max(Um);
minP   = min(P);  maxP  = max(P);

axsize = 400;
figure('Position', [100 100 4*axsize 4*axsize]);
subplot(3,2,[1 2])
quiver(nodes(:,1), nodes(:,2), Ux, Uy);
axis('equal','tight');
xlabel('x');
ylabel('y');
set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');% tightfig;

subplot(3,2,3)
trisurf(adj, nodes(:,1), nodes(:,2), Ux);
shading interp;
view(0,90)
axis('equal','tight');
xlabel('x');
ylabel('y');
caxis([minUx,maxUx])
colormap(bone)
title('x-Velocity','fontsize',16);
set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');% tightfig;

subplot(3,2,4)
trisurf(adj, nodes(:,1), nodes(:,2), Um);
shading interp;
view(0,90)
axis('equal','tight');
xlabel('x');
ylabel('y');
caxis([minUm,maxUm])
colormap(bone)
title('Velocity Magnitude','fontsize',16);
set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');% tightfig;

subplot(3,2,5)
trisurf(adj, nodes(:,1), nodes(:,2), Uy);
shading interp;
view(0,90)
axis('equal','tight');
xlabel('x');
ylabel('y');
caxis([minUy,maxUy])
colormap(bone)
title('y-Velocity','fontsize',16);
set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');% tightfig;

subplot(3,2,6)
trisurf(adj, nodes(:,1), nodes(:,2), P);
shading interp;
view(0,90)
axis('equal','tight');
xlabel('x');
ylabel('y');
caxis([minP,maxP])
colormap(bone)
title('Pressure','fontsize',16);
set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');% tightfig;

printFig = true;
if printFig
  addpath ~/research/export_fig-master
  
  figure,
  quiver(nodes(:,1), nodes(:,2), Ux, Uy, 0.5);
  axis equal
  axis tight
  xlim([-2,6])
  ylim([-2,2])
  set(gca,'XTickLabelMode','manual','XTickLabel',[])
  set(gca,'YTickLabelMode','manual','YTickLabel',[])
  set(gca,'Visible','off')
  set(gcf,'color','white')
  print('-depsc2','initial_velocity_quiver.eps')

  figure,
  trisurf(adj, nodes(:,1), nodes(:,2), Ux);
  shading interp;
  view(0,90)
  axis equal
  axis tight
  xlim([-5,15])
  ylim([-5,5])
  set(gca,'XTickLabelMode','manual','XTickLabel',[])
  set(gca,'YTickLabelMode','manual','YTickLabel',[])
  set(gca,'Visible','off')
  set(gcf,'color','white')
  colormap(bone)
  export_fig('initial_velocity_x','-png','-painters','-m3','-a3');

  figure,
  trisurf(adj, nodes(:,1), nodes(:,2), Uy);
  shading interp;
  view(0,90)
  axis equal
  axis tight
  xlim([-5,15])
  ylim([-5,5])
  set(gca,'XTickLabelMode','manual','XTickLabel',[])
  set(gca,'YTickLabelMode','manual','YTickLabel',[])
  set(gca,'Visible','off')
  set(gcf,'color','white')
  colormap(bone)
  export_fig('initial_velocity_y','-png','-painters','-m3','-a3');

  figure,
  trisurf(adj, nodes(:,1), nodes(:,2), Um);
  shading interp;
  view(0,90)
  axis equal
  axis tight
  xlim([-5,15])
  ylim([-5,5])
  set(gca,'XTickLabelMode','manual','XTickLabel',[])
  set(gca,'YTickLabelMode','manual','YTickLabel',[])
  set(gca,'Visible','off')
  set(gcf,'color','white')
  colormap(bone)
  export_fig('initial_velocity_mag','-png','-painters','-m3','-a3');

  figure,
  trisurf(adj, nodes(:,1), nodes(:,2), P);
  shading interp;
  view(0,90)
  axis equal
  axis tight
  xlim([-5,15])
  ylim([-5,5])
  set(gca,'XTickLabelMode','manual','XTickLabel',[])
  set(gca,'YTickLabelMode','manual','YTickLabel',[])
  set(gca,'Visible','off')
  set(gcf,'color','white')
  colormap(bone)
  export_fig('initial_pressure','-png','-painters','-m3','-a3');

  x = linspace(-30,30,500).';
  y = linspace(-15,15,250).';
  [X,Y] = meshgrid(x,y);
  Fx = scatteredInterpolant(nodes(:,1), nodes(:,2), Ux, 'natural');
  Fy = scatteredInterpolant(nodes(:,1), nodes(:,2), Uy, 'natural');
  Vx = Fx(X,Y);
  Vy = Fy(X,Y);
  [cz,~] = curl(X, Y, Vx, Vy);
  Cz = interp2(X, Y, cz, nodes(:,1), nodes(:,2));
  figure,
  trisurf(adj, nodes(:,1), nodes(:,2), Cz)
  %contour(X, Y, Cz, 30)
  shading interp;
  view(0,90)
  axis equal
  axis tight
  xlim([-5,15])
  ylim([-5,5])
  set(gca,'XTickLabelMode','manual','XTickLabel',[])
  set(gca,'YTickLabelMode','manual','YTickLabel',[])
  set(gca,'Visible','off')
  %box on
  set(gcf,'color','white')
  %print('-depsc2','initial_vorticity.eps')
  colormap(bone)
  export_fig('initial_vorticity','-png','-painters','-m3','-a3');
end

