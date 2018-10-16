
function plotFlow(name)

adj   = load('cell_to_node_quad.txt') + 1;       %% load node adjacency, increment by 1 for 1-based indexing
nodes = load('nodes.txt');                       %% load node coordinates
N     = 3*length(nodes);                         %% determine number of nodes

%% Read in state
state  = importdata([name,'.txt']);
[M,nt] = size(state);
map    = importdata(['map_',name,'.txt']);
map    = map(1:2:end)+1;
[tmp, perm] = sort(map);
state  = state(perm,:);  %% we need to permute the state according to parallel maps

Ux     = state(1:3:N,:);
Uy     = state(2:3:N,:);
P      = state(3:3:N,:);
Um     = sqrt(Ux.^2 + Uy.^2);

minUx  = min(min(Ux)); maxUx = max(max(Ux));
minUy  = min(min(Uy)); maxUy = max(max(Uy));
minUm  = min(min(Um)); maxUm = max(max(Um));
minP   = min(min(P));  maxP  = max(max(P));

axsize = 400;
figure('Position', [100 100 4*axsize 4*axsize]);
for i=1:nt
  subplot(3,2,[1 2])
  quiver(nodes(:,1), nodes(:,2), Ux(:,i), Uy(:,i));
  axis('equal','tight');
  xlabel('x');
  ylabel('y');
  title(['Time/T = ',num2str((i-1)/(nt-1))],'fontsize',16);
  set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');% tightfig;

  subplot(3,2,3)
  trisurf(adj, nodes(:,1), nodes(:,2), Ux(:,i));
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
  trisurf(adj, nodes(:,1), nodes(:,2), Um(:,i));
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
  trisurf(adj, nodes(:,1), nodes(:,2), Uy(:,i));
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
  trisurf(adj, nodes(:,1), nodes(:,2), P(:,i));
  shading interp;
  view(0,90)
  axis('equal','tight');
  xlabel('x');
  ylabel('y');
  caxis([minP,maxP])
  colormap(bone)
  title('Pressure','fontsize',16);
  set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');% tightfig;

  drawnow
end
