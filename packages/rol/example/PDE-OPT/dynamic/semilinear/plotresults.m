
a = dir('control.*.txt');
nt = numel(a);

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

axsize = 400;

figure('Position', [100 100 2*axsize 2*axsize]);
for i=0:nt-1
  ctrl = importdata(['control.',int2str(i),'.txt'], ' ', 2);  %% we need to skip the first two lines
  trisurf(adj, nodes(:,1), nodes(:,2), ctrl.data);
  shading interp;
  view(2);
  axis equal;
  axis tight;
  title(['Time t = ',num2str(i/(nt-1))],'fontsize',16);
  xlabel('x');
  ylabel('y');
  colorbar
  caxis([0,10])
  set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');% tightfig;
  drawnow
end

figure('Position', [100 100 2*axsize 2*axsize]);
for i=0:nt-1
  state = importdata(['state.',int2str(i),'.txt'], ' ', 2);  %% we need to skip the first two lines
  trisurf(adj, nodes(:,1), nodes(:,2), state.data);
  shading interp;
  view(2);
  axis equal;
  axis tight;
  title(['Time t = ',num2str(i/(nt-1))],'fontsize',16);
  xlabel('x');
  ylabel('y');
  colorbar
  caxis([0,1.25])
  set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');% tightfig;
  drawnow
end
