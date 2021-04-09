function plotresultsPC1(nx,ny)

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

publish_results = 0;

sensor = load('sensor.txt');
figure,
scatter(sensor(:,1),sensor(:,2),30,sensor(:,3),'filled','o')
axis equal
xlim([0,2]);
ylim([0,1])

%axsize = 400;
%figure('Position', [100 100 3.3*axsize 1.6*axsize])
%
%subplot(1,2,1)
%data_obj = importdata(['state_',int2str(order),'.txt'], ' ', 2);  %% we need to skip the first two lines
%state = data_obj.data;
%trisurf(adj, nodes(:,1), nodes(:,2), state);
%hold on
%plot3(sensor(:,1),sensor(:,2),max(state)*ones(size(sensor(:,1))),'k.','MarkerSize',10)
%hold off
%shading interp;
%view(2);
%axis equal;
%axis tight;
%xlabel('x');
%ylabel('y');
%title('Computed state');
%colorbar

%subplot(1,2,2)
control = load([int2str(nx),'x',int2str(ny),'.sol.txt']).';
X       = load(['X_',int2str(nx),'_',int2str(ny),'.txt']);
Y       = load(['Y_',int2str(nx),'_',int2str(ny),'.txt']);
tmp = control;
for i = 1:nx
  for j = 1:ny
    control(j+(i-1)*ny) = tmp(i+(j-1)*nx);
  end
end
scatter3(2*X/nx,Y/ny,control,100,control,'filled')
view(2)
%imagesc(control)
axis equal;
axis tight;
box on;
xlabel('x');
ylabel('y');
title('Computed controls');
colorbar

set(findall(gcf,'-property','FontSize'),'FontSize',16); set(gcf, 'Color', 'White');% tightfig;
if publish_results
  myaa(2);
  myaa('publish');
end

end
