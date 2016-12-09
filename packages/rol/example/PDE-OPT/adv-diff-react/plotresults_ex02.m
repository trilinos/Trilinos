
adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

axsize = 400;

figure('Position', [100 100 2*axsize 2*axsize])
data_obj = importdata('mean_state.txt', ' ', 2);  %% we need to skip the first two lines
state = data_obj.data;
trisurf(adj, nodes(:,1), nodes(:,2), state);
shading interp;
view(2);
axis square;
xlabel('x');
ylabel('y');
colorbar
set(gca,'FontSize',16); set(gcf, 'Color', 'White');% tightfig;

figure('Position', [100 100 2*axsize 2*axsize]);
control = load('control.txt');
data = control(1)*exp(-0.5*( (nodes(:,1)-0.25).^2./(0.05)^2 + (nodes(:,2)-0.25).^2./(0.05)^2)) + ...
       control(2)*exp(-0.5*( (nodes(:,1)-0.50).^2./(0.05)^2 + (nodes(:,2)-0.25).^2./(0.05)^2)) + ...
       control(3)*exp(-0.5*( (nodes(:,1)-0.75).^2./(0.05)^2 + (nodes(:,2)-0.25).^2./(0.05)^2)) + ...
       control(4)*exp(-0.5*( (nodes(:,1)-0.25).^2./(0.05)^2 + (nodes(:,2)-0.50).^2./(0.05)^2)) + ...
       control(5)*exp(-0.5*( (nodes(:,1)-0.50).^2./(0.05)^2 + (nodes(:,2)-0.50).^2./(0.05)^2)) + ...
       control(6)*exp(-0.5*( (nodes(:,1)-0.75).^2./(0.05)^2 + (nodes(:,2)-0.50).^2./(0.05)^2)) + ...
       control(7)*exp(-0.5*( (nodes(:,1)-0.25).^2./(0.05)^2 + (nodes(:,2)-0.75).^2./(0.05)^2)) + ...
       control(8)*exp(-0.5*( (nodes(:,1)-0.50).^2./(0.05)^2 + (nodes(:,2)-0.75).^2./(0.05)^2)) + ...
       control(9)*exp(-0.5*( (nodes(:,1)-0.75).^2./(0.05)^2 + (nodes(:,2)-0.75).^2./(0.05)^2));
trisurf(adj, nodes(:,1), nodes(:,2), data);
shading interp;
view(2);
axis square;
xlabel('x');
ylabel('y');
colorbar
set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');% tightfig;
