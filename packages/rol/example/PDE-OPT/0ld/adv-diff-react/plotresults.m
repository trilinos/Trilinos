function plotresults

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing
nodes = load('nodes.txt');  %% load node coordinates
publish_results = 0;

axsize = 400;

figure('Position', [100 100 3.3*axsize 1.6*axsize])
localplot('state.txt', 'control.txt', 'weights.txt', 'All sensors', adj, nodes, publish_results, axsize);

if exist('weightsOED.txt', 'file') == 2
  figure('Position', [150 150 3.3*axsize 1.6*axsize])
  localplot('stateOED.txt', 'controlOED.txt', 'weightsOED.txt', 'OED sensors', adj, nodes, publish_results, axsize);
end


if exist('weightsRandom.txt', 'file') == 2
  figure('Position', [200 200 3.3*axsize 1.6*axsize])
  localplot('stateRandom.txt', 'controlRandom.txt', 'weightsRandom.txt', 'Random sensors', adj, nodes, publish_results, axsize);
end

end % plotresults


function localplot(statefile, controlfile, weightsfile, weightscaption, adj, nodes, publish_results, axsize)

subplot(2,3,4)
data_obj = importdata('data.txt', ' ', 2);  %% we need to skip the first two lines
data = data_obj.data;
trisurf(adj, nodes(:,1), nodes(:,2), data);
shading interp;
view(2);
axis square;
xlabel('x');
ylabel('y');
title('Data (true state)');

subplot(2,3,1)
data_obj = importdata('sources.txt', ' ', 2);  %% we need to skip the first two lines
sources = data_obj.data;
% normalize source data to 1
maxsource = max(abs(sources));
sources = sources * (1/maxsource);
trisurf(adj, nodes(:,1), nodes(:,2), sources);
shading interp;
view(2);
axis square;
xlabel('x');
ylabel('y');
title('True sources');

subplot(2,3,5)
data_obj = importdata(statefile, ' ', 2);  %% we need to skip the first two lines
state = data_obj.data;
trisurf(adj, nodes(:,1), nodes(:,2), state);
shading interp;
view(2);
axis square;
xlabel('x');
ylabel('y');
title('Computed state');

subplot(2,3,2)
data_obj = importdata(controlfile, ' ', 2);  %% we need to skip the first two lines
control = data_obj.data;
trisurf(adj, nodes(:,1), nodes(:,2), control);
shading interp;
view(2);
axis square;
xlabel('x');
ylabel('y');
title('Computed sources');

subplot(2,3,3)
trisurf(adj, nodes(:,1), nodes(:,2), sources-control);
shading interp;
view(2);
axis square;
xlabel('x');
ylabel('y');
title('Source error (rel)');
colorbar

subplot(2,3,6)
trisurf(adj, nodes(:,1), nodes(:,2), (state-data)/max(abs(data)));
shading interp;
view(2);
axis square;
xlabel('x');
ylabel('y');
title('State error (rel)');
colorbar

set(findall(gcf,'-property','FontSize'),'FontSize',16); set(gcf, 'Color', 'White'); tightfig;
if publish_results
  myaa(2);
  myaa('publish');
end

fig1 = gcf;
sizevec = get(fig1, 'Position');
figure('Position', [sizevec(1)+sizevec(3)+15 sizevec(2) axsize axsize]);
data_obj = importdata('data.txt', ' ', 2);  %% we need to skip the first two lines
data = data_obj.data;
trisurf(adj, nodes(:,1), nodes(:,2), data);
shading interp;
view(2);
axis square;
hold on;
data_obj = importdata(weightsfile, ' ', 2);  %% we need to skip the first two lines
weights = data_obj.data;
sum(weights)
xloc = nodes(logical(weights),1);
yloc = nodes(logical(weights),2);
scatter3(xloc, yloc, 10*ones(size(xloc)), 22, 's', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'white', 'LineWidth', 2);
view(2);
axis square;
xlabel('x');
ylabel('y');
title(weightscaption);
set(findall(gcf,'-property','FontSize'),'FontSize',16); set(gcf, 'Color', 'White'); tightfig;
if publish_results
  myaa(2);
  myaa('publish');
end

end % localplot
