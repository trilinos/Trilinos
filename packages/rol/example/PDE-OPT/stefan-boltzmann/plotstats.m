dir = {'meanValue','riskAverse'};

adj = load([char(dir(1)),'/cell_to_node_quad.txt']) + 1;  %% load node adjacency table, increment by 1 for 1-based indexing
nodes = load([char(dir(1)),'/nodes.txt']);  %% load node coordinates

for i=1:length(dir)
  data_obj = importdata([char(dir(i)),'/mean_state.txt'], ' ', 2);  %% we need to skip the first two lines
  state = data_obj.data;
  figure(i)
  trisurf(adj, nodes(:,1), nodes(:,2), state);
  shading interp;
  view(0,90)
  axis equal
  axis tight
end

control = [];
for i = 1:length(dir)
  data_obj = importdata([char(dir(i)),'/control.txt'], ' ', 2);  %% we need to skip the first two lines
  control = [control, data_obj.data];
end
figure(length(dir)+1)
ind = find(nodes(:,2)==0);
plot(nodes(ind,1),control(ind,:),'linewidth',3)
legend(char(dir));
axis square

data = [];
for i = 1:length(dir)
  A    = load([char(dir(i)),'/obj_samples.txt']);
  data = [data, A(:,end-1)];
end

M = 1000;
x = linspace(min(min(data))-50,max(max(data))+50,M).';
cdf = zeros(M,length(dir));
for i = 1:length(dir)
  for j = 1:M
    ind = find(data(:,i) <= x(j));
    cdf(j,i) = length(ind)/length(data(:,i));
  end
end
figure(length(dir)+2)
plot(x,cdf,'linewidth',3)
legend(char(dir))
axis square
