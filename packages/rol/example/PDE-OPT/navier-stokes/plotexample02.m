adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

filenames = {'mean_state.txt', 'state_m1_m1.txt', 'state_m1_p1.txt', 'state_p1_m1.txt', 'state_p1_p1.txt'};

for i=1:5

statefile = filenames{i};

data_obj = importdata(statefile, ' ', 2);  %% we need to skip the first two lines
state = data_obj.data;

N = 3*length(nodes);
% Extract x-velocity
Ux = state(1:3:N);
% Extract y-velocity
Uy = state(2:3:N);
% Extract pressure
P  = state(3:3:N);

figure,
trisurf(adj, nodes(:,1), nodes(:,2), P);
shading interp;
view(0,90)
axis equal
axis tight

figure,
quiver(nodes(:,1), nodes(:,2), Ux, Uy);
axis equal
axis tight

end

data_obj = importdata('control.txt', ' ', 2);  %% we need to skip the first two lines
control = data_obj.data;

% Extract x-velocity
Zx = control(1:3:N);
% Extract y-velocity
Zy = control(2:3:N);
% Extract pressure
P  = control(3:3:N);

figure,
idx = find((nodes(:,1) == 1) .* (nodes(:,2) <= 0.5));
ZX = zeros(size(Zx)); ZX(idx) = Zx(idx);
ZY = zeros(size(Zy)); ZY(idx) = Zy(idx);
quiver(nodes(:,1), nodes(:,2), ZX, ZY);
axis equal
axis tight
xlim([0.6,1.2]);
ylim([-0.2,0.7]);

obj_samples = load('obj_samples.txt');

Nint = 100;
uval = obj_samples(:,3);
zval = obj_samples(:,4);
t = linspace(min(uval), max(uval), Nint);
dist = [];
for i=1:Nint
  val = sum(uval<=t(i))/length(uval);
  dist = [dist; val];
end
figure
plot(t, dist, 'LineWidth', 3)

