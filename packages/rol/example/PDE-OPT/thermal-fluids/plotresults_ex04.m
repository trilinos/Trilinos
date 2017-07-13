printFig = false;

addpath ~/research/export_fig-master

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

filenames = {'mean_uncontrolled_state.txt', 'mean_state.txt'};
names = {'mean_uncontrolled_','mean_'};

for i = 1:2
  
  statefile = filenames{i};
  
  %% UNCONTROLLED.
  data_obj = importdata(statefile, ' ', 2);  %% we need to skip the first two lines
  state = data_obj.data;
  
  n = length(nodes);
  N = 4*n;
  % Extract x-velocity
  Ux = state(1:4:N);
  % Extract y-velocity
  Uy = state(2:4:N);
  % Extract pressure
  P  = state(3:4:N);
  % Extract temperature
  T  = state(4:4:N);
  
  figure,
  quiver(nodes(:,1), nodes(:,2), Ux, Uy);
  axis equal
  axis tight
  set(gca,'XTickLabelMode','manual','XTickLabel',[])
  set(gca,'YTickLabelMode','manual','YTickLabel',[])
  if (printFig)
    print('-depsc2',[names{i},'velocity.eps'])
  end

  XX = reshape(nodes(:,1),sqrt(n),sqrt(n)).';
  YY = reshape(nodes(:,2),sqrt(n),sqrt(n)).';
  UXX = reshape(Ux,sqrt(n),sqrt(n)).';
  UYY = reshape(Uy,sqrt(n),sqrt(n)).';
  PP  = reshape(P,sqrt(n),sqrt(n)).';
  TT  = reshape(T,sqrt(n),sqrt(n)).';

  figure, streamslice(XX,YY,UXX,UYY)
  axis equal
  axis tight
  set(gca,'XTickLabelMode','manual','XTickLabel',[])
  set(gca,'YTickLabelMode','manual','YTickLabel',[])
  if (printFig)
    print('-depsc2',[names{i},'velocity_vort.eps'])
  end
  
  figure,
  trisurf(adj, nodes(:,1), nodes(:,2), P);
  shading interp;
  view(0,90)
  axis equal
  axis tight
  set(gca,'XTickLabelMode','manual','XTickLabel',[])
  set(gca,'YTickLabelMode','manual','YTickLabel',[])
  set(gcf,'color','white')
  if (printFig)
    export_fig([names{i},'pressure'],'-png','-painters','-m3','-a3');
  end

  figure,
  contourf(XX,YY,PP,15)
  axis equal
  axis tight
  set(gca,'XTickLabelMode','manual','XTickLabel',[])
  set(gca,'YTickLabelMode','manual','YTickLabel',[])
  set(gcf,'color','white')
  if (printFig)
    export_fig([names{i},'pressure_contour'],'-png','-painters','-m3','-a3');
  end
  
  figure,
  trisurf(adj, nodes(:,1), nodes(:,2), T);
  shading interp;
  view(0,90)
  axis equal
  axis tight
  set(gca,'XTickLabelMode','manual','XTickLabel',[])
  set(gca,'YTickLabelMode','manual','YTickLabel',[])
  set(gcf,'color','white')
  if (printFig)
    export_fig([names{i},'temperature'],'-png','-painters','-m3','-a3');
  end

  figure,
  contourf(XX,YY,TT,15)
  axis equal
  axis tight
  set(gca,'XTickLabelMode','manual','XTickLabel',[])
  set(gca,'YTickLabelMode','manual','YTickLabel',[])
  set(gcf,'color','white')
  if (printFig)
    export_fig([names{i},'temperature_contour'],'-png','-painters','-m3','-a3');
  end
  
end

data_obj = importdata('control.txt', ' ', 2);  %% we need to skip the first two lines
control = data_obj.data;

% Extract x-velocity
Zx = control(1:4:N);
% Extract y-velocity
Zy = control(2:4:N);
% Extract pressure
P  = control(3:4:N);
% Extract Temperature
T  = control(4:4:N);

Xmin = min(nodes(:,2));
Xmax = max(nodes(:,2));

figure,
idx = find(nodes(:,1) == 0.0);
plot(nodes(idx,2), T(idx));
title('Left Controls')
xlim([Xmin,Xmax]);
if (printFig)
  print('-depsc2','left_control.txt');
end

figure,
idx = find(nodes(:,1) == 1.0);
plot(nodes(idx,2), T(idx));
title('Right Controls')
xlim([Xmin,Xmax]);
if (printFig)
  print('-depsc2','right_control.txt');
end
