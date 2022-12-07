saveplot = false;

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

data_obj = importdata('state.txt', ' ', 2);  %% we need to skip the first two lines
state = data_obj.data;
data_obj = importdata('map_state.txt', ' ', 9);  %% we need to skip the first 9 lines
map_state = data_obj.data;
map_state = map_state(1:2:end)+1;
[tmp, state_permute] = sort(map_state);
state = state(state_permute);  %% we need to permute the state according to parallel maps

N = 3*length(nodes);
% Extract x-velocity
Ux = state(1:3:N);
% Extract y-velocity
Uy = state(2:3:N);
% Extract pressure
P  = state(3:3:N);

figure,
trisurf(adj, nodes(:,1), nodes(:,2), Ux);
shading interp;
title('X-Velocity')
view(0,90)
axis equal
axis tight
%set(gca,'XTickLabelMode','manual','XTickLabel',[])
%set(gca,'YTickLabelMode','manual','YTickLabel',[])
%set(gca,'Visible','off')
set(gcf,'color','white')
colormap(jet), colorbar
if saveplot, export_fig('velocity_x','-png','-painters','-m3','-a3'), end;

figure,
trisurf(adj, nodes(:,1), nodes(:,2), Uy);
shading interp;
title('Y-Velocity')
view(0,90)
axis equal
axis tight
%set(gca,'XTickLabelMode','manual','XTickLabel',[])
%set(gca,'YTickLabelMode','manual','YTickLabel',[])
%set(gca,'Visible','off')
set(gcf,'color','white')
colormap(jet), colorbar
if saveplot, export_fig('velocity_y','-png','-painters','-m3','-a3'), end;

figure,
trisurf(adj, nodes(:,1), nodes(:,2), sqrt(Ux.^2 + Uy.^2));
shading interp;
title('Speed')
view(0,90)
axis equal
axis tight
%set(gca,'XTickLabelMode','manual','XTickLabel',[])
%set(gca,'YTickLabelMode','manual','YTickLabel',[])
%set(gca,'Visible','off')
set(gcf,'color','white')
colormap(jet), colorbar
if saveplot, export_fig('speed','-png','-painters','-m3','-a3'), end;

figure,
trisurf(adj, nodes(:,1), nodes(:,2), P);
shading interp;
title('Pressure')
view(0,90)
axis equal
axis tight
%set(gca,'XTickLabelMode','manual','XTickLabel',[])
%set(gca,'YTickLabelMode','manual','YTickLabel',[])
%set(gca,'Visible','off')
set(gcf,'color','white')
colormap(jet), colorbar
if saveplot, export_fig('pressure','-png','-painters','-m3','-a3'), end;

figure,
quiver(nodes(:,1), nodes(:,2), Ux, Uy);
title('Velocity Quiver Plot')
axis equal
axis tight
%set(gca,'XTickLabelMode','manual','XTickLabel',[])
%set(gca,'YTickLabelMode','manual','YTickLabel',[])
%set(gca,'Visible','off')
set(gcf,'color','white')
if saveplot, print('-depsc2','velocity.eps'), end

ztol = 1e-8;
ind =( (nodes(:,2)>-ztol) & (nodes(:,2)<ztol) );
[X,ord] = sort(nodes(ind,1),'ascend');
UX = Ux(ind);
UX = UX(ord);
UY = Uy(ind);
UY = UY(ord);
figure, plot(X,UX,'linewidth',2)
title('X-Velocity at Target Surface')
axis('tight')
set(gcf,'color','white')
pbaspect([3 1 1])
if saveplot, print('-depsc2','velocity_x_out.eps'), end

figure, plot(X,UY,'linewidth',2)
title('Y-Velocity at Target Surface')
axis('tight')
set(gcf,'color','white')
pbaspect([3 1 1])
if saveplot, print('-depsc2','velocity_y_out.eps'), end

data_obj = importdata('control.txt', ' ', 2);  %% we need to skip the first two lines
control = data_obj.data;
data_obj = importdata('map_control.txt', ' ', 9);  %% we need to skip the first 9 lines
map_control = data_obj.data;
map_control = map_control(1:2:end)+1;
[tmp, control_permute] = sort(map_control);
control = control(control_permute);  %% we need to permute the control according to parallel maps

% Extract x-velocity
Zx = control(1:3:N);
% Extract y-velocity
Zy = control(2:3:N);
% Extract pressure
P  = control(3:3:N);

figure,
trisurf(adj, nodes(:,1), nodes(:,2), P);
shading interp;
title('Porosity')
view(0,90)
axis equal
axis tight
%set(gca,'XTickLabelMode','manual','XTickLabel',[])
%set(gca,'YTickLabelMode','manual','YTickLabel',[])
%set(gca,'Visible','off')
set(gcf,'color','white')
colormap(flipud(jet)), colorbar
if saveplot, export_fig('density','-png','-painters','-m3','-a3'), end;

a0 = 0;
a1 = 750;
q  = 10;
figure,
trisurf(adj, nodes(:,1), nodes(:,2), a0 + (a1-a0)*q*(1-P)./(q+P));
shading interp;
title('Impermeability')
view(0,90)
axis equal
axis tight
%set(gca,'XTickLabelMode','manual','XTickLabel',[])
%set(gca,'YTickLabelMode','manual','YTickLabel',[])
%set(gca,'Visible','off')
set(gcf,'color','white')
colormap(jet), colorbar
