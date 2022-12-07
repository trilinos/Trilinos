function plotDarcy(saveplot)
if nargin==0, saveplot = false; end

addpath('~/research/export_fig-master')

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates
N = length(nodes);

%% PLOT PRESSURE
data_obj = importdata('state.txt', ' ', 2);  %% we need to skip the first two lines
state = data_obj.data;
data_obj = importdata('map_state.txt', ' ', 9);  %% we need to skip the first 9 lines
map_state = data_obj.data;
map_state = map_state(1:2:end)+1;
[tmp, state_permute] = sort(map_state);
state = state(state_permute);  %% we need to permute the state according to parallel maps

figure,
trisurf(adj, 1e-1*nodes(:,1), 1e-1*nodes(:,2), 1e3*state(1:N));
shading interp;
title('Pressure (Pa)')
xlabel('r (cm)')
ylabel('z (cm)')
view(0,90)
axis equal
axis tight
%set(gca,'XTickLabelMode','manual','XTickLabel',[])
%set(gca,'YTickLabelMode','manual','YTickLabel',[])
%set(gca,'Visible','off')
set(gcf,'color','white')
colormap(jet), colorbar
if saveplot, export_fig('pressure','-png','-painters','-m3','-a3'), end;

%% PLOT VELOCITY
data = importdata('velocity.txt');
[npts,dim] = size(data);
dim  = dim/2;
pts  = data(:,1:dim);
vel  = data(:,dim+1:end);
mag  = sqrt(vel(:,1).^2+vel(:,2).^2);

UXinterp = scatteredInterpolant(pts(:,1),pts(:,2),vel(:,1),'natural');
UYinterp = scatteredInterpolant(pts(:,1),pts(:,2),vel(:,2),'natural');
UX = UXinterp(nodes(:,1),nodes(:,2));
UY = UYinterp(nodes(:,1),nodes(:,2));

figure,
trisurf(adj,1e-1*nodes(:,1), 1e-1*nodes(:,2), 1e-1*UX)
shading interp;
title('X-Velocity (cm/s)')
xlabel('r (cm)')
ylabel('z (cm)')
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
trisurf(adj,1e-1*nodes(:,1), 1e-1*nodes(:,2), 1e-1*UY)
shading interp;
title('Y-Velocity (cm/s)')
xlabel('r (cm)')
ylabel('z (cm)')
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
trisurf(adj,1e-1*nodes(:,1), 1e-1*nodes(:,2), log10(1e-1*sqrt(UX.^2+UY.^2)))
shading interp;
title('Log-Speed log(cm/s)')
xlabel('r (cm)')
ylabel('z (cm)')
view(0,90)
axis equal
axis tight
%set(gca,'XTickLabelMode','manual','XTickLabel',[])
%set(gca,'YTickLabelMode','manual','YTickLabel',[])
%set(gca,'Visible','off')
set(gcf,'color','white')
colormap(jet), colorbar
if saveplot, export_fig('speed','-png','-painters','-m3','-a3'), end;

%xlim([0,0.08])
%ylim([0.92,1])
%if saveplot, export_fig('speed-zoom','-png','-painters','-m3','-a3'), end;

figure,
quiver(1e-1*nodes(:,1), 1e-1*nodes(:,2), 1e-1*UX, 1e-1*UY);
title('Velocity Quiver Plot')
xlabel('r (cm)')
ylabel('z (cm)')
axis equal
axis tight
%set(gca,'XTickLabelMode','manual','XTickLabel',[])
%set(gca,'YTickLabelMode','manual','YTickLabel',[])
%set(gca,'Visible','off')
set(gcf,'color','white')
if saveplot, print('-depsc2','velocity.eps'), end

pause(5)

M = 1000;
x = linspace(0,10,M+1).';
y = zeros(M+1,1);
UX = UXinterp(x,y);
UY = UYinterp(x,y);

figure, plot(1e-1*x(3:end),1e-1*UX(3:end),'linewidth',2)
title('Radial Velocity at Target Surface')
xlabel('r (cm)')
ylabel('Velocity (cm/s)')
axis('tight')
set(gcf,'color','white')
pbaspect([3 1 1])
if saveplot, print('-depsc2','velocity_x_out.eps'), end

pause(5)

figure, plot(1e-1*x(3:end),1e-1*UY(3:end),'linewidth',2)
title('Axial Velocity at Target Surface')
xlabel('z (cm)')
ylabel('Velocity (cm/s)')
axis('tight')
set(gcf,'color','white')
pbaspect([3 1 1])
if saveplot, print('-depsc2','velocity_y_out.eps'), end

pause(5)

%% PLOT PERMEABILITY
data = importdata('permeability.txt');
[npts,dim] = size(data);
dim  = dim-1;
pts  = data(:,1:dim);
per  = data(:,dim+1:end);

Kinterp = scatteredInterpolant(pts(:,1),pts(:,2),per,'natural');
K = Kinterp(nodes(:,1),nodes(:,2));

figure,
trisurf(adj,1e-1*nodes(:,1), 1e-1*nodes(:,2), 1e-6*K)
shading interp;
title('Permeability (m^2)')
xlabel('r (cm)')
ylabel('z (cm)')
view(0,90)
axis equal
axis tight
%set(gca,'XTickLabelMode','manual','XTickLabel',[])
%set(gca,'YTickLabelMode','manual','YTickLabel',[])
%set(gca,'Visible','off')
set(gcf,'color','white')
colormap(jet), colorbar
if saveplot, export_fig('permeability','-png','-painters','-m3','-a3'), end;

%% OUTPUT INFORMATION
fprintf('\n');
fprintf('  MinPres = %9.8e    MaxPres = %9.8e\n',1e3*min(state),1e3*max(state));
fprintf('  MinVelo = %9.8e    MaxVelo = %9.8e\n',1e-1*min(mag),1e-1*max(mag));
fprintf('  MinPerm = %9.8e    MaxPerm = %9.8e\n',1e-6*min(per),1e-6*max(per));
fprintf('\n');

%Plot transit times
vel = @(t,p) [UXinterp(p(1),p(2)); UYinterp(p(1),p(2))];

options = odeset('Events',@myEventFcn,'MaxStep',1e-1);

tspan = [0,1e1];

p0  = [0.1;-10];
[t1,P1] = ode45(vel,tspan,p0,options);
tt1 = t1(end);

p0  = [0.2;-10];
[t2,P2] = ode45(vel,tspan,p0,options);
tt2 = t2(end);

p0  = [0.3;-10];
[t3,P3] = ode45(vel,tspan,p0,options);
tt3 = t3(end);

p0  = [0.4;-10];
[t4,P4] = ode45(vel,tspan,p0,options);
tt4 = t4(end);

p0  = [0.5;-10];
[t5,P5] = ode45(vel,tspan,p0,options);
tt5 = t5(end);

p0  = [0.6;-10];
[t6,P6] = ode45(vel,tspan,p0,options);
tt6 = t6(end);

format long
[tt1 tt2 tt3 tt4 tt5 tt6]

figure,
plot(P1(:,1),P1(:,2),'linewidth',3), hold on
plot(P2(:,1),P2(:,2),'linewidth',3)
plot(P3(:,1),P3(:,2),'linewidth',3)
plot(P4(:,1),P4(:,2),'linewidth',3)
plot(P5(:,1),P5(:,2),'linewidth',3)
plot(P6(:,1),P6(:,2),'linewidth',3), hold off
axis equal
xlim([0,10])
ylim([-10,10])
print('-depsc2','streamlines.eps')

pause(0.5)

X = linspace(0,0.6875,100).';
T = zeros(100,1);
for i = 1:100
  p0 = [X(i),-10];
  [t,P] = ode45(vel,tspan,p0,options);
  T(i) = t(end);
end

figure,
plot(X,smoothdata(T,'rloess'),'linewidth',3)
%axis equal
%axis tight
print('-depsc2','transitTime.eps')

rmpath('~/research/export_fig-master')
return

function [value,isterminal,direction] = myEventFcn(t,y)
  value      = y(2)-10;
  isterminal = 1;
  direction  = 0;
return
