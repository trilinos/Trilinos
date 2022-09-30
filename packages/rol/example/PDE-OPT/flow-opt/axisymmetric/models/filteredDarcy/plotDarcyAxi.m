function plotDarcyAxi(saveplot)
if nargin==0, saveplot = false; end

addpath('~/software/export_fig')

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates
N = length(nodes);

fprintf("\nNumber of nodes: %d\n", N);
fprintf("Min r: %f\n", min(nodes(:,1)));
fprintf("Max r: %f\n", max(nodes(:,1)));
fprintf("Min z: %f\n", min(nodes(:,2)));
fprintf("Max z: %f\n", max(nodes(:,2)));

minr = min(nodes(:,1));
maxr = max(nodes(:,1));
minz = min(nodes(:,2));
maxz = max(nodes(:,2));

%% LOAD VELOCITY CONSTANT
cc = load('target.txt');
c_const = cc(1);
fprintf("\nOptimized c scalar: %f\n", c_const);

%% PLOT PRESSURE
data_obj = importdata('state.txt', ' ', 2);      %% we need to skip the first two lines
state = data_obj.data;
data_obj = importdata('map_state.txt', ' ', 9);  %% we need to skip the first 9 lines
map_state = data_obj.data;
map_state = map_state(1:2:end)+1;
[tmp, state_permute] = sort(map_state);
state = state(state_permute);                    %% we need to permute the state according to parallel maps

figure,
trisurf(adj, 1e-1*nodes(:,1), 1e-1*nodes(:,2), 1e3*state(1:N));
shading interp;
title('Pressure (Pa)')
xlabel('r (cm)')
ylabel('z (cm)')
view(0,90)
axis equal
axis tight
grid on
grid minor
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

UXinterp = scatteredInterpolant(pts(:,1),pts(:,2),vel(:,1),'natural','nearest');
UYinterp = scatteredInterpolant(pts(:,1),pts(:,2),vel(:,2),'natural','nearest');
%UXinterp = @(x,y) griddata(pts(:,1),pts(:,2),vel(:,1),x,y,'cubic');
%UYinterp = @(x,y) griddata(pts(:,1),pts(:,2),vel(:,2),x,y,'cubic');
UX = UXinterp(nodes(:,1),nodes(:,2));
UY = UYinterp(nodes(:,1),nodes(:,2));

figure,
trisurf(adj,1e-1*nodes(:,1), 1e-1*nodes(:,2), 1e-1*UX)
shading interp;
title('r-Velocity (cm/s)')
xlabel('r (cm)')
ylabel('z (cm)')
view(0,90)
axis equal
axis tight
grid on
grid minor
%set(gca,'XTickLabelMode','manual','XTickLabel',[])
%set(gca,'YTickLabelMode','manual','YTickLabel',[])
%set(gca,'Visible','off')
set(gcf,'color','white')
colormap(jet), colorbar
if saveplot, export_fig('velocity_x','-png','-painters','-m3','-a3'), end;

figure,
trisurf(adj,1e-1*nodes(:,1), 1e-1*nodes(:,2), 1e-1*UY)
shading interp;
title('z-Velocity (cm/s)')
xlabel('r (cm)')
ylabel('z (cm)')
view(0,90)
axis equal
axis tight
grid on
grid minor
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
grid on
grid minor
%set(gca,'XTickLabelMode','manual','XTickLabel',[])
%set(gca,'YTickLabelMode','manual','YTickLabel',[])
%set(gca,'Visible','off')
set(gcf,'color','white')
colormap(jet), colorbar
if saveplot, export_fig('logspeed','-png','-painters','-m3','-a3'), end;

figure,
trisurf(adj,1e-1*nodes(:,1), 1e-1*nodes(:,2), (1e-1*sqrt(UX.^2+UY.^2)))
shading interp;
title('Speed (cm/s)')
xlabel('r (cm)')
ylabel('z (cm)')
view(0,90)
axis equal
axis tight
grid on
grid minor
%set(gca,'XTickLabelMode','manual','XTickLabel',[])
%set(gca,'YTickLabelMode','manual','YTickLabel',[])
%set(gca,'Visible','off')
set(gcf,'color','white')
colormap(jet), colorbar
ylimvec = get(gca, 'ylim');
if saveplot, export_fig('speed','-png','-painters','-m3','-a3'), end;

figure,
speedvec = 1e-1*sqrt(UX.^2+UY.^2);
speedvec = (1./(speedvec<5)).*speedvec;
trisurf(adj,1e-1*nodes(:,1), 1e-1*nodes(:,2), speedvec)
shading interp;
title('Speed (cm/s)')
xlabel('r (cm)')
ylabel('z (cm)')
view(0,90)
axis equal
axis tight
grid on
grid minor
%set(gca,'XTickLabelMode','manual','XTickLabel',[])
%set(gca,'YTickLabelMode','manual','YTickLabel',[])
%set(gca,'Visible','off')
set(gcf,'color','white')
colormap(jet), colorbar
set(gca, 'ylim', ylimvec)
if saveplot, export_fig('speed5','-png','-painters','-m3','-a3'), end;

figure,
speedvec = 1e-1*sqrt(UX.^2+UY.^2);
speedvec = (1./(speedvec<2)).*speedvec;
trisurf(adj,1e-1*nodes(:,1), 1e-1*nodes(:,2), speedvec)
shading interp;
title('Speed (cm/s)')
xlabel('r (cm)')
ylabel('z (cm)')
view(0,90)
axis equal
axis tight
grid on
grid minor
%set(gca,'XTickLabelMode','manual','XTickLabel',[])
%set(gca,'YTickLabelMode','manual','YTickLabel',[])
%set(gca,'Visible','off')
set(gcf,'color','white')
colormap(jet), colorbar
set(gca, 'ylim', ylimvec)
if saveplot, export_fig('speed2','-png','-painters','-m3','-a3'), end;

figure,
speedvec = 1e-1*sqrt(UX.^2+UY.^2);
speedmask = speedvec < 5;
quiver(1e-1*nodes(:,1), 1e-1*nodes(:,2), 1e-1*UX.*speedmask, 1e-1*UY.*speedmask);
title('Velocity Quiver Plot')
xlabel('r (cm)')
ylabel('z (cm)')
axis equal
axis tight
grid on
grid minor
%set(gca,'XTickLabelMode','manual','XTickLabel',[])
%set(gca,'YTickLabelMode','manual','YTickLabel',[])
%set(gca,'Visible','off')
set(gcf,'color','white')
if saveplot, print('-depsc2','velocity5.eps'), end

pause(0.5)

M = 100;
x = linspace(0.01,10,M+1).';
y = zeros(M+1,1);
UX = UXinterp(x,y);
UY = UYinterp(x,y);
xstart = 1;

figure, plot(1e-1*x(xstart:end),1e-1*UX(xstart:end),'linewidth',2)
title('r-Velocity at Midplane (z=0)')
xlabel('r (cm)')
ylabel('Velocity (cm/s)')
axis('tight')
grid on
grid minor
set(gcf,'color','white')
pbaspect([3 1 1])
if saveplot, print('-depsc2','velocity_x_out.eps'), end

pause(0.5)

figure, plot(1e-1*x(xstart:end),1e-1*UY(xstart:end),'linewidth',2)
title('z-Velocity at Midplane (z=0)')
xlabel('r (cm)')
ylabel('Velocity (cm/s)')
axis('tight')
grid on
grid minor
set(gcf,'color','white')
pbaspect([3 1 1])
if saveplot, print('-depsc2','velocity_y_out.eps'), end

pause(0.5)

%% PLOT PERMEABILITY

%%%% PLOT PERMEABILITY AT QUADRATURE POINTS
%data = importdata('permeability.txt');
%[npts,dim] = size(data);
%dim  = dim-1;
%pts  = data(:,1:dim);
%per  = data(:,dim+1:end);
%Kinterp = scatteredInterpolant(pts(:,1),pts(:,2),per,'natural');
%K = Kinterp(nodes(:,1),nodes(:,2));

%%%% PLOT FILTERED PERMEABILITY AT NODES
data_obj = importdata('filtered_control.txt', ' ', 2);      %% we need to skip the first two lines
filtered_control = data_obj.data;
data_obj = importdata('map_filtered_control.txt', ' ', 9);  %% we need to skip the first 9 lines
map_filter = data_obj.data;
map_filter = map_filter(1:2:end)+1;
[tmp, filter_permute] = sort(map_filter);
filtered_control = filtered_control(filter_permute);        %% we need to permute the state according to parallel maps

%% the bounds of 3e-7 and 3e-6 aree hard-coded
K = 3e-7 + filtered_control*(3e-6-3e-7);
outK = [nodes(:,1) nodes(:,2) K];
save('NodalPermeability.txt', 'outK', '-ascii', '-double');

figure,
trisurf(adj,1e-1*nodes(:,1), 1e-1*nodes(:,2), 1e-6*K)
shading interp;
title('Permeability (m^2)')
xlabel('r (cm)')
ylabel('z (cm)')
view(0,90)
axis equal
axis tight
grid on
grid minor
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
fprintf('  MinPerm = %9.8e    MaxPerm = %9.8e\n',1e-6*min(K),1e-6*max(K));
fprintf('\n');


%% PLOT STREAMLINES AND TRANSIT TIMES

%%%% DEFINE VELOCITY FUNCTION
%vel = @(t,p) myVel(t, p, maxr, maxz, UXinterp, UYinterp);
%vel = @(t,p) [UXinterp(p(1),p(2)); UYinterp(p(1),p(2))];
vel = @(p) [UXinterp(p(1),p(2)); UYinterp(p(1),p(2))];

Rfunc = @(z) Rofz(z, maxr, maxz);

overFcn = @(T, Y) myEventFcn(T, Y, maxz);

options = odeset('Events',overFcn,'MaxStep',1e-2);

tspan = [0,1];
start_y = minz+0e-2;
end_x = 0.99*Rfunc(start_y);
dzmax = (maxz-minz)/10000;
chopVel = false;
chopDom = false;

fprintf('\nGenerating streamlines ');

p0  = [end_x-0.5;start_y];
%[t1,P1] = ode23s(vel,tspan,p0,options);
[P1,t1] = solveODE(p0(1),vel,Rfunc,start_y,maxz,dzmax,chopVel,chopDom);
tt1 = t1(end);
fprintf('-');

p0  = [end_x-0.4;start_y];
%[t2,P2] = ode23s(vel,tspan,p0,options);
[P2,t2] = solveODE(p0(1),vel,Rfunc,start_y,maxz,dzmax,chopVel,chopDom);
tt2 = t2(end);
fprintf('-');

p0  = [end_x-0.3;start_y];
%[t3,P3] = ode23s(vel,tspan,p0,options);
[P3,t3] = solveODE(p0(1),vel,Rfunc,start_y,maxz,dzmax,chopVel,chopDom);
tt3 = t3(end);
fprintf('-');

p0  = [end_x-0.2;start_y];
%[t4,P4] = ode23s(vel,tspan,p0,options);
[P4,t4] = solveODE(p0(1),vel,Rfunc,start_y,maxz,dzmax,chopVel,chopDom);
tt4 = t4(end);
fprintf('-');

p0  = [end_x-0.1;start_y];
%[t5,P5] = ode23s(vel,tspan,p0,options);
[P5,t5] = solveODE(p0(1),vel,Rfunc,start_y,maxz,dzmax,chopVel,chopDom);
tt5 = t5(end);
fprintf('-');

p0  = [end_x;start_y];
%[t6,P6] = ode23s(vel,tspan,p0,options);
[P6,t6] = solveODE(p0(1),vel,Rfunc,start_y,maxz,dzmax,chopVel,chopDom);
tt6 = t6(end);
fprintf('-\n');

format long
fprintf('\nOptimized streamline transit times:\n');
[tt1 tt2 tt3 tt4 tt5 tt6]

figure
hold on
plot(P1(:,1),P1(:,2),'linewidth',3)
plot(P2(:,1),P2(:,2),'linewidth',3)
plot(P3(:,1),P3(:,2),'linewidth',3)
plot(P4(:,1),P4(:,2),'linewidth',3)
plot(P5(:,1),P5(:,2),'linewidth',3)
plot(P6(:,1),P6(:,2),'linewidth',3)

%%%% TARGET STREAMLINES
%vel = @(t,p) [targetUX(p(1),p(2),maxr,maxz,c_const); targetUY(p(1),p(2),maxr,maxz,c_const)];
vel = @(p) [targetUX(p(1),p(2),maxr,maxz,c_const); targetUY(p(1),p(2),maxr,maxz,c_const)];


fprintf('\nGenerating target streamlines ');

p0  = [end_x-0.5;start_y];
%[t1,P1] = ode23s(vel,tspan,p0,options);
[P1,t1] = solveODE(p0(1),vel,Rfunc,start_y,maxz,dzmax,chopVel,chopDom);
tt1 = t1(end);
fprintf('-');

p0  = [end_x-0.4;start_y];
%[t2,P2] = ode23s(vel,tspan,p0,options);
[P2,t2] = solveODE(p0(1),vel,Rfunc,start_y,maxz,dzmax,chopVel,chopDom);
tt2 = t2(end);
fprintf('-');

p0  = [end_x-0.3;start_y];
%[t3,P3] = ode23s(vel,tspan,p0,options);
[P3,t3] = solveODE(p0(1),vel,Rfunc,start_y,maxz,dzmax,chopVel,chopDom);
tt3 = t3(end);
fprintf('-');

p0  = [end_x-0.2;start_y];
%[t4,P4] = ode23s(vel,tspan,p0,options);
[P4,t4] = solveODE(p0(1),vel,Rfunc,start_y,maxz,dzmax,chopVel,chopDom);
tt4 = t4(end);
fprintf('-');

p0  = [end_x-0.1;start_y];
%[t5,P5] = ode23s(vel,tspan,p0,options);
[P5,t5] = solveODE(p0(1),vel,Rfunc,start_y,maxz,dzmax,chopVel,chopDom);
tt5 = t5(end);
fprintf('-');

p0  = [end_x;start_y];
%[t6,P6] = ode23s(vel,tspan,p0,options);
[P6,t6] = solveODE(p0(1),vel,Rfunc,start_y,maxz,dzmax,chopVel,chopDom);
tt6 = t6(end);
fprintf('-\n');

fprintf('\nTarget streamline transit times:\n');
format long
[tt1 tt2 tt3 tt4 tt5 tt6]

plot(P1(:,1),P1(:,2),'linewidth',2,'color','black','linestyle',':')
plot(P2(:,1),P2(:,2),'linewidth',2,'color','black','linestyle',':')
plot(P3(:,1),P3(:,2),'linewidth',2,'color','black','linestyle',':')
plot(P4(:,1),P4(:,2),'linewidth',2,'color','black','linestyle',':')
plot(P5(:,1),P5(:,2),'linewidth',2,'color','black','linestyle',':')
plot(P6(:,1),P6(:,2),'linewidth',2,'color','black','linestyle',':')
hold off

axis equal
h = gca;
h.XTickMode = 'manual';
h.XTickLabel = h.XTick / 10;
h.YTickMode = 'manual';
h.YTickLabel = h.YTick / 10;
xlim([0,max(max(P6(:,1)),maxr)])
ylim([minz,maxz])
xlabel 'r (cm)'
ylabel 'z (cm)'
title 'Streamlines'
set(gcf,'color','white')
grid on
grid minor
if saveplot, export_fig('streamlines','-png','-painters','-m3','-a3'), end;
if saveplot, export_fig('streamlines','-eps','-painters'), end;
%print('-depsc2','streamlines.eps')

pause(0.5)

%vel = @(t,p) myVel(t, p, maxr, maxz, UXinterp, UYinterp);
%vel = @(t,p) [UXinterp(p(1),p(2)); UYinterp(p(1),p(2))];
vel = @(p) [UXinterp(p(1),p(2)); UYinterp(p(1),p(2))];

nPts = 100;
X = linspace(0,end_x,nPts).';

fprintf('\nGenerating transit times ');
T = zeros(nPts,1);
for i = 1:nPts
  p0 = [X(i),start_y];
  %[t,P] = ode23s(vel,tspan,p0,options);
  [P,t] = solveODE(p0(1),vel,Rfunc,start_y,maxz,dzmax,chopVel,chopDom);
  T(i) = t(end);
  fprintf('-')
end
fprintf('\n')

figure
plot(X,smoothdata(T,'rloess'),'linewidth',3)
title 'Transit time for inlet particles'
xlabel 'Particle position at inlet (mm)'
ylabel 'Time (s)'
ylim([0.4,1.6])
f = gcf;
f.Position = [120 120 700 300];
set(gcf,'color','white')
grid on
grid minor
print('-depsc2','transitTimeFullz.eps')


figure
hold on
plot(X,smoothdata(T,'rloess'),'linewidth',3)

%vel = @(t,p) [targetUX(p(1),p(2),maxr,maxz,c_const); targetUY(p(1),p(2),maxr,maxz,c_const)];
vel = @(p) [targetUX(p(1),p(2),maxr,maxz,c_const); targetUY(p(1),p(2),maxr,maxz,c_const)];

fprintf('\nGenerating target transit times ');
T = zeros(nPts,1);
for i = 1:nPts
  p0 = [X(i),start_y];
  %[t,P] = ode23s(vel,tspan,p0,options);
  [P,t] = solveODE(p0(1),vel,Rfunc,start_y,maxz,dzmax,chopVel,chopDom);
  T(i) = t(end);
  fprintf('-')
end
fprintf('\n')

plot(X,smoothdata(T,'rloess'),'linewidth',1,'color','black')
hold off

title 'Transit time for inlet particles'
xlabel 'Particle position at inlet (mm)'
ylabel 'Time (s)'
f = gcf;
f.Position = [100 100 700 300];
set(gcf,'color','white')
grid on
grid minor
print('-depsc2','transitTime.eps')

rmpath('~/software/export_fig')
end

function [value,isterminal,direction] = myEventFcn(t,y,top)
  %value      = y(2)-9.976339;
  %value      = y(2)-14;
  value      = y(2) - (top - 0e-2);
  isterminal = 1;
  direction  = 0;
end

function [UX] = targetUX(r, z, rmax, zmax, c)
  if (zmax < 10.5)
    UX = c * ( -r.*z./(rmax^2-z.^2).^2 );
  else
    UX = zeros(size(r));
    rvec = zeros(6,1); zvec = zeros(6,1);
    rvec(1) = 0.6875; zvec(1) = -14.001;
    rvec(2) = 3.0;    zvec(2) = -9.0;
    rvec(3) = 10.0;   zvec(3) = -3.0;
    rvec(4) = 10.0;   zvec(4) = 3.0;
    rvec(5) = 3.0;    zvec(5) = 9.0;
    rvec(6) = 0.6875; zvec(6) = 14.001;
    for i=1:5
      zmask = double(((z >= zvec(i)) & (z < zvec(i+1))));
      slope = (rvec(i+1) - rvec(i)) / (zvec(i+1) - zvec(i));
      UX    = UX + c * slope*zmask.*r./(rvec(i) + slope*(z-zvec(i))).^3;
    end
  end
end

function [UY] = targetUY(r, z, rmax, zmax, c)
  if (zmax < 10.5)
    UY = c * ( 1./(rmax^2-z.^2) );
  else
    UY = zeros(size(r));
    rvec = zeros(6,1); zvec = zeros(6,1);
    rvec(1) = 0.6875; zvec(1) = -14.001;
    rvec(2) = 3.0;    zvec(2) = -9.0;
    rvec(3) = 10.0;   zvec(3) = -3.0;
    rvec(4) = 10.0;   zvec(4) = 3.0;
    rvec(5) = 3.0;    zvec(5) = 9.0;
    rvec(6) = 0.6875; zvec(6) = 14.001;
    for i=1:5
      zmask = double(((z >= zvec(i)) & (z < zvec(i+1))));
      slope = (rvec(i+1) - rvec(i)) / (zvec(i+1) - zvec(i));
      UY    = UY + c * zmask./(rvec(i) + slope*(z-zvec(i))).^2;
    end
  end
end

function [R] = Rofz(z, rmax, zmax)
  if (zmax < 10.5)
    R = sqrt(rmax^2-z.^2);
  else
    R = zeros(size(z));
    rvec = zeros(6,1); zvec = zeros(6,1);
    rvec(1) = 0.6875; zvec(1) = -14.001;
    rvec(2) = 3.0;    zvec(2) = -9.0;
    rvec(3) = 10.0;   zvec(3) = -3.0;
    rvec(4) = 10.0;   zvec(4) = 3.0;
    rvec(5) = 3.0;    zvec(5) = 9.0;
    rvec(6) = 0.6875; zvec(6) = 14.001;
    for i=1:5
      zmask = double(((z >= zvec(i)) & (z < zvec(i+1))));
      slope = (rvec(i+1) - rvec(i)) / (zvec(i+1) - zvec(i));
      R     = R + zmask.*(rvec(i) + slope*(z-zvec(i)));
    end
  end
end

function [vel] = myVel(t, p, rmax, zmax, UXinterp, UYinterp)
  mask = (p(1) < Rofz(p(2),rmax,zmax)-1e-3) & (p(2) < zmax-1e-3);
  if mask
    vel = [UXinterp(p(1), p(2)); UYinterp(p(1), p(2))];
  else
    vel = [0; 0];
  end
end

function [P,T] = solveODE(r0,vel,R,z0,z1,dzmax,chopVelocity,chopDomain)
  Pnew = [r0,z0];
  V    = evaluateVel(Pnew,vel,R,chopVelocity);
  P    = [Pnew];
  T    = [0];
  while Pnew(2) < z1
    Pold     = Pnew;
    dz       = dzmax;
    if Pold(2) + dz > z1
      dz     = z1-Pold(2);
    end
    Pnew     = [Pold(1) + dz * V(1)/V(2), Pold(2) + dz];
    if (chopDomain)
      while abs(Pnew(1)) > R(Pnew(2))
        dz     = 0.5 * dz;
        Pnew   = [Pold(1) + dz * V(1)/V(2), Pold(2) + dz];
      end
    end
    V        = evaluateVel(Pnew,vel,R,chopVelocity);
    P        = [P; Pnew];
    T        = [T; T(end) + dz / V(2)];
  end
end

function V = evaluateVel(P,vel,R,chopVelocity)
  if (chopVelocity)
    V = [0;0];
    if P(1) >= 0 && P(1) < R(P(2))
      V  = vel(P);
    elseif P(1) < 0 && P(1) > -R(P(2))
      V0 = vel([-P(1),P(2)]);
      V  = [-V0(1),V0(2)];
    end
  else
    V = vel(P);
  end
end

