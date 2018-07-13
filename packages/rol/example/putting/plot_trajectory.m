
%% Read in trajectory
T = load('position.txt');
x = T(:,1);
y = T(:,2);
z = T(:,3);

%% Build mesh grid for plotting green surface
t = linspace(-4,4,100).';
[X,Y] = meshgrid(t,t);
Z = -0.3*atan(Y) + 0.05*(X+Y);

%% Plot trajectory on green
figure, hold on
% Putting Surface
  h = surf(X,Y,Z);
  % Apply shading to green surface
  shading interp
  view(355,20)
  lightangle(5,45)
  h.FaceLighting = 'gouraud';
  h.AmbientStrength = 0.3;
  h.DiffuseStrength = 0.8;
  h.SpecularStrength = 0.9;
  h.SpecularExponent = 25;
  h.BackFaceLighting = 'unlit';
  h.FaceColor = [0.2,1,0.2];
% Pin and flag
  nx = 0.05;                       % Surface normal x-component
  ny = -0.3./(1+y(end).^2) + 0.05; % Surface normal y-component
  nz = 1;                          % Surface normal z-component
  % Pin end points
  xpin = [x(end), x(end)+nx];
  ypin = [y(end), y(end)+ny];
  zpin = [z(end), z(end)+nz];
  % Flag vertices
  xflag = [x(end)+  5*nx, x(end)+    nx, x(end)+nx, x(end)+5*nx];
  yflag = [y(end)+  5*ny, y(end)+    ny, y(end)+ny, y(end)+5*ny];
  zflag = [z(end)+0.8*nz, z(end)+0.8*nz, z(end)+nz, z(end)+  nz];
  plot3(xpin,ypin,zpin,'Color',[0.7 0.7 0.7],'LineWidth',1.5)
  fill3(xflag,yflag,zflag,[1 0 0])
% Ball
  plot3(x(1),y(1),z(1),'.w','MarkerSize',15)
% Hole
  plot3(x(end),y(end),z(end),'.k','MarkerSize',30)
% Ball Trajectory
  plot3(x,y,z,'c','LineWidth',2)
hold off

%% Set axes and add perspective
axis('tight','square','equal');
set(gca,'Projection','Perspective');
