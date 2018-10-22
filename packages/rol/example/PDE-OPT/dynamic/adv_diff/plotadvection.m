
nx = 25;
nt = 200;
xmesh = linspace(0,1,nx).';
tmesh = linspace(0,1,nt).';
[X,Y] = meshgrid(xmesh);
U0 = @(t,x,y)  cos((x-t)*2*pi).*sin((y-t)*2*pi).*exp(-2*t)/(2*pi) + 5e-2*(7.5 - 2.5*x).*exp(20*(t-1));
V0 = @(t,x,y) -sin((x-t)*2*pi).*cos((y-t)*2*pi).*exp(-2*t)/(2*pi) + 5e-2*(      2.5*y).*exp(20*(t-1));

figure('Position', [100 100 2*axsize 2*axsize]);
for i = 1:nt
  U = U0(tmesh(i),X,Y);
  V = V0(tmesh(i),X,Y);
  quiver(X,Y,U,V,0,'LineWidth',1.25,'MaxHeadSize',0.75/nx)
  title(['Time t = ',num2str((i-1)/(nt-1))],'fontsize',16);
  set(gca,'fontsize',16)
  axis square
  xlim([0,1])
  ylim([0,1])
  drawnow
end
