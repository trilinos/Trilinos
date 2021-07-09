function plotresultsPEBBL(nx,ny)

control = load([int2str(nx),'x',int2str(ny),'.sol.txt']).';
X       = load(['X_',int2str(nx),'_',int2str(ny),'.txt']);
Y       = load(['Y_',int2str(nx),'_',int2str(ny),'.txt']);
tmp = control;
for i = 1:nx
  for j = 1:ny
    control(j+(i-1)*ny) = tmp(i+(j-1)*nx);
  end
end
figure,
scatter3(2*X/nx,Y/ny,control,100,control,'filled')
view(2)
axis equal;
axis tight;
box on;
xlabel('x');
ylabel('y');
title('Computed controls');
colorbar
print('-depsc2',['control_',int2str(nx),'x',int2str(ny),'.eps']);

end
