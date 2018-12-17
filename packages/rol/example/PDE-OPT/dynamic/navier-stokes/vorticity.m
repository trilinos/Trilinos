function Cz = vorticity(nodes, Ux, Uy)
[nx,nt] = size(Ux);

x = linspace(-30,30,500).';
y = linspace(-15,15,250).';
[X,Y] = meshgrid(x,y);
Cz = zeros(nx,nt);
for i=1:nt
  Fx = scatteredInterpolant(nodes(:,1), nodes(:,2), Ux(:,i), 'natural');
  Fy = scatteredInterpolant(nodes(:,1), nodes(:,2), Uy(:,i), 'natural');
  Vx = Fx(X,Y);
  Vy = Fy(X,Y);
  [tmp,~] = curl(X, Y, Vx, Vy);
  Cz(:,i) = interp2(X, Y, tmp, nodes(:,1), nodes(:,2));
end
