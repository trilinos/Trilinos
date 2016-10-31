function plotdensity(filename, xsize, ysize)

data_obj = importdata(filename, ' ', 2);  %% we need to skip the first two lines
density = data_obj.data;
nodes = load('nodes.txt');

density = reshape(density, xsize, ysize);

xnodes  = reshape(nodes(:,1), xsize+1, ysize+1);
ynodes  = reshape(nodes(:,2), xsize+1, ysize+1);

xmax = max(max(xnodes));
ymax = max(max(ynodes));

X = []; Y = []; D = [];
for i = 1:xsize
  for j = 1:ysize
    x = xnodes(i:i+1,j:j+1);
    y = ynodes(i:i+1,j:j+1);
    x = x(:);
    y = y(:);
    x(3:4) = x(4:-1:3);
    d = density(i,j)^3;
    X = [X,x];
    Y = [Y,y];
    D = [D,d];
  end
end

figure
fill(X,Y,D)
xlim([0,xmax])
ylim([0,ymax])
colorbar
%view(2)
axis('equal','tight');
colormap(flipud(bone))

%figure
%imagesc(density.^3);
%colorbar
%colormap('bone')

end
