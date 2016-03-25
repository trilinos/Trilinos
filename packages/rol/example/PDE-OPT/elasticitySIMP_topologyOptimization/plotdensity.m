data_obj = importdata('density.txt', ' ', 2);  %% we need to skip the first two lines
density = data_obj.data;
nodes = load('nodes.txt');

xsize = round(sqrt(size(density, 1)));
ysize = xsize;

density = reshape(density, xsize, ysize);
Xnodes  = reshape(nodes(:,1), xsize+1, ysize+1);
Ynodes  = reshape(nodes(:,2), xsize+1, ysize+1);

Xnodes  = 0.5*(Xnodes(1:end-1,1:end-1)+Xnodes(2:end,1:end-1));
Ynodes  = 0.5*(Ynodes(1:end-1,1:end-1)+Ynodes(1:end-1,2:end));

figure
fill(Xnodes,Ynodes,density.^3)
colorbar
%view(2)
axis('equal','tight');
colormap('bone')

figure
imagesc(density.^3);
colorbar
colormap('bone')
