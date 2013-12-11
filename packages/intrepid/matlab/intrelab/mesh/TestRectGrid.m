addpath '../metis'

clear all;

xmin = 0; xmax = 1;
ymin = 0; ymax = 1;

nsdx = 2;
nsdy = 3;

nx = nsdx*2;
ny = nsdy*2;


[mesh] = RectGrid(xmin, xmax, ymin, ymax, nx, ny);

figure(1); clf
plot_mesh2(mesh, [1,1])

figure(2); clf
[dd_mesh] = RectGridDD(nsdx, nsdy, mesh.t, mesh.p, mesh.e);


rmpath '../metis'


for i = 1:nsdx*nsdy
    ifedge = dd_mesh{i}.ifedge;
    fprintf(1,'subdomain  %d \n ', i)
    for j = 1:length(ifedge);
        % Compute the outward unit normal.
        evec = mesh.p(ifedge(j,2),:)-mesh.p(ifedge(j,1),:);
        % The normal points to the right of the initial direction.
        nvec = (1/norm(evec))*[evec(2);-evec(1)];
        
        fprintf(1,'interface edge with nodes  %4d %4d \n ', ifedge(j,:))
        fprintf(1,'outward unit normal    %4.2f %4.2f \n ', nvec)
       
    end
end 
