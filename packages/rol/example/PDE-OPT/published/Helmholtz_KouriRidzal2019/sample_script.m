
Adata = importdata('Amatrix.txt',' ',2); Adata = Adata.data;
Bdata = importdata('Bmatrix.txt',' ',2); Bdata = Bdata.data;
Ldata = importdata('Lmatrix.txt',' ',2); Ldata = Ldata.data;
Cdata = importdata('Cmatrix.txt',' ',2); Cdata = Cdata.data;
Rdata = importdata('Rmatrix.txt',' ',2); Rdata = Rdata.data;
Mdata = importdata('Mmatrix.txt',' ',2); Mdata = Mdata.data;

wrvec = importdata('WRvector.txt',' ',2); wrvec = wrvec.data;
wivec = importdata('WIvector.txt',' ',2); wivec = wivec.data;

Amat = sparse(Adata(:,1),Adata(:,2),Adata(:,3));
Bmat = sparse(Bdata(:,1),Bdata(:,2),Bdata(:,3));
Lmat = sparse(Ldata(:,1),Ldata(:,2),Ldata(:,3));
Cmat = sparse(Cdata(:,1),Cdata(:,2),Cdata(:,3));
Rmat = sparse(Rdata(:,1),Rdata(:,2),Rdata(:,3));
Mmat = sparse(Mdata(:,1),Mdata(:,2),Mdata(:,3));

Kmat = Amat - sqrt(-1)*Bmat;

beta = 1e-4;

Zmat = 0*Kmat;
kkt  = [ Cmat      Zmat Kmat';
         Zmat beta*Rmat Lmat';
         Kmat      Lmat Zmat ];

w0   = wrvec + sqrt(-1)*wivec;
zvec = 0*w0;
rhs  = [ w0;
         zvec;
         zvec ];

sol  = kkt \ rhs;

adj   = load('cell_to_node_quad.txt') + 1;       %% load node adjacency, increment by 1 for 1-based indexing
nodes = load('nodes.txt');                       %% load node coordinates
map   = importdata('map.txt');
map   = map.data(1:2:end)+1;
[tmp, perm] = sort(map);
n     = length(w0);
u     = sol(1:n);
u     = u(perm);
z     = sol(n+1:2*n);
z     = z(perm);

figure, trisurf(adj, nodes(:,1), nodes(:,2), real(u));
shading interp;
view(0,90)
axis('equal','tight');
xlabel('x');
ylabel('y');
title('State: Real Part');
axis square

figure, trisurf(adj, nodes(:,1), nodes(:,2), imag(u));
shading interp;
view(0,90)
axis('equal','tight');
xlabel('x');
ylabel('y');
title('State: Imaginary Part');
axis square

figure, trisurf(adj, nodes(:,1), nodes(:,2), real(z));
shading interp;
view(0,90)
axis('equal','tight');
xlabel('x');
ylabel('y');
title('Control: Real Part');
axis square

figure, trisurf(adj, nodes(:,1), nodes(:,2), imag(z));
shading interp;
view(0,90)
axis('equal','tight');
xlabel('x');
ylabel('y');
title('Control: Imaginary Part');
axis square
