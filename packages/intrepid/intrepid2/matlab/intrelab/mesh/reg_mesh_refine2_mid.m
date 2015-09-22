function [mesh1, I1] = reg_mesh_refine2_mid(mesh)
%FIX THE COMMENTS LATER
%  [mesh1,I1]=reg_mesh_refine2(mesh)
%
%       Compute a regular refinement of a mesh given by p, e, t.
%       Each triangle of the original mesh is refined into four triangles
%       by dividing each edge into two.
%
%  Input
%     mesh  structure with the fields mesh.p, mesh.e, mesh.t
%           If size(mesh.t, 2) == 3
%               mesh.p(i, 1:2) x-, y-coordinates of the i-th vertex
%
%               mesh.e(i, 1:2)  indices of vertices of boundary edge i.
%               mesh.e(i, 3)    boundary marker of edge i
%        
%               mesh.t(i, 1:3)  indices of vertices in triangle i.
%           If size(mesh.t, 2) == 6
%               mesh.p(i, 1:2) x-, y-coordinates of the i-th vertex
%
%               mesh.e(i, 1:2)  indices of vertices of boundary edge i.
%               mesh.e(i, 3)    indices of the mid-point of boundary edge i.
%               mesh.e(i, 4)    boundary marker of edge i
%        
%               mesh.t(i, 1:3)  indices of vertices in triangle i.
%               mesh.t(i, 4:6)  indices of mid-points in triangle i.
%
%  Output
%     mesh1 structure with the fields mesh1.p, mesh1.e, mesh1.t
%           If size(mesh.t, 2) == 3
%               mesh1.p(i, 1:2) x-, y-coordinates of the i-th vertex in 
%               the refined mesh.  
%               If np = size(mesh.p, 1), then mesh1.p(1:np, :) = mesh.p.
%
%               mesh1.e(i, 1:2)  indices of vertices of boundary edge i of 
%               the refined mesh
%               mesh1.e(i, 3)    boundary marker of edge i of the refined mesh.
%               If ne = size(mesh.e, 1), then mesh1.e(1:ne, :) = mesh.e.
%        
%               mesh1.t(i, 1:3)  indices of vertices in triangle i of the 
%               refined mesh.
%               If nt = size(mesh.t, 1), then mesh1.t(1:nt, :) = mesh.t.
%           If size(mesh.t, 2) == 6
%               mesh1.p(i, 1:2) x-, y-coordinates of the i-th vertex in 
%               the refined mesh.  
%               If np = size(mesh.p, 1), then mesh1.p(1:np, :) = mesh.p.
%
%               mesh1.e(i, 1:2)  indices of vertices of boundary edge i of 
%               the refined mesh
%               mesh1.e(i, 3)    indices of mid-point of boundary edge i of 
%               the refined mesh
%               mesh1.e(i, 4)    boundary marker of edge i of the refined
%               mesh.
%        
%               mesh1.t(i, 1:3)  indices of vertices in triangle i of the 
%               refined mesh.
%               mesh1.t(i, 4:6)  indices of mid-points in triangle i of the 
%               refined mesh.
%
%     I1    Interpolation matrix of size np1 x np, where np1 is the number
%           of vertices in the refined mesh and np is the number of
%           vertices in the original mesh. 
%           If u is the vector of function values of a piecewise linear function 
%           on the original mesh (u(i) is the function value at node
%           mesh.p(i,:)), then u1 = I1*u is the vector of function values of the 
%           piecewise linear function on the refined mesh (u1(i) is the function 
%           value at node mesh1.p(i,:)).
%    
%  reg_mesh_refine2_mid is based on a modification of refinemesh.m in the Matlab PDEToolbox. 
%
%  AUTHOR:  Patricia A. Howard
%           Department of Computational and Applied Mathematics
%           Rice University
%           March 3, 2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       

np = size(mesh.p, 1);
nt = 0;
ne = 0;

% Cannot use matrix indices that exceeds the size of a signed int
[comp, maxsize] = computer;
indexproblem = np^2 > maxsize;

% Make a connectivity matrix, with edges to be refined.
% -1 means no point is yet allocated
it = [];     
ip1 = [];
ip2 = [];
ip3 = []; 
ip4 = [];
ip5 = [];
ip6 = [];
if size(mesh.t, 2) == 3
    ne = size(mesh.e, 1);
    nt = size(mesh.t, 1);
    it =(1:nt);    
    ip1 = mesh.t(it, 1);
    ip2 = mesh.t(it, 2);
    ip3 = mesh.t(it, 3);
elseif size(mesh.t, 2) == 6
    ip1 = mesh.t(:, 1);
    ip2 = mesh.t(:, 2);
    ip3 = mesh.t(:, 3);
    ip4 = mesh.t(:, 4);
    ip5 = mesh.t(:, 5);
    ip6 = mesh.t(:, 6);
    % Construct new triangles with existing mid-point information
    mesh1.t = [ip1 ip6 ip5;
               ip2 ip4 ip6;
               ip4 ip5 ip6;
               ip3 ip5 ip4];
    nt = size(mesh1.t, 1);
    it =(1:nt);    
    ip1 = mesh1.t(it, 1);
    ip2 = mesh1.t(it, 2);
    ip3 = mesh1.t(it, 3);
    ne = 2 * size(mesh.e, 1);
    % Construct new edges with existing mid-point information
    mesh1.e = [mesh.e(:, 1) mesh.e(:, 3) (-1 * ones(length(mesh.e(:, 1)), 1)) mesh.e(:, 4);
               mesh.e(:, 3) mesh.e(:, 2) (-1 * ones(length(mesh.e(:, 2)), 1)) mesh.e(:, 4)];
end

A = sparse(ip1, ip2, -1, np, np);
A = A + sparse(ip2, ip3, -1, np, np);
A = A + sparse(ip3, ip1, -1, np, np);
A = -((A + A.') < 0);

% Generate points on (interior and boundary) edges
% Determine (interior and boundary) edges
[i1, i2] = find(A == -1 & A.' == -1);
i = find(i2 > i1);
i1 = i1(i);
i2 = i2(i);
% Edges (interior and boundary) have vertices i1 i2.
% Compute midpoints
mesh1.p = [mesh.p; ((mesh.p(i1,1:2) + mesh.p(i2,1:2))/2)];

% Fill in the new points. 
% The midpoint of edge with vertices i1 and i2 gets index ip
ip = ((np + 1):(np + length(i)))';
if ~indexproblem
  A(i1 + np * (i2 - 1)) = ip;
  A(i2 + np * (i1 - 1)) = ip;
else
  A=l_assign(A, [i1 i2], [i2 i1], [ip ip]);
end

% Construct interpolation matrix
I1 = [];
if size(mesh.t, 2) == 3
    irowI1 = [(1:np)'];
    icolI1 = [(1:np)'];
    sI1    = [ones(np, 1)];

    irowI1 = [irowI1; ip; ip];
    icolI1 = [icolI1; i1; i2];
    sI1    = [sI1; 0.5 * ones(2 * length(i), 1)];

    I1 = sparse(irowI1, icolI1, sI1);
elseif size(mesh.t, 2) == 6
    % Constructs an interpolation matrix for linear elements
    % WITHOUT mid-point (as is required in source-inv)
    nv = max(max(mesh.t(:, 1:3)));
    row_col_val = [(1:nv)' (1:nv)' ones(nv, 1)];
    
    row_col_val = [row_col_val;
                   ip6 mesh.t(:, 1) 0.5 * ones(length(ip6), 1);
                   ip6 mesh.t(:, 2) 0.5 * ones(length(ip6), 1);
                   ip4 mesh.t(:, 2) 0.5 * ones(length(ip4), 1);
                   ip4 mesh.t(:, 3) 0.5 * ones(length(ip4), 1);
                   ip5 mesh.t(:, 3) 0.5 * ones(length(ip5), 1);
                   ip6 mesh.t(:, 1) 0.5 * ones(length(ip5), 1)];
                   
    row_col_val = unique(row_col_val, 'rows');
    
    I1 = sparse(row_col_val(:, 1), row_col_val(:, 2), row_col_val(:, 3));
end

% Form the new triangles (if size(mesh.t, 2) = 3) or 
% new mid-points (if size(mesh.t, 2) = 6)
if ~indexproblem
    mp1 = full(A(ip2 + np * (ip3 - 1)));  % mp1 is index of midpoint of edge opposite of vectex 1
    mp2 = full(A(ip3 + np * (ip1 - 1)));  % mp2 is index of midpoint of edge opposite of vectex 2
    mp3 = full(A(ip1 + np * (ip2 - 1)));  % mp3 is index of midpoint of edge opposite of vectex 3
else
    mp1 = l_extract(A, ip2, ip3);
    mp2 = l_extract(A, ip3, ip1);
    mp3 = l_extract(A, ip1, ip2);
end

if size(mesh.t, 2) == 3
    mesh1.t = zeros(4 * nt, 3);
    i  = (1:nt);
    nt1 = 0;
    mesh1.t((nt1 + 1):(nt1 + length(i)), :)=[mesh.t(it(i), 1) mp3(i) mp2(i)];
    nt1 = nt1 + length(i);
    mesh1.t((nt1 + 1):(nt1 + length(i)), :) = [mesh.t(it(i), 2) mp1(i) mp3(i)];
    nt1 = nt1 + length(i);
    mesh1.t((nt1 + 1):(nt1 + length(i)), :) = [mesh.t(it(i), 3) mp2(i) mp1(i)];
    nt1 = nt1 + length(i);
    mesh1.t((nt1 + 1):(nt1 + length(i)), :) = [mp1(i) mp2(i) mp3(i)];
    nt1 = nt1 + length(i);
elseif size(mesh.t, 2) == 6
    i  = (1:nt);
    mesh1.t(it, 4) = mp1(i);
    mesh1.t(it, 5) = mp2(i);
    mesh1.t(it, 6) = mp3(i);
end


% Form new edges (if size(mesh.t, 2) = 3) or add mid-points (if
% size(mesh.t, 2) = 6)
ie = (1:ne);
if size(mesh.t, 2) == 3
    if ~indexproblem
        mp1 = full(A(mesh.e(ie, 1) + np * (mesh.e(ie, 2) - 1)));  %mp1 is index of midpoint of edge ie
    else
        mp1 = l_extract(A, mesh.e(ie, 1), mesh.e(ie, 2));
    end

    % Create new edges
    mesh1.e = [[mesh.e(ie,1) mp1 mesh.e(ie,3)]; ...
               [mp1 mesh.e(ie,2) mesh.e(ie,3)]];
elseif size(mesh.t, 2) == 6
    if ~indexproblem
        mp1 = full(A(mesh1.e(ie, 1) + np * (mesh1.e(ie, 2) - 1)));  %mp1 is index of midpoint of edge ie
    else
        mp1 = l_extract(A, mesh1.e(ie, 1), mesh1.e(ie, 2));
    end
    
    mesh1.e(:, 3) = mp1;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function k=l_extract(A,i,j)

if numel(i)~=numel(j)
  error('PDE:refinemesh:ijNumel', 'i and j must have the same number of elements.')
end

k=zeros(size(i));

for l=1:numel(i)
  k(l)=A(i(l),j(l));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A=l_assign(A,i,j,k)

if numel(i)~=numel(j) || numel(i)~=numel(k) 
  error('PDE:refinemesh:ijkNumel', 'i, j, and k must have the same number of elements.')
end

for l=1:numel(i)
  A(i(l),j(l))=k(l);
end
