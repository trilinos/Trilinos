
function [mesh] = RectGridQuad(xmin, xmax, ymin, ymax, nx, ny)
%
%  [mesh] = RectGridQuad(xmin, xmax, ymin, ymax, nx, ny)
%
%RECTGRID   sets up the grid for piecewise linear elements
%  in a rectangular domain on quadrilateral cells.
%  
%  The grid is constructed by subdividing the  x-interval into
%  nx subintervals and the  y-interval into ny subintervals.
%  This generates a grid with nx*ny rectangular quadrilaterals.
%  
%
%
%  Input
%         xmin, xmax  size of the rectangle
%         ymin, ymax
%         nx          number of subintervals on x-interval
%         ny          number of subintervals on y-interval
%
%  Output
%         mesh       structure array with the following fields
%
%         mesh.p     Real num_nodes x 2
%                    array containing the x- and y- coordinates
%                    of the nodes
%
%         mesh.t     Integer num_cells x 4   
%                    t(i,1:4) contains the indices of the vertices of
%                    quadrilateral i. 
%
%         mesh.e     Integer num_edges x 3
%                    e(i,1:2) contains the indices of the vertices of
%                    edge i.
%                    edge(i,3) contains the boundary marker of edge i.
%                    Currently set to one.
%                    e(i,3) = 1  Dirichlet bdry conds are imposed on edge i
%                    e(i,3) = 2  Neumann bdry conds are imposed on edge i
%                    e(i,3) = 3  Robin bdry conds are imposed on edge i
%
%         mesh.sidesPerCell  Integer = 4 (Quadrilateral)
%
%         mesh.elem_ss{i}    = global element IDs on sideset (boundary) i
%
%         mesh.side_ss{i}    = local side (subcell) IDs on sideset (boundary) i
%
%         mesh.cellType      = string for cell type (= 'Quadrilateral')
%
%         mesh.sideType      = string for side subcell type (= 'Line')
%     
%
%  Vertical ordering:
%  The quadrilaterals are ordered column wise, for instance:
% 
%    03 -------- 06 -------- 09 -------- 12
%     |           |           |           |
%     |     2     |     4     |     6     |
%     |           |           |           |
%    02 -------- 05 -------- 08 -------- 11
%     |           |           |           |      
%     |     1     |     3     |     5     |
%     |           |           |           |
%    01 -------- 04 -------- 07 -------- 10
%
%  The vertices in a quadrilateral are numbered
%  counterclockwise, for example
%          quadrilateral 3: (04, 07, 08, 05)
%          quadrilateral 6: (08, 11, 12, 09)
%
%  (Usually, the local grid.node numbering should not be important.)
%
%  number of quadrilaterals: nx*ny,
%  number of vertices:  (nx+1)*(ny+1), 
%
%  AUTHOR: 
%            Denis Ridzal
%            Sandia National Laboratories


nt = nx*ny;
np = (nx+1)*(ny+1);

% Create arrays
mesh.t = zeros(nt,4);
mesh.p = zeros(np,2);

nxp1 = nx + 1;
nyp1 = ny + 1;

% Create quadrilaterals
nt  = 0;
for ix = 1:nx
   for iy = 1:ny
      
      iv  = (ix-1)*nyp1 + iy;
      iv1 = iv + nyp1;
      
      nt = nt + 1;
      mesh.t(nt,1) = iv;
      mesh.t(nt,2) = iv1;
      mesh.t(nt,3) = iv1+1;
      mesh.t(nt,4) = iv+1;

  end
end

% Create vertex coodinates

hx   = (xmax-xmin)/nx;
hy   = (ymax-ymin)/ny;
x    = xmin;

for ix = 1:nx

  % set coordinates for vertices with fixed 
  % x-coordinate at x
  i1 = (ix-1)*(ny+1)+1;
  i2 = ix*(ny+1);
  mesh.p(i1:i2,1) = x*ones(nyp1,1);
  mesh.p(i1:i2,2) = (ymin:hy:ymax)';
   
  x = x + hx;
end

% set coordinates for vertices with fixed 
% x-coordinate at xmax
i1 = nx*(ny+1)+1;
i2 = (nx+1)*(ny+1);
mesh.p(i1:i2,1) = xmax*ones(nyp1,1);
mesh.p(i1:i2,2) = (ymin:hy:ymax)';
  
% Set grid.edge (edges are numbered counter clock wise starting
% at lower left end).

mesh.e = ones(2*(nx+ny),3);

% edges on left on left boundary
mesh.e(1:ny,1) = (1:ny)';
mesh.e(1:ny,2) = (2:ny+1)';

% edges on top boundary
mesh.e(ny+1:nx+ny,1) = (ny+1:ny+1:np-1)';
mesh.e(ny+1:nx+ny,2) = (2*(ny+1):ny+1:np)';

% edges on right boundary
mesh.e(nx+ny+1:nx+2*ny,1) = (np-ny:np-1)';
mesh.e(nx+ny+1:nx+2*ny,2) = (np-ny+1:np)';

% edges on lower boundary
mesh.e(nx+2*ny+1:2*(nx+ny),1) = (1:ny+1:np-2*ny-1)';
mesh.e(nx+2*ny+1:2*(nx+ny),2) = (ny+2:ny+1:np-ny)';

% sides per cell
mesh.sidesPerCell = 4;

% cell type
mesh.cellType     = 'Quadrilateral';

% side type
mesh.sideType     = 'Line';

% sidesets
xvec = [1:nx]';
yvec = [1:ny]';
mesh.elem_ss{1}   = ny*(xvec-1) + 1;             % bottom
mesh.elem_ss{2}   = ny*(nx-1) + (yvec-1) + 1;    % right
mesh.elem_ss{3}   = ny*xvec;                     % top
mesh.elem_ss{4}   = yvec;                        % left

% local side ids for sidesets
mesh.side_ss{1}   = 0*ones(nx,1);  % bottom 
mesh.side_ss{2}   = 1*ones(ny,1);  % right 
mesh.side_ss{3}   = 2*ones(nx,1);  % top 
mesh.side_ss{4}   = 3*ones(ny,1);  % left

