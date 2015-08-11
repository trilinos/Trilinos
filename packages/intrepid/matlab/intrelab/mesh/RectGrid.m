function [mesh] = RectGrid(xmin, xmax, ymin, ymax, nx, ny, celltype)
%
%  [mesh] = RectGrid(xmin, xmax, ymin, ymax, nx, ny, celltype)
%
%RECTGRID   sets up the grid in a rectangular domain.
%  
%  The grid is constructed by subdividing the  x-interval into
%  nx subintervals and the  y-interval into ny subintervals.
%  This generates a grid with nx*ny rectangles.
%  If celltype = 'Quadrilateral', these are the mesh cells.
%  If celltype = 'Triangle' each rectangle is subdivided into
%  two triangles by cutting the rectangle from bottom left to
%  top right.
%  
%
%
%  Input
%         xmin, xmax  size of the rectangular domain
%         ymin, ymax
%         nx          number of subintervals on x-interval
%         ny          number of subintervals on y-interval
%         celltype    'Triangle' or 'Quadrilateral'
%
%  Output
%         mesh       structure array with the following fields
%
%         mesh.p     Real num_nodes x 2
%                    array containing the x- and y- coordinates
%                    of the nodes
%
%         mesh.t     Integer num_cells x 3 or 4   
%                    - if celltype = 'Triangle' then t(i,1:3) contains
%                    the indices of the vertices of triangle i
%                    - if celltype = 'Quadrilateral' then t(i,1:4) contains
%                    the indices of the vertices of quadrilateral i
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
%         mesh.sidesPerCell  Integer = 3 (triangle) or 4 (quadrilateral)
%
%         mesh.elem_ss{i}    = global element IDs on sideset (boundary) i
%
%         mesh.side_ss{i}    = local side (subcell) IDs on sideset (boundary) i
%
%         mesh.cellType      = string for cell type (= 'Triangle' or 'Quadrilateral')
%
%         mesh.sideType      = string for side subcell type (= 'Line')
%     
%
%  Vertical ordering:
%  The triangles are ordered column wise, for instance:
% 
%    03 -------- 06 -------- 09 -------- 12
%     |  4     /  |  8     /  | 12     /  |
%     |     /     |     /     |     /     |
%     |  /    3   |  /    7   |  /    11  |
%    02 -------- 05 -------- 08 -------- 11
%     |  2     /  |  6     /  | 10     /  |      
%     |     /     |     /     |     /     |
%     |  /    1   |  /    5   |  /     9  |
%    01 -------- 04 -------- 07 -------- 10
%
%  The vertices in a triangle are numbered
%  counterclockwise, for example
%          triangle 7: (05, 08, 09)
%          triangle 8: (05, 09, 06)
%
%  number of triangles: 2*nx*ny,
%  number of vertices:  (nx+1)*(ny+1), 
%
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
%  number of quadrilaterals: nx*ny,
%  number of vertices:  (nx+1)*(ny+1), 
%
%  (Usually, the local grid.node numbering should not be important.)
%
%  AUTHORS:
%            Denis Ridzal
%            Sandia National Laboratories

% Numbers of points
np = (nx+1)*(ny+1);
nxp1 = nx + 1;
nyp1 = ny + 1;

if strcmp(lower(celltype), 'triangle')
  nt = 2*nx*ny;
  % Create connectivity array
  mesh.t = zeros(nt,3);

  % Create triangles
  nt  = 0;
  iyvec = [1:ny]';
  for ix = 1:nx
    mesh.t(2*(iyvec+(ix-1)*ny)-1,1) = iyvec+(ix-1)*(ny+1);
    mesh.t(2*(iyvec+(ix-1)*ny)-1,2) = iyvec+ix*(ny+1);
    mesh.t(2*(iyvec+(ix-1)*ny)-1,3) = iyvec+ix*(ny+1)+1;
    mesh.t(2*(iyvec+(ix-1)*ny),1)   = iyvec+(ix-1)*(ny+1);
    mesh.t(2*(iyvec+(ix-1)*ny),2)   = iyvec+ix*(ny+1)+1;
    mesh.t(2*(iyvec+(ix-1)*ny),3)   = iyvec+(ix-1)*(ny+1)+1;
  end
elseif strcmp(lower(celltype), 'quadrilateral')
  nt = nx*ny;
  % Create connectivity array
  mesh.t = zeros(nt,4);

  % Create quadrilaterals
  nt  = 0;
  iyvec = [1:ny]';
  for ix = 1:nx
    mesh.t(iyvec+(ix-1)*ny,1) = iyvec+(ix-1)*(ny+1);
    mesh.t(iyvec+(ix-1)*ny,2) = iyvec+ix*(ny+1);
    mesh.t(iyvec+(ix-1)*ny,3) = iyvec+ix*(ny+1)+1;
    mesh.t(iyvec+(ix-1)*ny,4) = iyvec+(ix-1)*(ny+1)+1;
  end
end

%mesh.t

% Create point array
mesh.p = zeros(np,2);

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

if strcmp(lower(celltype), 'triangle')
  % sides per cell
  mesh.sidesPerCell = 3;
  % cell type
  mesh.cellType     = 'Triangle';
  % side type
  mesh.sideType     = 'Line';

  % sidesets
  xvec = [1:nx]';
  yvec = [1:ny]';
  mesh.elem_ss{1}   = 2*ny*(xvec-1) + 1;             % bottom
  mesh.elem_ss{2}   = 2*ny*(nx-1) + 2*(yvec-1) + 1;  % right
  mesh.elem_ss{3}   = 2*ny*xvec;                     % top
  mesh.elem_ss{4}   = 2*yvec;                        % left

  % local side ids for sidesets
  mesh.side_ss{1}   = 0*ones(nx,1);  % bottom 
  mesh.side_ss{2}   = 1*ones(ny,1);  % right 
  mesh.side_ss{3}   = 1*ones(nx,1);  % top 
  mesh.side_ss{4}   = 2*ones(ny,1);  % left
elseif strcmp(lower(celltype), 'quadrilateral')
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
end


