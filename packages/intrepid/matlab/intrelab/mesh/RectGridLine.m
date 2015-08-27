function [mesh] = RectGridLine(xmin, xmax, nx)
%
%  [mesh] = RectGridLine(xmin, xmax, nx)
%
%RECTGRIDLINE   sets up the grid on a line.
%  
%  The grid is constructed by subdividing the  x-interval into
%  nx subintervals.
%
%
%  Input
%         xmin, xmax  size of the domain
%         nx          number of subintervals on x-interval
%
%  Output
%         mesh       structure array with the following fields
%
%         mesh.p     Real num_nodes x 1
%                    array containing the x-coordinates
%                    of the nodes
%
%         mesh.t     Integer num_cells x 2   
%                    - t(i,1:2) contains the indices of the
%                    vertices of line segment i
%
%         mesh.e     empty
%
%         mesh.sidesPerCell  Integer = 2
%
%         mesh.elem_ss{i}    = global element IDs on sideset (boundary) i
%
%         mesh.side_ss{i}    = local side (subcell) IDs on sideset (boundary) i
%
%         mesh.cellType      = string for cell type (= 'Line')
%
%         mesh.sideType      = string for side subcell type (= 'Vertex')
%     
%
%  AUTHORS:
%            Denis Ridzal
%            Sandia National Laboratories

np = nx+1;

nt = nx;
% Create connectivity array
mesh.t = zeros(nt,2);
mesh.t(:,1) = [1:nx]';
mesh.t(:,2) = [2:nx+1]';

% Create point array
mesh.p = zeros(np,1);

% Create vertex coodinates
hx   = (xmax-xmin)/nx;
mesh.p(:,1) = [xmin:hx:xmax]';

mesh.e = [];

% sides per cell
mesh.sidesPerCell = 2;
% cell type
mesh.cellType     = 'Line';
% side type
mesh.sideType     = 'Vertex';

% sidesets
mesh.elem_ss{1}   = 1;    % left
mesh.elem_ss{2}   = nx;   % right

% local side ids for sidesets
mesh.side_ss{1}   = 0;  % left 
mesh.side_ss{2}   = 1;  % right 

end
