function plot_mesh2(mesh, varargin)
%
%   plot_mesh2(mesh, option)
%
%   plot the mesh
%
%   Input:
%         mesh       structure array with the following fields
%
%         mesh.p     Real np x 2
%                    array containing the x- and y- coordinates
%                    of the nodes
%
%         mesh.t     Integer nt x 3    or   nt x 6
%                    t(i,1:4) contains the indices of the vertices of
%                    triangle i. 
%                    If t is a nt x 6 array, t(i,4:6) contains 
%                    the indices of the edge mid-points of triangle i.
%
%         mesh.e     Integer ne x 3
%                    e(i,1:2) contains the indices of the vertices of
%                    edge i.
%
%         option     option = [i,j]  controls the labels of the plot:
%                    i = 0   no triangle numbers
%                    i = 1   print triangle numbers
%                    j = 0   no node numbers
%                    j = 1   print node numbers
%                    If GridPlot is called without option, then no labels are plotted
%
%  AUTHOR:  Matthias Heinkenschloss
%           Department of Computational and Applied Mathematics
%           Rice University
%           November 23, 2005


% the number of nodes per triangle is either 3 or 6. 
% If there are 6 nodes per triangle, then the first 
% three are the triangle vertices and the last three 
% are the midpoints of the edges.

% get number of triangles, nt, and 
% number of nodes per traingle, nn.
[nt,nn] = size(mesh.t);

if( nargin == 1 )
    nodelabel   = 0;
    trianglabel = 0;
else
    nodelabel   = varargin{1}(2);
    trianglabel = varargin{1}(1);
end

% Plot grid.
% Only the vertices are connected. The vertex numbers are
% are stored in t(it,1), t(it,2), t(it,3)
xx = zeros(1,4);
yy = zeros(1,4);
for it = 1:nt
    xx    = mesh.p(mesh.t(it,1:3),1);
    yy    = mesh.p(mesh.t(it,1:3),2);
    xx(4) = mesh.p(mesh.t(it,1),1);
    yy(4) = mesh.p(mesh.t(it,1),2);
    plot(xx,yy,'r'); hold on
end
axis( [min(mesh.p(:,1))  max(mesh.p(:,1))  min(mesh.p(:,2))  max(mesh.p(:,2))] );


if( nodelabel == 1 )
%   Label nodes (in green) 
    for i = 1:nt
        plot(mesh.p(mesh.t(i,(1:nn)),1),...
             mesh.p(mesh.t(i,(1:nn)),2),'go')
        for j = 1:nn
            x = 1.02*mesh.p(mesh.t(i,j),1);
            y = 1.02*mesh.p(mesh.t(i,j),2);
            h = text(x,y,int2str(mesh.t(i,j)));
            set(h,'Color','k','FontSize',12);
        end
    end
end

if( trianglabel == 1 )
%   Label triangles  (in blue)
    for i = 1:nt
        x = sum(mesh.p(mesh.t(i,1:3),1)) / 3;
        y = sum(mesh.p(mesh.t(i,1:3),2)) / 3;
        h = text(x,y,int2str(i));
        set(h,'Color','b');
    end
end
hold off

drawnow
