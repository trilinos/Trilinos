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
%         mesh.t     Integer nt x nn
%                    t(i,1:nn) contains the indices of the vertices of
%                    cell i. 
%
 %         option     option = [i,j]  controls the labels of the plot:
%                    i = 0   no cell numbers
%                    i = 1   print cell numbers
%                    j = 0   no node numbers
%                    j = 1   print node numbers
%                    If GridPlot is called without option, then no labels are plotted
%
%  AUTHOR:  Matthias Heinkenschloss and Denis Ridzal
%           Department of Computational and Applied Mathematics
%           Rice University
%           November 23, 2005


% get number of cells, nt, and 
% number of nodes per cell, nn.
[nt,nn] = size(mesh.t)

if( nargin == 1 )
    nodelabel   = 0;
    trianglabel = 0;
else
    nodelabel   = varargin{1}(2);
    trianglabel = varargin{1}(1);
end

% Plot grid.
xx = zeros(1,4);
yy = zeros(1,4);
for it = 1:nt
    xx    = mesh.p(mesh.t(it,1:nn),1);
    yy    = mesh.p(mesh.t(it,1:nn),2);
    xx(nn+1) = mesh.p(mesh.t(it,1),1);
    yy(nn+1) = mesh.p(mesh.t(it,1),2);
    plot(xx,yy,'r'); hold on
end
axis( [min(mesh.p(:,1))  max(mesh.p(:,1))  min(mesh.p(:,2))  max(mesh.p(:,2))] );


if( nodelabel == 1 )
%   Label nodes (in green) 
    for i = 1:nt
        plot(mesh.p(mesh.t(i,(1:nn)),1),...
             mesh.p(mesh.t(i,(1:nn)),2),'go')
        for j = 1:nn
            x = mesh.p(mesh.t(i,j),1);
            y = mesh.p(mesh.t(i,j),2);
            h = text(x,y,int2str(mesh.t(i,j)));
            set(h,'Color','k','FontSize',11,'FontWeight','bold');
        end
    end
end

if( trianglabel == 1 )
%   Label cells  (in blue)
    for i = 1:nt
        x = sum(mesh.p(mesh.t(i,1:nn),1)) / nn;
        y = sum(mesh.p(mesh.t(i,1:nn),2)) / nn;
        h = text(x,y,int2str(i));
        set(h,'Color','b','FontSize',11);
    end
end
hold off

drawnow
