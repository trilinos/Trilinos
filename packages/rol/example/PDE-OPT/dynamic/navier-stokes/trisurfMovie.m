function trisurfMovie(adj, nodes, data, name, T, xlimit, ylimit)
[nx,nt] = size(data);

minD = min(min(data));
maxD = max(max(data));

writerObj = VideoWriter([name,'.avi']);
writerObj.FrameRate = nt/T;
open(writerObj);

%% Set up figure handle
axsize      = 400;
h.fig       = figure('Position', [100 100 4*axsize 4*axsize], 'Color', 'White');
h.patch     = trisurf(adj, nodes(:,1), nodes(:,2), data(:,1));
h.ax        = h.fig.CurrentAxes;
h.EdgeColor = 'none';
%% Update patch and store movie
for i=1:nt
  h.patch = handle(trisurf(adj, nodes(:,1), nodes(:,2), data(:,i)));
  axis(h.ax, 'equal', 'tight');
  shading(h.ax, 'interp');
  colormap(h.ax, 'bone');
  h.ax.XTickMode      = 'Manual';
  h.ax.YTickMode      = 'Manual';
  h.ax.XTickLabelMode = 'Manual';
  h.ax.YTickLabelMode = 'Manual';
  h.ax.XLim           = xlimit;
  h.ax.YLim           = ylimit;
  h.ax.CLim           = [minD, maxD];
  h.ax.View           = [0, 90];
  h.ax.Visible        = 'off';

  drawnow();

  writeVideo(writerObj, getframe(gcf));

  if (mod(i,10)==0)
    fprintf('Frame %d  %1.8e\n',i,norm(data(:,i)));
  end
end
close(writerObj);
