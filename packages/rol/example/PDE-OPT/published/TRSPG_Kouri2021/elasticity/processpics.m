addpath ~/software/export_fig
plotdensity
set(gcf,'color','white')
set(gca,'fontsize',18)
export_fig 'density' -png -painters -m3 -a3
