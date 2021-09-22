
design = load('D_optimal_design_0.txt');

ind = design(:,2) > 0.1;

figure,
plot(design(ind,1),design(ind,2),'o',...
     'MarkerFaceColor', [1 0.5 0], ...
     'MarkerEdgeColor', [1 0.5 0], ...
     'MarkerSize', 8), hold on
plot(design(~ind,1),0*design(~ind,2),'kx'), hold off
xlim([-1 1])
ylim([0 1])
xlabel('$x$', 'fontsize', 16, 'interpreter', 'latex')
ylabel('$p$', 'fontsize', 16, 'interpreter', 'latex')
set(gca, 'fontsize', 14)

x = linspace(-1,1,1000).';
sig = @(t,a) exp(a * abs(t)) / exp(a);
figure,
plot(x, sig(x,5), 'k', 'linewidth', 3);
xlabel('$x$', 'fontsize', 16, 'interpreter', 'latex')
ylabel('$\sigma(x)$', 'fontsize', 16, 'interpreter', 'latex')
xlim([-1 1])
ylim([0 1])
set(gca, 'fontsize', 14)
