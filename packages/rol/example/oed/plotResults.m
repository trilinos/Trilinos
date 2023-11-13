
tol = 1e-6;

figure,
% D
design = load('D_optimal_design_0.txt');
ind = design(:,2) > tol;
plot(design(ind,1),design(ind,2),'o',...
     'MarkerFaceColor', [1 0.5 0], ...
     'MarkerEdgeColor', [1 0.5 0], ...
     'MarkerSize', 8), hold on

% I
design = load('I_optimal_design_0.txt');
ind = design(:,2) > tol;
plot(design(ind,1),design(ind,2),'o',...
     'MarkerFaceColor', [0 0 1], ...
     'MarkerEdgeColor', [0 0 1], ...
     'MarkerSize', 8)

% A
design = load('A_optimal_design_0.txt');
ind = design(:,2) > tol;
plot(design(ind,1),design(ind,2),'o',...
     'MarkerFaceColor', [0 0.5 0], ...
     'MarkerEdgeColor', [0 0.5 0], ...
     'MarkerSize', 8)

plot(design(:,1),0*design(:,2),'kx'), hold off
xlim([-1 1])
ylim([0 1])
xlabel('$x$', 'fontsize', 16, 'interpreter', 'latex')
ylabel('$p$', 'fontsize', 16, 'interpreter', 'latex')
legend('D','I','A','fontsize',16)
set(gca, 'fontsize', 14)
print('-depsc2','optimalDesign.eps')

x = linspace(-1,1,1000).';
sig = @(t,a) exp(a * abs(t)) / exp(a);
figure,
plot(x, sig(x,5), 'k', 'linewidth', 3);
xlabel('$x$', 'fontsize', 16, 'interpreter', 'latex')
ylabel('$\sigma(x)$', 'fontsize', 16, 'interpreter', 'latex')
xlim([-1 1])
ylim([0 1])
set(gca, 'fontsize', 14)
print('-depsc2','noise.eps')
