% Driver for Poisson problem.

clear all;
clc;

set(0, 'defaultaxesfontsize',12,'defaultaxeslinewidth',0.7,...
    'defaultlinelinewidth',0.8,'defaultpatchlinewidth',0.7,...
    'defaulttextfontsize',12);

% mesh utilities directory
addpath ../../mesh/

fprintf('\n*** Poisson Equation Solver ***\n');

% number of intervals in x and y direction on a rectangular grid
nxint = 256;
nyint = 256;

% generate rhs_fn = @generateSine;
rhs_fn=@(usr_par)rhs(usr_par);

% set up PDE
fprintf('\nSetting up problem and assembling FEM operators ...\n');
tic
[usr_par] = problemsetup(nxint, nyint, rhs_fn);
toc

% solve PDE
fprintf('\nSolving PDE ...\n');
tic
[u] = pdesolve(usr_par);
toc

figure(1)
trisurf(usr_par.mesh.t, usr_par.mesh.p(:,1), usr_par.mesh.p(:,2), ...
        u, 'facecolor','interp')
view(10,40);
shading interp;
title('PDE Solution')
xlabel('x')
ylabel('y')

figure(2)
trisurf(usr_par.mesh.t, usr_par.mesh.p(:,1), usr_par.mesh.p(:,2), ...
        usr_par.k, 'facecolor','interp')
view(10,40);
shading flat;
title('Diffusivity')
xlabel('x')
ylabel('y')
