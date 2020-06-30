%% function [niter] = run(beta,wave,fullObs,fullCtrl,solver)
%% Inputs: beta     - control penalty parameter
%%         wave     - wavenumber
%%         fullObs  - enable full observations, true/false
%%         fullCtrl - enable full controls, true/false
%%         solver   - solver:
%%                    'DIRECT'   for direct KKT solve
%%                    'PRDW'     for (perturbed) Rees-Dollar-Wathen preconditioner with SYMMLQ
%%                    'IMSHIFT'  for imaginary shift preconditioner with SYMMLQ
%% Output: niter - number of solver iterations
function [niter,sol] = run(beta,wave,fullObs,fullCtrl,solver)

% Global data structure.
global GLB_PROB;

% Solver tolerance.
tol = min(1e-4, 0.1*beta*wave^2/(1+wave^2));

% Solver maximum number of iterations and Krylov subspace size.
maxit = 600;

% Plot solution?
iplot = false;

% Print output?
iprint = true;

force_full_obs  = fullObs;
force_full_ctrl = fullCtrl;
compute_reference = false;
GLB_PROB.reuse_factors = true;

%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Load matrices from files %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iprint, fprintf('Reading problem matrices from file ...\n'), end;

Adata = importdata('Amatrix.txt',' ',2); Adata = Adata.data;
Amat = sparse(Adata(:,1),Adata(:,2),Adata(:,3));
[s1,s2] = size(Amat);
clear Adata;
if iprint, fprintf('A matrix is %d by %d.\n', s1, s2), end;
m = s1;

Bdata = importdata('Bmatrix.txt',' ',2); Bdata = Bdata.data;
Bmat = sparse(Bdata(:,1),Bdata(:,2),Bdata(:,3));
[s1,s2] = size(Bmat);
clear Bdata;
if iprint, fprintf('B matrix is %d by %d.\n', s1, s2), end;

Ldata = importdata('Lmatrix.txt',' ',2); Ldata = Ldata.data;
Lmat = sparse(Ldata(:,1),Ldata(:,2),Ldata(:,3));
nzcols = any(Lmat,1);
ctrlIdx = find(nzcols)';
Lmat(:, ~any(Lmat,1)) = []; % remove zero columns
[s1,s2] = size(Lmat);
clear Ldata;
if iprint, fprintf('L matrix is %d by %d.\n', s1, s2), end;
n = s2;

Cdata = importdata('Cmatrix.txt',' ',2); Cdata = Cdata.data;
Cmat = sparse(Cdata(:,1),Cdata(:,2),Cdata(:,3));
[s1,s2] = size(Cmat);
clear Cdata;
if iprint, fprintf('C matrix is %d by %d.\n', s1, s2), end;

Rdata = importdata('Rmatrix.txt',' ',2); Rdata = Rdata.data;
Rmat = sparse(Rdata(:,1),Rdata(:,2),Rdata(:,3));
Jmat = Rmat;
Rmat(:, ~any(Rmat,1)) = []; % remove zero columns
Rmat(~any(Rmat,2), :) = []; % remove zero rows
[s1,s2] = size(Rmat);
clear Rdata;
if iprint, fprintf('R matrix is %d by %d.\n', s1, s2), end;

Mdata = importdata('Mmatrix.txt',' ',2); Mdata = Mdata.data;
Mmat = sparse(Mdata(:,1),Mdata(:,2),Mdata(:,3));
[s1,s2] = size(Mmat);
clear Mdata;
if iprint, fprintf('M matrix is %d by %d.\n', s1, s2), end;

if (force_full_obs)
  Cmat = Mmat;
end

if (force_full_ctrl)
  Rmat = Mmat;
  Lmat = Mmat;
  n = m;
  ctrlIdx = [1:m]';
end

wrvec = importdata('WRvector.txt',' ',2); wrvec = wrvec.data;
wivec = importdata('WIvector.txt',' ',2); wivec = wivec.data;

%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Build optimality system %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Kmat = Amat - sqrt(-1)*Bmat;

kkt  = [ Cmat         sparse(m,n)  Kmat';
         sparse(n,m)  beta*Rmat    Lmat';
         Kmat         Lmat         sparse(m,m) ];

w    = wrvec + sqrt(-1)*wivec;
rhs  = [ w;
         zeros(n,1);
         zeros(m,1)];

GLB_PROB.ct = 0;
GLB_PROB.beta = beta;
GLB_PROB.nu = m;
GLB_PROB.nz = n;
GLB_PROB.K = Kmat;
GLB_PROB.C = Cmat;
GLB_PROB.R = Rmat;
GLB_PROB.J = Jmat;
GLB_PROB.M = Mmat;
GLB_PROB.L = Lmat;

%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Solve optimality system %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GLB_PROB.factors_computed = false;

if compute_reference
  fprintf("\nComputing reference solution using backslash...\n");
  tic
  ref_sol = kkt \ rhs;
  toc
  ref_z   = ref_sol(m+1:m+n, 1);
end

if (strcmp(lower(solver),'direct'))

  fprintf("\nComputing direct solution using backslash...\n");
  tic
  sol  = kkt \ rhs;
  toc

elseif (strcmp(lower(solver),'imshift') || strcmp(lower(solver),'prdw'))

  fprintf("\nComputing solution using iterative solver ...\n");

  % Set initial guess.
  t0  = zeros(2*m+n,1);

  % Choose method.
  if strcmp(lower(solver),'imshift')
    prec = @imshift;
  else
    prec = @prdw;
  end

  % MATLAB's SYMMLQ
  myprec = prec;
  tic
  [sol, flag, relres, niter, resvec] ...
    = symmlq(kkt, rhs, tol, maxit, myprec, [], t0, GLB_PROB);
  fprintf("\n");
  toc

  if iplot
    semilogy(resvec/resvec(1), 'LineWidth', 2)
  end
  % Compute final relative residual and control error, if enabled.
  fprintf('Relative residual = %12.10e\n', norm(kkt*sol - rhs)/norm(rhs));
  if compute_reference
    fprintf('Relative control error = %12.10e\n', norm(sol(m+1:m+n, 1) - ref_z)/norm(ref_z));
  end

end

%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Plot optimal solution %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (iplot)

  adj   = load('cell_to_node_quad.txt') + 1;       %% load node adjacency, increment by 1 for 1-based indexing
  nodes = load('nodes.txt');                       %% load node coordinates
  map   = importdata('map.txt');
  map   = map.data(1:2:end)+1;
  [tmp, perm] = sort(map);
  sol_plot = NaN(3*m, 1)+sqrt(-1)*NaN(3*m, 1);
  sol_plot(1:m, 1) = sol(1:m, 1);
  sol_plot(m+ctrlIdx, 1) = sol(m+1:m+n, 1);
  sol_plot(2*m+1:end, 1) = sol(m+n+1:end, 1);
  u     = sol_plot(1:m);
  u     = u(perm);
  z     = sol_plot(m+1:2*m);
  z     = z(perm);

  figure, trisurf(adj, nodes(:,1), nodes(:,2), real(u));
  shading interp;
  view(0,90)
  axis('equal','tight');
  xlabel('x');
  ylabel('y');
  title('State: Real Part');
  axis square
  x = 0; y = 0; r = 2;
  hold on
  th = [0:pi/50:2*pi]';
  xunit = r * cos(th) + x;
  yunit = r * sin(th) + y;
  h = plot3(xunit, yunit, 2*ones(size(xunit)), 'k-', 'LineWidth', 2);
  hold off

  figure, trisurf(adj, nodes(:,1), nodes(:,2), imag(u));
  shading interp;
  view(0,90)
  axis('equal','tight');
  xlabel('x');
  ylabel('y');
  title('State: Imaginary Part');
  axis square
  hold on
  h = plot3(xunit, yunit, 2*ones(size(xunit)), 'k-', 'LineWidth', 2);
  hold off

  figure, trisurf(adj, nodes(:,1), nodes(:,2), real(z));
  shading interp;
  view(0,90)
  axis('equal','tight');
  xlabel('x');
  ylabel('y');
  title('Control: Real Part');
  axis square
  hold on
  h = plot3(xunit, yunit, 2*ones(size(xunit)), 'k-', 'LineWidth', 2);
  hold off

  figure, trisurf(adj, nodes(:,1), nodes(:,2), imag(z));
  shading interp;
  view(0,90)
  axis('equal','tight');
  xlabel('x');
  ylabel('y');
  title('Control: Imaginary Part');
  axis square
  hold on
  h = plot3(xunit, yunit, 2*ones(size(xunit)), 'k-', 'LineWidth', 2);
  hold off

end
