% plotSolution.m
%
% Plot the solution.
%
% Input:
%   sol     solution
%   dim     spatial dimension of the problem
%

function [] = plotSolution(sol, dim)

% create figure
figure(1);

switch (dim)
  case 1
    plot(sol);
  case 2
    nx = sqrt(size(sol,1));
    
    if (nx^2 ~= size(sol,1))
      error('Missmatch in dimensions.');
    end
    
    [X,Y] = meshgrid(1:nx,1:nx);
    Z = zeros(nx,nx);
    for j = 1:nx % rows of the mesh
      for i = 1:nx % cols of the mesh
        Z(j,i) = sol(nx*(j-1) + i);
      end
    end
    surf(X,Y,Z);
  otherwise
    error("Can't plot %dD problems. Needs to be implemented.", dim);
end

end