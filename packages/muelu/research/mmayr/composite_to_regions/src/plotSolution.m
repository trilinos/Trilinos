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
    [X,Y] = meshgrid(1:7,1:7);
    Z = zeros(7,7);
    for j = 1:7 % rows of the mesh
      for i = 1:7 % cols of the mesh
        Z(j,i) = sol(7*(j-1) + i);
      end
    end
    surf(X,Y,Z);
  otherwise
    error("Can't plot %dD problems. Needs to be implemented.", dim);
end

end