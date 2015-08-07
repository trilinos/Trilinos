%Prototype function to test evolution strength-of-connection
%
% reference: "A new perspective on strength measures in algebraic multigrid"
%            Luke Olson, Jacob Schroder, and Ray Tuminaro
%            Numerical Linear Algebra with Applications, Vol. 17, p 713-733, 2010
%            doi 10.1002/nla.669
%
% Notes:  Dofs per node is hardwired to 1
%
function [Graph, DofsPerNode] = evolutionSoC(A)
  fprintf('entering evolutionSoC\n');

  % for now, just return the matrix graph to see that it works
  Graph.edges = logical(A);
  Graph.boundaryNodes = int32([]);
  DofsPerNode = 1;

  fprintf('leaving evolutionSoC\n');
end 
