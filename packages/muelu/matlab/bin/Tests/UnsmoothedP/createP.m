%Prototype function to generate unsmoothed P from aggregates
function [Ptent, Nullspace] = createP(agg)
  Ptent = sparse(double(agg.nVertices), double(agg.nAggregates));
  for i = 1:agg.nVertices
    Ptent(i, 1 + agg.vertexToAggID(i)) = 1;
  end
  Nullspace = ones(1, agg.nVertices);
end 
