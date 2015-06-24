%Prototype function to generate unsmoothed P from aggregates
function Ptent = matlabPtentTest(agg)
  Ptent = double(sparse(false(agg.nVertices, agg.nAggregates)));
  for i = 1:nVertices
    Ptent(i, agg.vertexToAggID(i)) = 1;
  end
end 
