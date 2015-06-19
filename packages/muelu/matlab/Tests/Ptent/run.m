%Prototype function to generate unsmoothed P from aggregates
function Ptent = matlabPtentTest(A, agg)
  %really only need agg for this
  Ptent = double(sparse(false(agg.nVertices, agg.nAggregates)));
  for i = 1:nVertices
    Ptent(i, agg.vertexToAggID(i)) = 1;
  end
end

%then call this to compare what MueLu generated and what above generated
function verify(mueluP, matlabP)
  if mueluP == matlabP:
    disp('The prolongators match.');
  else:
    disp('The prolongators do not match.');
  end
end
