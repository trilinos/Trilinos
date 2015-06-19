%Prototype function to generate unsmoothed P from aggregates
function Ptent = matlabPtentTest(A, agg)
  
end

%then call this to compare what MueLu generated and what above generated
function verify(mueluP, matlabP)
  if mueluP == matlabP:
    disp('The prolongators match.');
  else:
    disp('The prolongators do not match.');
  end
end
