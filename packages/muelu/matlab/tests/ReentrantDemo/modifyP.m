%Function to demonstrate reentrant capability in muemex
%Take a P from a default MueLu hierarchy and modify it
%so that all nonzero values only have one significant figure
function P = createP(A)
  %Create a separate hierarchy to get a fresh, default P
  disp('Generating inner hierarchy...');
  newHier = muelu('setup', A, 'coarse: max size', 100, 'verbosity', 'none');
  P = muelu('get', newHier, 1, 'P'); %This always gives the next prolongator, no matter what the current level # is
  disp('Done generating inner hierarchy.');
  muelu('cleanup', newHier); %If not cleaned up, the inner hierarchy will stay in list of problems along with outer one
  [i, j, s] = find(P);
  for k = 1:length(s)
    placeVal = 10^floor(log10(s(k))); %If P(i) is 0.0513, placeVal will be 0.01.
    s(k) = placeVal * ceil(s(k) / placeVal); %Dividing by placeVal shifts the decimal representation so that exactly one digit is left of decimal place.
  end
  P = sparse(i, j, s);
end 
