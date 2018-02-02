nNodes = 25;
refA = oneDimensionalLaplace(nNodes);
refb = zeros(nNodes,1);
refx = refb;
refb(end) = 1.0;


maxIter = 80;
omega = 0.9;

for iter = 1:maxIter
  refr = refb - refA*refx;
  refx = refx + omega* inv(diag(diag(refA)))*refr;
end
