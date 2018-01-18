% Do a 2-level V-cycle in composite form for caseFour
%
% Perform all operations in composite form to allow for comparison to the region
% algorithm. For simplicity, we skip postsmooting.
%
% Note: Initial guess and forcing vector need to be hard-coded in main.cpp.
%
fineA = mmread('Amat.mm');
finex = [0 1 4 9 16 25 36 49 64 81 0 1 4 9 16 25 36 49 64 0 1 4 9 16 25 36 49 64]';
fineb = zeros(size(finex));
fineb(end) = 1e-3;
findd = diag(fineA);
omega = .67;

% construct 'template P' with plain-aggregation aggregates
P = zeros(length(finex)+2,10);
for i = 1:10
  P(3*i-2,i) = 1.0;
  P(3*i-1,i) = 1.0;
  P(3*i,i) = 1.0;
end

% extract tentative P and form coarse level operator
Phat = P(2:end-1,:);
Rhat = Phat';
AH = Rhat*fineA*Phat;
coarsed = diag(AH);

 
for i=1:1
   finex =  finex + omega*(fineb - fineA*finex)./findd;
   finer = fineb - fineA*finex;
   fprintf('r is %20.13e\n',norm(finer));
   coarseb = Rhat*finer;
   coarsex = zeros(size(coarseb,1),1);
   for k=1:1
      coarsex =  coarsex + omega*(coarseb - AH*coarsex)./coarsed;
   end
   finex   = finex + Phat*coarsex;
end
Phat*coarsex



