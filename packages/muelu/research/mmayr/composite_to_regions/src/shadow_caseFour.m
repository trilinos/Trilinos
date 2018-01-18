% Do a 2-level V-cycle in composite form for caseFour
%
% Perform all operations in composite form to allow for comparison to the region
% algorithm. For simplicity, we skip postsmooting.
%
% Note: Initial guess and forcing vector need to be hard-coded in main.cpp.
%
fineA = mmread('Amat.mm');
finex = [0 1 4 9 0 1 4 0 1 4 9 16 25 36 49 64 0 1 4 0 1 0 1 0 1]';
fineb = zeros(25,1);
fineb(25) = 1e-3;
findd = diag(fineA);
omega = .67;
load P ; P(:,1) = P(:,1)+1; P(:,2) = P(:,2)+1; P = spconvert(P);
Pcombo = P;
Pcombo(:, 6+1) = Pcombo(:, 6+1) + Pcombo(:,26+1);
Pcombo(:,15+1) = Pcombo(:,15+1) + Pcombo(:,27+1);
Pcombo(26:end,:) = 0;
Pcombo(:,26:end) = 0;
nonEmptyRows = find(Pcombo *ones(28,1));
nonEmptyCols = find(Pcombo'*ones(28,1));
Phat = Pcombo(nonEmptyRows,nonEmptyCols);
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
Phat*coarsex,



