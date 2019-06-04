% Do a 2-level V-cycle in composite form for caseTen
%
% Perform all operations in composite form to allow for comparison to the region
% algorithm. For simplicity, we skip postsmooting.
%
% Note: 
% - Initial guess and forcing vector need to be hard-coded in main.cpp.
% - Disable QR in tentative P
%
fineA = mmread('Amat.mm');
finex = [0 1 4 9 0 1 4 9 16 25 0 1 4]';
fineb = zeros(size(finex));
fineb(end) = 1e-3;
findd = diag(fineA);
omega = .67;

% load P ; P(:,1) = P(:,1)+1; P(:,2) = P(:,2)+1; P = spconvert(P);
% Pcombo = P;
% Pcombo(:, 6+1) = Pcombo(:, 6+1) + Pcombo(:,26+1);
% Pcombo(:,15+1) = Pcombo(:,15+1) + Pcombo(:,27+1);
% Pcombo(26:end,:) = 0;
% Pcombo(:,26:end) = 0;
% nonEmptyRows = find(Pcombo *ones(28,1));
% nonEmptyCols = find(Pcombo'*ones(28,1));
% Phat = Pcombo(nonEmptyRows,nonEmptyCols);

Phat = [1 0 0 0 0;
        1 0 0 0 0;
        0 1 0 0 0;
        0 1 0 0 0;
        0 1 0 0 0;
        0 0 1 0 0;
        0 0 1 0 0;
        0 0 1 0 0;
        0 0 0 1 0;
        0 0 0 1 0;
        0 0 0 1 0;
        0 0 0 0 1;
        0 0 0 0 1];
    
Rhat = Phat';
AH = Rhat*fineA*Phat;
coarsed = diag(AH);

onlyJacobi = false;

if (onlyJacobi == true)
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
else
  for i=1:1
    finex =  finex + omega*(fineb - fineA*finex)./findd;
    finer = fineb - fineA*finex;
    fprintf('r is %20.13e\n',norm(finer));
    coarseb = Rhat*finer;
    coarsex = AH \ coarseb;
  end
end
Phat*coarsex,



