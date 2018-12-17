% Do a 2-level V-cycle in composite form for caseFifteen
%
% Perform all operations in composite form to allow for comparison to the region
% algorithm. For simplicity, we skip postsmooting.
%
% Note: 
% - Initial guess and forcing vector need to be hard-coded in main.cpp.
% - Disable QR in tentative P
%
fineA = twoDimensionalLaplace(49);
finex = [0 1 4 9 16 25 36 49 64 81 100 121 144 169 196 225 0 1 4 9 16 25 36 49 64 81 100 121 0 1 4 9 16 25 36 49 64 81 100 121 0 1 4 9 16 25 36 49 64]';
fineb = zeros(size(finex));
fineb(end) = 1e-3;
findd = diag(fineA);
omega = .67;

Phat = [1 0 0 0 0   0 0 0 0; % 0
        1 0 0 0 0   0 0 0 0; % 1
        0 1 0 0 0   0 0 0 0; % 2
        0 1 0 0 0   0 0 0 0; % 3
        0 1 0 0 0   0 0 0 0; % 4
        0 0 1 0 0   0 0 0 0; % 5
        0 0 1 0 0   0 0 0 0; % 6
        1 0 0 0 0   0 0 0 0; % 7
        1 0 0 0 0   0 0 0 0; % 8
        0 1 0 0 0   0 0 0 0; % 9
        0 1 0 0 0   0 0 0 0; % 10
        0 1 0 0 0   0 0 0 0; % 11
        0 0 1 0 0   0 0 0 0; % 12
        0 0 1 0 0   0 0 0 0; % 13
        0 0 0 1 0   0 0 0 0; % 14
        0 0 0 1 0   0 0 0 0; % 15
        0 0 0 0 1   0 0 0 0; % 16
        0 0 0 0 1   0 0 0 0; % 17
        0 0 0 0 1   0 0 0 0; % 18
        0 0 0 0 0   1 0 0 0; % 19
        0 0 0 0 0   1 0 0 0; % 20
        0 0 0 1 0   0 0 0 0; % 21
        0 0 0 1 0   0 0 0 0; % 22
        0 0 0 0 1   0 0 0 0; % 23
        0 0 0 0 1   0 0 0 0; % 24
        0 0 0 0 1   0 0 0 0; % 25
        0 0 0 0 0   1 0 0 0; % 26
        0 0 0 0 0   1 0 0 0; % 27
        0 0 0 1 0   0 0 0 0; % 28
        0 0 0 1 0   0 0 0 0; % 29
        0 0 0 0 1   0 0 0 0; % 30
        0 0 0 0 1   0 0 0 0; % 31
        0 0 0 0 1   0 0 0 0; % 32
        0 0 0 0 0   1 0 0 0; % 33
        0 0 0 0 0   1 0 0 0; % 34
        0 0 0 0 0   0 1 0 0; % 35
        0 0 0 0 0   0 1 0 0; % 36
        0 0 0 0 0   0 0 1 0; % 37
        0 0 0 0 0   0 0 1 0; % 38
        0 0 0 0 0   0 0 1 0; % 39
        0 0 0 0 0   0 0 0 1; % 40
        0 0 0 0 0   0 0 0 1; % 41
        0 0 0 0 0   0 1 0 0; % 42
        0 0 0 0 0   0 1 0 0; % 43
        0 0 0 0 0   0 0 1 0; % 44
        0 0 0 0 0   0 0 1 0; % 45
        0 0 0 0 0   0 0 1 0; % 46
        0 0 0 0 0   0 0 0 1; % 47
        0 0 0 0 0   0 0 0 1;]; % 48
    
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



