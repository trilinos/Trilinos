function [L,U,P] = ilu_left_schur (S, domain, partr)
%ILU_LEFT_SCHUR left-looking ILU factorization.
% Use for factoring the Schur complement. Should be used after LUfactor.
% 
% Example:
%   [L,U,P] = ilu_left_schur (S, domain, partr)
%
%   S - input matrix, (if LUdata = LUfactor(A, k); then S = LUdata.S{k+1})
%   domain - connectivity of the partitions (LUdata.domain)
%   partr  - partition numbers for each of the rows of S (LUdata.Spartr)
%
%

% Based on Tim Davis's lu_left in CSparse

A = S ;
n = size (A,1) ;
P = eye (n) ;
L = zeros (n) ;
U = zeros (n) ;
for k = 1:n
    x = [ L(:,1:k-1) [ zeros(k-1,n-k+1) ; eye(n-k+1) ]] \ (P * A (:,k));
    [a i] = max (abs (x (k:n)));        % find the pivot row i
    i = i + k - 1 ;

    % subdomains connected to partr(i) (subdomain of the pivot row)
    di = find(domain(:, partr(i)));
    pxi = partr(find(x));     % subdomain of each non zero.
    dropd = setdiff(pxi, di); % subdomains in x that need to be dropped
    dropi = [];               % indices to drop 
    % workaround bug in setdiff, dindex is 0-by-1 and empty if no match in 
    % pxi and di
    if (~isempty(dropd))
        for dindex = dropd'
            % TODO: using partr instead of pxi for ease.
            dummy = find(partr == dindex);
            dropi = [dropi; dummy];
        end
    end
    x(dropi) = 0.0;

    U (1:k-1,k) = x (1:k-1) ;           % the column of U
    L ([i k],:) = L ([k i], :) ;        % swap rows i and k of L, P, and x
    P ([i k],:) = P ([k i], :) ;
    x ([i k]) = x ([k i]) ;
    partr ([i k]) = partr([k i]);
    U (k,k) = x (k) ;
    L (k,k) = 1 ;
    L (k+1:n,k) = x (k+1:n) / x (k) ;   % divide the pivot column by U(k,k)
end
