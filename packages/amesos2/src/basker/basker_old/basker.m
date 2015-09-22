function [L, U, p] = basker (A, lnnz, unnz) 
% Finds the LU factorization of A such that A=L*U.
% This version doesn't use partial pivoting.
%
% [L, U] = basker(A) will use fully dense L and U to start with.
%
% If nnz of L and U are known then the usage is 
% [L, U] = basker(A, lnnz, unnz)
%

[ m, n ] = size(A) ;
if (nargin == 1)
    %fprintf(' Warning : Assuming maximum size for L and U \n ') ;
    [m, n] = size(A) ;
    annz = nnz(A);
    lnnz = annz;
    %lnnz= ((m*n)/2)/2+ 2*n;
    %lnnz= ((m*n)/2) + 2*n;
    unnz = lnnz;
end

[L, U, p] = baskermex(A, lnnz, unnz) ;

% sort L and U.
L=(L')' ;
U=(U')' ;
p=p' ;

