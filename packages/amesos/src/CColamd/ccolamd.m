function [p,stats] = ccolamd (S, cset, knobs) 
%COLAMD Column approximate minimum degree permutation.
%    P = COLAMD (S) returns the column approximate minimum degree permutation
%    vector for the sparse matrix S.  For a non-symmetric matrix S, S (:,P)
%    tends to have sparser LU factors than S.  The Cholesky factorization of
%    S (:,P)' * S (:,P) also tends to be sparser than that of S'*S.  COLAMD
%    tends to be faster than COLMMD and tends to return a better ordering.
%
%    See also COLMMD, COLPERM, SPPARMS, SYMAMD, SYMMMD, SYMRCM.
%
%    Usage:  P = colamd (S)
%            P = colamd (S, knobs)
%            [P, stats] = colamd (S)
%            [P, stats] = colamd (S, knobs)
%
%    knobs is an optional two-element input vector.  If S is m-by-n, then
%    rows with more than (knobs (1))*n entries are ignored.  Columns with
%    more than (knobs (2))*m entries are removed prior to ordering, and
%    ordered last in the output permutation P.  If the knobs parameter is not
%    present, then 0.5 is used instead, for both knobs (1) and knobs (2). 
%    knobs (3) controls the printing of statistics and error messages.
%
%    stats is an optional 20-element output vector that provides data about the
%    ordering and the validity of the input matrix S.  Ordering statistics are
%    in stats (1:3).  stats (1) and stats (2) are the number of dense or empty
%    rows and columns ignored by COLAMD and stats (3) is the number of
%    garbage collections performed on the internal data structure used by
%    COLAMD (roughly of size 2.2*nnz(S) + 4*m + 7*n integers).
%
%    MATLAB built-in functions are intended to generate valid sparse matrices,
%    with no duplicate entries, with ascending row indices of the nonzeros
%    in each column, with a non-negative number of entries in each column (!)
%    and so on.  If a matrix is invalid, then COLAMD may or may not be able
%    to continue.  If there are duplicate entries (a row index appears two or
%    more times in the same column) or if the row indices in a column are out
%    of order, then COLAMD can correct these errors by ignoring the duplicate
%    entries and sorting each column of its internal copy of the matrix S (the
%    input matrix S is not repaired, however).  If a matrix is invalid in other
%    ways then COLAMD cannot continue, an error message is printed, and no
%    output arguments (P or stats) are returned.  COLAMD is thus a simple way
%    to check a sparse matrix to see if it's valid.
%
%    stats (4:7) provide information if COLAMD was able to continue.  The
%    matrix is OK if stats (4) is zero, or 1 if invalid.  stats (5) is the
%    rightmost column index that is unsorted or contains duplicate entries,
%    or zero if no such column exists.  stats (6) is the last seen duplicate
%    or out-of-order row index in the column index given by stats (5), or zero
%    if no such row index exists.  stats (7) is the number of duplicate or
%    out-of-order row indices.
%
%    stats (8:20) is always zero in the current version of COLAMD (reserved
%    for future use).
%
%    The ordering is followed by a column elimination tree post-ordering.
%
% CCOLAMD version 0.1, May 13, 2005
% Copyright (C) 2005, Univ. of Florida.  Authors: Timothy A. Davis
% and Sivasankaran Rajamanickam.  Closely based on COLAMD by Davis,
% Stefan Larimore, Esmond Ng, and John Gilbert.
% See License.txt for the Version 2.1 of the GNU Lesser General Public License
% http://www.cise.ufl.edu/research/sparse

error ('ccolamd: mexFunction not found') ;

%-------------------------------------------------------------------------------
% Perform the colamd ordering:
%-------------------------------------------------------------------------------

if (nargout <= 1 & nargin == 2)
    p = ccolamdmex (S, cset) ;
elseif (nargout <= 1 & nargin == 3)
    p = ccolamdmex (S, cset, knobs) ;
elseif (nargout == 2 & nargin == 2)
    [p, stats] = ccolamdmex (S, cset) ;
elseif (nargout == 2 & nargin == 3)
    [p, stats] = ccolamdmex (S, cset, knobs) ;
else
    error ('ccolamd:  incorrect number of input and/or output arguments') ;
end

%-------------------------------------------------------------------------------
% column elimination tree post-ordering:
%-------------------------------------------------------------------------------

% Commented for ccolamd
%[ignore, q] = sparsfun ('coletree', S (:,p)) ;
%p = p (q) ;

