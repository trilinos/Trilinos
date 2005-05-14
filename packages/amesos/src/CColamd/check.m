function check (p)
% check the validity of a permutation vector

% CCOLAMD version 0.1, May 13, 2005
% Copyright (C) 2005, Univ. of Florida.  Authors: Timothy A. Davis
% and Sivasankaran Rajamanickam
% See License.txt for the Version 2.1 of the GNU Lesser General Public License
% http://www.cise.ufl.edu/research/sparse


n = length (p) ;
fprintf ('check perm, length: %d\n') ;

pinv = -ones (1,n) ;

for k = 1:n
    i = p (k) ;
    if (i > 0 & i <= n)
	pinv (i) = k ;
    end
end

lost = find (pinv == -1)

lost0 = lost-1
