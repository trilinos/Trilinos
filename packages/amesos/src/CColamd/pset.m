function pset (cset)
% print a cset


% CCOLAMD version 0.1, May 13, 2005
% Copyright (C) 2005, Univ. of Florida.  Authors: Timothy A. Davis
% and Sivasankaran Rajamanickam
% See License.txt for the Version 2.1 of the GNU Lesser General Public License
% http://www.cise.ufl.edu/research/sparse


nsets = max (cset) ;
n = length (cset)

for k = 1:nsets
    s = find (cset == k) ;
    fprintf ('set %2d: size: %2d [', k, length (s)) ;
    fprintf (' %2d', s) ;
    fprintf (']\n') ;
end

