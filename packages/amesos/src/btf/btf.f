C===============================================================================
C=== btf =======================================================================
C===============================================================================

	subroutine btf (Ai0, Ap0, nz, n, Ap, Ai, Perm, Lenc, Cp, Bp,
     *	    Cperm, Rperm, Work, nblocks)
	integer nz, n
	integer Ai (nz), Ap (n+1), Ai0 (nz), Ap0 (n+1), Perm (n), Lenc (n),
     *	    Work (4*n), Cp (n), Bp (n), Cperm (n), Rperm (n)
	integer p, col, nzdiag, k, newcol, oldcol

c	copy the matrix and convert to 1-based indexing
	do 10 col = 1,n+1
	    Ap (col) = Ap0 (col) + 1
10	continue
	do 20 p = 1,nz
	    Ai (p) = Ai0 (p) + 1
20	continue

c	compute the lengths of each column
	do 30 col = 1,n
	    Lenc (col) = Ap (col+1) - Ap (col)
30	continue

cc	print the matrix
c	print *, n, nz
c	do 130 col = 1,n
c	    print *,' col: ', col, ' len: ', Lenc (col)
c	    do 140 p = Ap (col), Ap (col+1)-1
c		print *, '    row: ', Ai (p)
c140	    continue
c130	continue

c	get a zero-free diagonal
	call mc21a (n, Ai, nz, Ap, Lenc, Perm, nzdiag, Work)

	if (nzdiag .ne. n) then
c	    matrix is structurally singular - do not continue
	    nblocks = 0
	    return
	endif

c	Perm (newcol) = oldcol is now the zero-free permutation.  In MATLAB
c	notation, A (:,Perm) has a zero-free diagonal

c	permute the columns of A.  Use Cp as workspace for column pointers
	do 50 newcol = 1, n
	    oldcol = Perm (newcol)
	    Cp (newcol) = Ap (oldcol)
	    Lenc (newcol) = Ap (oldcol+1) - Ap (oldcol)
50	continue

c	permutation to block-triangular-form
	call mc13d (n, Ai, nz, Cp, Lenc, Rperm, Bp, nblocks, Work)

c	find the composite column permutation
	do 60 k = 1, n
	    Cperm (k) = Perm (Rperm (k))
60	continue

c	do 110 k = 1,n
c	    print *, 'k:', k, ' P:', Rperm (k), ' Q:', Cperm (k)
c110	continue
c	do 111 k = 1,nblocks
c	    print *, 'block:', k, ' Bp:', Bp (k)
c111	continue

	return
	end


