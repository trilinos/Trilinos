      subroutine scscmm1( m, n, val, indx, pntr, x, y, nrhs)
*
*     Performs the matrix-vector operation
*
*                               y = A*x
*
*     where x and y are vectors and A is a sparse matrix stored
*     in compress column format.
*
*     ----------------------------
*     Specifications for arguments
*     ----------------------------
      implicit none
      integer m, n, nrhs, indx(*), pntr(*)
      real*8 val(*), x(*), y(*)
*
*     ----------------------------------
*     Specifications for local variables
*     ----------------------------------
      integer i,j, jj, jbgn, jend, indxi, rem, chunk, k
      real*8 vali, xj
*
*     --------------------------
*     First executable statement
*     --------------------------
*
c.....initialize soln 
      do 10 i = 1, m*nrhs
         y(i) = 0.0D0
 10   continue
c
c.....do a series of SPAXPYs (sparse saxpys)
      do 20 j = 1, n
C         write(*,*)'In scscmm1 Column ',j
         jbgn = pntr(j)
         jend = pntr(j+1) - 1
         jj = (j-1)*nrhs
         do 40 k=1, nrhs
            xj = x(jj+k)
            do 30 i = jbgn, jend
               indxi = (indx(i)-1)*nrhs
               y(indxi+k) = y(indxi+k) + val(i) * xj
 30         continue
 40      continue
 20   continue
C      do 11 i=1,n
C         write(*,*)'In scscmm1 y(',i,') = ',y(i)
C 11   continue
      return
      end
