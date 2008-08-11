      subroutine scscmm2( m, n, val, indx, pntr, x, y, nrhs)
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
      real*8 val(*), x(*), y(*), xj
*
*     ----------------------------------
*     Specifications for local variables
*     ----------------------------------
      integer i,j,jbgn, jend, indxi, rem, chunk, k, incn
      real*8 vali
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
C         write(*,*)'In scscmm2 Column ',j
         jbgn = pntr(j)
         jend = pntr(j+1) - 1
         incn = 0
         do 40 k=1, nrhs
            xj = x(j+incn)
            do 30 i = jbgn, jend
               y(indx(i)+incn) = y(indx(i)+incn) + val(i) * xj
 30         continue
            incn = incn + n
 40      continue
 20   continue
C      do 11 i=1,n
C         write(*,*)'In scscmm2 y(',i,') = ',y(i)
C 11   continue
      return
      end
