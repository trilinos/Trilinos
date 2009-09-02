      subroutine scsrmm2( m, n, val, indx, pntr, x, y, nrhs)
*
*     Performs the matrix-vector operation
*
*                               y = A*x
*
*     where x and y are vectors and A is a sparse matrix stored
*     in compress row format.
*
*     ----------------------------
*     Specifications for arguments
*     ----------------------------
      implicit none
      integer m, n, nrhs, indx(*), pntr(m+1)
      real*8 val(*), x(n), y(m)
*
*     ----------------------------------
*     Specifications for local variables
*     ----------------------------------
      integer i,j,jbgn, jend, indxi, k, incn
      real*8 vali

*
*     --------------------------
*     First executable statement
*     --------------------------
*
c.....do sequence of SPDOTs (sparse sdots)
      do 10 j = 1, n
C         write(*,*)'In scsrmm1 Row ',j
         jbgn = pntr(j)
         jend = pntr(j+1) - 1
         incn = 0
         do 15 k = 1, nrhs
            y(j+incn) = 0.0
            incn = incn + n
 15      continue
         do 20 i = jbgn, jend
            indxi = indx(i)
            vali = val(i)
            incn = 0
            do 30 k = 1, nrhs
               y(j+incn) = y(j+incn) + vali * x(indxi+incn)
               incn = incn + n
 30         continue
 20      continue
 10   continue
C      do 11 i=1,n
C         write(*,*)'In scsrmm2 y(',i,') = ',y(i)
C 11   continue
      return
      end
