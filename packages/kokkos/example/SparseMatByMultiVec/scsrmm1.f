      subroutine scsrmm1( m, n, val, indx, pntr, x, y, nrhs)
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
      integer m, n, nrhs, indx(*), pntr(*)
      real*8 val(*), x(*), y(*)
*
*     ----------------------------------
*     Specifications for local variables
*     ----------------------------------
      integer i,j, jj, jbgn, jend, indxi, rem, chunk, k
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
         jj = (j-1)*nrhs
         do 15 k = 1, nrhs
            y(jj+k) = 0.0
 15      continue
         do 20 i = jbgn, jend
            indxi = (indx(i)-1)*nrhs
            vali = val(i)
C            write(*,*)'In scsrmm1 val(',i,') = ',val(i)
C            write(*,*)'In scsrmm1 indx(',i,') = ',indx(i)
            do 30 k = 1, nrhs
              y(jj+k) = y(jj+k) + vali * x(indxi+k)
 30        continue
 20      continue
 10   continue
C      do 11 i=1,n
C         write(*,*)'In scsrmm1 y(',i,') = ',y(i)
C 11   continue
      return
      end
