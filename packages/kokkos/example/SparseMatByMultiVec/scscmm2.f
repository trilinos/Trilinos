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
      real*8 val(*), x(*), y(*)
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

      chunk = 4
      rem = mod(nrhs,chunk)
c
c.....do a series of SPAXPYs (sparse saxpys)
      do 20 j = 1, n
C         write(*,*)'In scscmm2 Column ',j
         jbgn = pntr(j)
         jend = pntr(j+1) - 1
         do 30 i = jbgn, jend
            indxi = indx(i)
            vali = val(i)
C            write(*,*)'In scscmm2 val(',i,') = ',val(i)
C            write(*,*)'In scscmm2 indx(',i,') = ',indx(i)
            incn = 0
            do 40 k=1, nrhs
               y(indxi+incn) = y(indxi+incn) + vali * x(j+incn)
               incn = incn + n
 40         continue
 30      continue
 20   continue
C      do 11 i=1,n
C         write(*,*)'In scscmm2 y(',i,') = ',y(i)
C 11   continue
      return
      end
