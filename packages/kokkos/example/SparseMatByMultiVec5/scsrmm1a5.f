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
      real*8 vali, yj1, yj2, yj3, yj4, yj5

*
*     --------------------------
*     First executable statement
*     --------------------------
*

      if (nrhs.eq.1) then

      do 110 j = 1, n
         jbgn = pntr(j)
         jend = pntr(j+1) - 1
            yj1 = 0.0
            do 120 i = jbgn, jend
               indxi = (indx(i)-1)*nrhs
               yj1 = yj1 + val(i) * x(indxi+1)
 120        continue
            y(j) = yj1
 110  continue

      else if (nrhs.eq.2) then

      do 210 j = 1, n
         jbgn = pntr(j)
         jend = pntr(j+1) - 1
         jj = (j-1)*nrhs
            yj1 = 0.0
            yj2 = 0.0
            do 220 i = jbgn, jend
               indxi = (indx(i)-1)*nrhs
               yj1 = yj1 + val(i) * x(indxi+1)
               yj2 = yj2 + val(i) * x(indxi+2)
 220        continue
            y(jj+1) = yj1
            y(jj+2) = yj2
 210  continue

      else if (nrhs.eq.3) then

      do 310 j = 1, n
         jbgn = pntr(j)
         jend = pntr(j+1) - 1
         jj = (j-1)*nrhs
            yj1 = 0.0
            yj2 = 0.0
            yj3 = 0.0
            do 320 i = jbgn, jend
               indxi = (indx(i)-1)*nrhs
               yj1 = yj1 + val(i) * x(indxi+1)
               yj2 = yj2 + val(i) * x(indxi+2)
               yj3 = yj3 + val(i) * x(indxi+3)
 320        continue
            y(jj+1) = yj1
            y(jj+2) = yj2
            y(jj+3) = yj3
 310  continue

      else if (nrhs.eq.4) then

      do 410 j = 1, n
         jbgn = pntr(j)
         jend = pntr(j+1) - 1
         jj = (j-1)*nrhs
            yj1 = 0.0
            yj2 = 0.0
            yj3 = 0.0
            yj4 = 0.0
            do 420 i = jbgn, jend
               indxi = (indx(i)-1)*nrhs
               yj1 = yj1 + val(i) * x(indxi+1)
               yj2 = yj2 + val(i) * x(indxi+2)
               yj3 = yj3 + val(i) * x(indxi+3)
               yj4 = yj4 + val(i) * x(indxi+4)
 420        continue
            y(jj+1) = yj1
            y(jj+2) = yj2
            y(jj+3) = yj3
            y(jj+4) = yj4
 410  continue

      else if (nrhs.eq.5) then

      do 510 j = 1, n
         jbgn = pntr(j)
         jend = pntr(j+1) - 1
         jj = (j-1)*nrhs
            yj1 = 0.0
            yj2 = 0.0
            yj3 = 0.0
            yj4 = 0.0
            yj5 = 0.0
            do 520 i = jbgn, jend
               indxi = (indx(i)-1)*nrhs
               yj1 = yj1 + val(i) * x(indxi+1)
               yj2 = yj2 + val(i) * x(indxi+2)
               yj3 = yj3 + val(i) * x(indxi+3)
               yj4 = yj4 + val(i) * x(indxi+4)
               yj5 = yj5 + val(i) * x(indxi+5)
 520        continue
            y(jj+1) = yj1
            y(jj+2) = yj2
            y(jj+3) = yj3
            y(jj+4) = yj4
            y(jj+5) = yj5
 510  continue

      endif
      return
      end
