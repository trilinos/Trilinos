      subroutine scsrmm2( m, n, val, indx, profile, x, y, nrhs)
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
      integer m, n, nrhs, indx(*), profile(m)
      real*8 val(*), x(0:n), y(m)
*
*     ----------------------------------
*     Specifications for local variables
*     ----------------------------------
      integer i,j,jbgn, jend, indxi, incx, incy
      real*8 vali, yj1, yj2, yj3, yj4, yj5

*
*     --------------------------
*     First executable statement
*     --------------------------
*
      jend = 0
      if (nrhs.eq.1) then
      do 110 j = 1, n
         jbgn = jend + 1
         jend = jbgn + profile(j) - 1
         yj1 = 0.0
         do 120 i = jbgn, jend
           yj1 = yj1 + val(i) * x(indx(i))
 120     continue
         y(j) = yj1
 110    continue


      else if (nrhs.eq.2) then

      do 210 j = 1, n
         jbgn = jend + 1
         jend = jbgn + profile(j) - 1
         yj1 = 0.0
         yj2 = 0.0
         do 220 i = jbgn, jend
           indxi = indx(i)
           vali = val(i)
           incx = indx(i)
           yj1 = yj1 + val(i) * x(incx)
           incx = incx + n
           yj2 = yj2 + val(i) * x(incx)
 220     continue
         incy = j
         y(incy) = yj1
         incy = incy + n
         y(incy) = yj2
 210  continue

      else if (nrhs.eq.3) then

      do 310 j = 1, n
         jbgn = jend + 1
         jend = jbgn + profile(j) - 1
         incx = 0
         yj1 = 0.0
         yj2 = 0.0
         yj3 = 0.0
         do 320 i = jbgn, jend
           indxi = indx(i)
           vali = val(i)
           incx = indx(i)
           yj1 = yj1 + val(i) * x(incx)
           incx = incx + n
           yj2 = yj2 + val(i) * x(incx)
           incx = incx + n
           yj3 = yj3 + val(i) * x(incx)
 320     continue
         incy = j
         y(incy) = yj1
         incy = incy + n
         y(incy) = yj2
         incy = incy + n
         y(incy) = yj3
 310  continue

      else if (nrhs.eq.4) then

      do 410 j = 1, n
         jbgn = jend + 1
         jend = jbgn + profile(j) - 1
         yj1 = 0.0
         yj2 = 0.0
         yj3 = 0.0
         yj4 = 0.0
         do 420 i = jbgn, jend
           indxi = indx(i)
           vali = val(i)
           incx = indx(i)
           yj1 = yj1 + val(i) * x(incx)
           incx = incx + n
           yj2 = yj2 + val(i) * x(incx)
           incx = incx + n
           yj3 = yj3 + val(i) * x(incx)
           incx = incx + n
           yj4 = yj4 + val(i) * x(incx)
 420     continue
         incy = j
         y(incy) = yj1
         incy = incy + n
         y(incy) = yj2
         incy = incy + n
         y(incy) = yj3
         incy = incy + n
         y(incy) = yj4
 410  continue

      else if (nrhs.eq.5) then

      do 510 j = 1, n
         jbgn = jend + 1
         jend = jbgn + profile(j) - 1
         yj1 = 0.0
         yj2 = 0.0
         yj3 = 0.0
         yj4 = 0.0
         yj5 = 0.0
         do 520 i = jbgn, jend
           indxi = indx(i)
           vali = val(i)
           incx = indx(i)
           yj1 = yj1 + val(i) * x(incx)
           incx = incx + n
           yj2 = yj2 + val(i) * x(incx)
           incx = incx + n
           yj3 = yj3 + val(i) * x(incx)
           incx = incx + n
           yj4 = yj4 + val(i) * x(incx)
           incx = incx + n
           yj5 = yj5 + val(i) * x(incx)
 520     continue
         incy = j
         y(incy) = yj1
         incy = incy + n
         y(incy) = yj2
         incy = incy + n
         y(incy) = yj3
         incy = incy + n
         y(incy) = yj4
         incy = incy + n
         y(incy) = yj5
 510  continue


      endif
      return
      end
