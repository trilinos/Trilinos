*----------------------------------------------------------------------
*
*
*\Documentation
*
*\Name: KOKKOS_DCRSMM
*
*\Description:
*	KOKKOS_DCRSMM performs one of the matrix-multivector operations
*
*    y = A*x or y = A'*x.
*
*    where x, y are vectors, A is an
*    m-by-n matrix, and A' is the transpose of A.  The matrix
*    A is stored in a special one dimensional column-pointer
*    format where the zeroes in each column are discarded and
*    the non-zeroes are stored contiguously.  A corresponding
*    integer array, pntr, holds pointers indicating the start of
*    each row in A.  Each element of the array rowind contains
*    the row index of the corresponding element of A.
*
*\Usage:
*     call DCRSMM(  itrans, udiag, m, n, val, indx, profile, x, ldx, y, ldy, nrhs )
*
*
*    itrans Integer (input)
*           On entry, itrans specifies the operation to be performed.
*           If itrans = 0, z := A*x + beta*y.
*           If itrans = 1, z := A'*x + beta*y.
*           If itrans is any other value, then no operation is
*           performed.
*           The itrans argument is unchanged on exit.
*
*    udiag  Integer (input)
*           On entry, udiag specifies whether or not the matrix should
*           be assumed to have a unit diagonal.
*           If udiag = 1, add x to the result of y as though unit 
*           diagonal were present.
*           If udiag is any other value, then no unit diagonal is assumed.
*           The udiag argument is unchanged on exit.
*
*    m      Integer (input)
*           On entry, m specifies the number of rows of
*           the matrix a.  m must be at least 0.  The m argument
*           is unchanged on exit.
*
*    n      Integer (input)
*           On entry, n specifies the number of columns of
*           the matrix a.  n must be at least 0.  The n argument
*           is unchanged on exit.
*
*    val    real*8 array (input)
*           On entry, val holds the values of matrix A in packed form as
*           described above.
*           The array val is unchanged on exit.
*
*    indx   Integer array (input)
*           On entry, indx holds the row indices of the non-zero
*           elements in A.  indx must have length > = nnz.
*           The array indx is unchanged on exit.
*
*    profile Integer array (input)
*           On entry, profile(j) contains the number of entries in the jth row
*           profile must have length > = n.
*           profile is unchanged on exit.
*
*    x      real*8 array (input)
*           Real array of dimension at least n by nrhs when itrans = 0
*           and at least m by nrhs otherwise.  Before entry, the array x
*           must contain the vector operand x.
*           The argument x is unchanged on exit.
*
*    y      real*8 array (output)
*           Real array of dimension at least m by nrhs when itrans = 0 and
*           at least n by nrhs otherwise.  On exit, y contains the result
*           vector.
*
*\Remarks:
*    1.  x and y cannot be the same vector.
*        Unpredictable results will occur if x and y are the same.
*
*    2.  Although the example below stores the elements of each
*        column in natural order, this routine makes no assumption
*        about the order of the non-zeroes within a column.
*
*\Examples:
*    If the original matrix is
*
*                   | 11  12   0  14   0 |
*                   | 21  22  23   0   0 |
*                   |  0  32  33   0   0 |
*                   |  0   0   0   0  45 |
*                   | 51  52   0  54   0 |
*
*    then the matrix is assumed to be store as
*
*    a = ( 11  21  51  12  22  32  52  23  33  14  54  45 )
*
*    with the corresponding pointer arrays
*
*   indx = (  1   2   5   1   2   3   5   2   3   1   5   4 ).
*
*                profile = (3  3  2  1  3)
*
*    Thus, indx(i) indicates the row position of the ith element
*    of a and pntr(j) points to the first element of jth column.
*
*\Enddoc

      subroutine kokkos_dcrsmm( itrans, udiag, m, n, val, indx, profile, 
     &                          x, ldx, y, ldy, nrhs)

*     ----------------------------
*     Specifications for arguments
*     ----------------------------
      implicit none
      integer itrans, udiag, m, n, ldx, ldy
      integer nrhs, indx(0:*), profile(0:*)
      real*8 val(0:*), x(0:*), y(0:*)
*
*     ----------------------------------
*     Specifications for local variables
*     ----------------------------------
      integer irhs, nrhs1, ix, iy

*     Strip mine nrhs by fives

      nrhs1 = mod(nrhs, 5)
      if (nrhs1.eq.0) nrhs1 = 5
      ix = 0
      iy = 0
      do 10 irhs = 1, nrhs, 5
         if (itrans.eq.0) then
            call kokkos_scrsmm5(udiag, m, n, val, indx, profile,
     &                          x(ix), ldx, y(iy), ldy, nrhs1)
         else
            call kokkos_sccsmm5(udiag, m, n, val, indx, profile,
     &                          x(ix), ldx, y(iy), ldy, nrhs1)
         endif
         ix = ix + nrhs1*ldx
         iy = iy + nrhs1*ldy
         nrhs1 = nrhs1 - 5
 10   continue

      return
      end

      subroutine kokkos_sccsmm5( udiag, m, n, val, indx, profile, 
     &                           x, ldx, y, ldy, nrhs)
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
      integer udiag, m, n, nrhs, ldx, ldy
      integer indx(0:*), profile(0:*)
      real*8 val(0:*), x(0:*), y(0:*)
*
*     ----------------------------------
*     Specifications for local variables
*     ----------------------------------
      integer i, j, k, jbgn, jend, incx, incy
      integer iy, ix
      real*8 vali, xj1, xj2, xj3, xj4, xj5

*
*     --------------------------
*     First executable statement
*     --------------------------
*
c.....initialize soln
*
*     Implicit unit diagonal

      if (udiag.eq.1) then
         if (ldy.eq.n.and.ldx.eq.m.and.m.eq.n) then
            do 1 k=0, n*nrhs-1
               y(k) = x(k)
 1          continue
         else
            do 5 k=1, nrhs
               iy = ldy*(k-1)
               ix = ldx*(k-1)
               do 10 i = 0, min(n-1, m-1)
                  y(i+iy) = x(i+ix)
 10            continue
 5          continue
         endif

* No implicit unit diagonal

      else
         if (ldy.eq.n) then
            do 11 k=0, n*nrhs-1
               y(k) = 0.0
 11         continue
         else
            do 15 k=1, nrhs
               iy = ldy*(k-1)
               do 20 i = 0, n-1
                  y(i+iy) = 0.0D0
 20            continue
 15         continue
         endif
      endif
c
c.....do a series of SPAXPYs (sparse saxpys)
      if (nrhs.eq.1) then
         jend = 0
         do 120 j = 0, n-1
            jbgn = jend
            jend = jbgn + profile(j)
            xj1 = x(j)
            do 130 i = jbgn, jend-1
               y(indx(i)) = y(indx(i)) + val(i) * xj1
 130        continue
 120     continue
         
      else if (nrhs.eq.2) then
         
      jend = 0
      do 220 j = 0, n-1
         jbgn = jend
         jend = jbgn + profile(j)
         incx = j
         xj1 = x(incx)
         incx = incx + ldx
         xj2 = x(incx)
         do 230 i = jbgn, jend-1
            incy = indx(i)
            vali = val(i)
            y(incy) = y(incy) + vali * xj1
            incy = incy + ldy
            y(incy) = y(incy) + vali * xj2
 230     continue
 220  continue
         
      else if (nrhs.eq.3) then
         
      jend = 0
      do 320 j = 0, n-1
         jbgn = jend
         jend = jbgn + profile(j)
         incx = j
         xj1 = x(incx)
         incx = incx + ldx
         xj2 = x(incx)
         incx = incx + ldx
         xj3 = x(incx)
         do 330 i = jbgn, jend-1
            incy = indx(i)
            vali = val(i)
            y(incy) = y(incy) + vali * xj1
            incy = incy + ldy
            y(incy) = y(incy) + vali * xj2
            incy = incy + ldy
            y(incy) = y(incy) + vali * xj3
 330     continue
 320  continue
      
      else if (nrhs.eq.4) then
         
      jend = 0
      do 420 j = 0, n-1
         jbgn = jend
         jend = jbgn + profile(j)
         incx = j
         xj1 = x(incx)
         incx = incx + ldx
         xj2 = x(incx)
         incx = incx + ldx
         xj3 = x(incx)
         incx = incx + ldx
         xj4 = x(incx)
         do 430 i = jbgn, jend-1
            incy = indx(i)
            vali = val(i)
            y(incy) = y(incy) + vali * xj1
            incy = incy + ldy
            y(incy) = y(incy) + vali * xj2
            incy = incy + ldy
            y(incy) = y(incy) + vali * xj3
            incy = incy + ldy
            y(incy) = y(incy) + vali * xj4
 430     continue
 420  continue
      
      else if (nrhs.eq.5) then
         
      jend = 0
      do 520 j = 0, n-1
         jbgn = jend
         jend = jbgn + profile(j)
         incx = j
         xj1 = x(incx)
         incx = incx + ldx
         xj2 = x(incx)
         incx = incx + ldx
         xj3 = x(incx)
         incx = incx + ldx
         xj4 = x(incx)
         incx = incx + ldx
         xj5 = x(incx)
         do 530 i = jbgn, jend-1
            incy = indx(i)
            vali = val(i)
            y(incy) = y(incy) + vali * xj1
            incy = incy + ldy
            y(incy) = y(incy) + vali * xj2
            incy = incy + ldy
            y(incy) = y(incy) + vali * xj3
            incy = incy + ldy
            y(incy) = y(incy) + vali * xj4
            incy = incy + ldy
            y(incy) = y(incy) + vali * xj5
 530     continue
 520  continue
      
      endif
      
      return
      end

      subroutine kokkos_scrsmm5( udiag, m, n, val, indx, profile, 
     &                           x, ldx, y, ldy, nrhs)
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
      integer udiag, m, n, nrhs, ldx, ldy
      integer indx(0:*), profile(0:*)
      real*8 val(0:*), x(0:*), y(0:*)
*
*     ----------------------------------
*     Specifications for local variables
*     ----------------------------------
      integer i,j,jbgn, jend, incx, incy
      real*8 vali, yj1, yj2, yj3, yj4, yj5

*
*     --------------------------
*     First executable statement
*     --------------------------
*
      if (nrhs.eq.1) then
         jend = 0
         do 110 j = 0, m-1
            jbgn = jend
            jend = jbgn + profile(j)
            yj1 = 0.0
            do 120 i = jbgn, jend-1
               yj1 = yj1 + val(i) * x(indx(i))
 120        continue
            if (udiag.eq.1) then
               y(j) = yj1 + x(j)
            else
               y(j) = yj1
            endif
 110     continue
         

      else if (nrhs.eq.2) then

         jend = 0
         do 210 j = 0, m-1
            jbgn = jend
            jend = jbgn + profile(j)
            yj1 = 0.0
            yj2 = 0.0
            do 220 i = jbgn, jend-1
               vali = val(i)
               incx = indx(i)
               yj1 = yj1 + vali * x(incx)
               incx = incx + ldx
               yj2 = yj2 + vali * x(incx)
 220        continue
            if (udiag.eq.1) then
               incy = j
               incx = j
               y(incy) = yj1 + x(incy)
               incy = incy + ldy
               incx = incx + ldx
               y(incy) = yj2 + x(incx)
            else
               incy = j
               y(incy) = yj1
               incy = incy + ldy
               y(incy) = yj2
            endif
 210     continue
         
      else if (nrhs.eq.3) then

         jend = 0
         do 310 j = 0, m-1
            jbgn = jend
            jend = jbgn + profile(j)
            incx = 0
            yj1 = 0.0
            yj2 = 0.0
            yj3 = 0.0
            do 320 i = jbgn, jend-1
               vali = val(i)
               incx = indx(i)
               yj1 = yj1 + vali * x(incx)
               incx = incx + ldx
               yj2 = yj2 + vali * x(incx)
               incx = incx + ldx
               yj3 = yj3 + vali * x(incx)
 320        continue
            if (udiag.eq.1) then
               incy = j
               incx = j
               y(incy) = yj1 + x(incy)
               incy = incy + ldy
               incx = incx + ldx
               y(incy) = yj2 + x(incx)
               incy = incy + ldy
               incx = incx + ldx
               y(incy) = yj3 + x(incx)
            else
               incy = j
               y(incy) = yj1
               incy = incy + ldy
               y(incy) = yj2
               incy = incy + ldy
               y(incy) = yj3
            endif
 310     continue

      else if (nrhs.eq.4) then

         jend = 0
         do 410 j = 0, m-1
            jbgn = jend
            jend = jbgn + profile(j)
            yj1 = 0.0
            yj2 = 0.0
            yj3 = 0.0
            yj4 = 0.0
            do 420 i = jbgn, jend-1
               vali = val(i)
               incx = indx(i)
               yj1 = yj1 + vali * x(incx)
               incx = incx + ldx
               yj2 = yj2 + vali * x(incx)
               incx = incx + ldx
               yj3 = yj3 + vali * x(incx)
               incx = incx + ldx
               yj4 = yj4 + vali * x(incx)
 420        continue
            if (udiag.eq.1) then
               incy = j
               incx = j
               y(incy) = yj1 + x(incy)
               incy = incy + ldy
               incx = incx + ldx
               y(incy) = yj2 + x(incx)
               incy = incy + ldy
               incx = incx + ldx
               y(incy) = yj3 + x(incx)
               incy = incy + ldy
               incx = incx + ldx
               y(incy) = yj4 + x(incx)
            else
               incy = j
               y(incy) = yj1
               incy = incy + ldy
               y(incy) = yj2
               incy = incy + ldy
               y(incy) = yj3
               incy = incy + ldy
               y(incy) = yj4
            endif
 410     continue

      else if (nrhs.eq.5) then

         jend = 0
         do 510 j = 0, m-1
            jbgn = jend
            jend = jbgn + profile(j)
            yj1 = 0.0
            yj2 = 0.0
            yj3 = 0.0
            yj4 = 0.0
            yj5 = 0.0
            do 520 i = jbgn, jend-1
               vali = val(i)
               incx = indx(i)
               yj1 = yj1 + vali * x(incx)
               incx = incx + ldx
               yj2 = yj2 + vali * x(incx)
               incx = incx + ldx
               yj3 = yj3 + vali * x(incx)
               incx = incx + ldx
               yj4 = yj4 + vali * x(incx)
               incx = incx + ldx
               yj5 = yj5 + vali * x(incx)
 520        continue
            if (udiag.eq.1) then
               incy = j
               incx = j
               y(incy) = yj1 + x(incy)
               incy = incy + ldy
               incx = incx + ldx
               y(incy) = yj2 + x(incx)
               incy = incy + ldy
               incx = incx + ldx
               y(incy) = yj3 + x(incx)
               incy = incy + ldy
               incx = incx + ldx
               y(incy) = yj4 + x(incx)
               incy = incy + ldy
               incx = incx + ldx
               y(incy) = yj5 + x(incx)
            else
               incy = j
               y(incy) = yj1
               incy = incy + ldy
               y(incy) = yj2
               incy = incy + ldy
               y(incy) = yj3
               incy = incy + ldy
               y(incy) = yj4
               incy = incy + ldy
               y(incy) = yj5
            endif
 510     continue


      endif
      return
      end
