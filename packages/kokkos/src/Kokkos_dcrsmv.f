*----------------------------------------------------------------------
*
*
*\Documentation
*
*\Name: KOKKOS_DCRSMV
*
*\Description:
*	KOKKOS_DCRSMV performs one of the matrix-vector operations
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
*     call DCRSMV( itrans, udiag, m, n, val, indx, profile, x, y )
*
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
*           Real array of dimension at least n when itrans = 0
*           and at least m otherwise.  Before entry, the array x
*           must contain the vector operand x.
*           The argument x is unchanged on exit.
*
*    y      real*8 array (output)
*           Real array of dimension at least m when itrans = 0 and
*           at least n otherwise.  On exit, y contains the result
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

      subroutine kokkos_dcrsmv( itrans, udiag, m, n, 
     &                          val, indx, profile, x, y)
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
      integer itrans, udiag, m, n, indx(0:*), profile(0:*)
      real*8 val(0:*), x(0:*), y(0:*)
*
*     ----------------------------------
*     Specifications for local variables
*     ----------------------------------
      integer i,j, jbgn, jend
      real*8 sum

*
*     --------------------------
*     First executable statement
*     --------------------------
*
      if (itrans.eq.0) then
c.....do sequence of SPDOTs (sparse sdots)
         jend = 0
         do 10 j = 0, m-1
            jbgn = jend
            jend = jbgn + profile(j)
            sum = 0.0
            do 20 i = jbgn, jend-1
               sum = sum + val(i) * x(indx(i))
 20         continue
            if (udiag.eq.1) then
               y(j) = sum + x(j)
            else
               y(j) = sum
            endif
 10      continue
      else 
         if (udiag.eq.1) then
            do 110 i = 0, n-1
               y(i) = x(i)
 110        continue
         else
            do 111 i = 0, n-1
               y(i) = 0.0d0
 111           continue
         endif
c     
c.....do a series of SPAXPYs (sparse daxpys)
         jend = 0
         do 120 j = 0, n-1
            jbgn = jend
            jend = jbgn + profile(j)
            do 130 i = jbgn, jend-1
               y(indx(i)) = y(indx(i)) + val(i)*x(j)
 130        continue
 120     continue
      endif
      return
      end
