C @HEADER
C ************************************************************************
C 
C               Epetra: Linear Algebra Services Package 
C                 Copyright 2011 Sandia Corporation
C 
C Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C the U.S. Government retains certain rights in this software.
C
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C
C 1. Redistributions of source code must retain the above copyright
C notice, this list of conditions and the following disclaimer.
C
C 2. Redistributions in binary form must reproduce the above copyright
C notice, this list of conditions and the following disclaimer in the
C documentation and/or other materials provided with the distribution.
C
C 3. Neither the name of the Corporation nor the names of the
C contributors may be used to endorse or promote products derived from
C this software without specific prior written permission.
C
C THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
C EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
C IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
C PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
C CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
C EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
C PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
C PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
C LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
C NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
C SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C
C Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
C 
C ************************************************************************
C @HEADER

*----------------------------------------------------------------------
*
*
*\Documentation
*
*\Name: EPETRA_DCRSMM
*
*\Description:
*	EPETRA_DCRSMM performs one of the matrix-multivector operations
*
*    y = A*x or y = A'*x.
*
*    where x, y are multivectors, A is an
*    m-by-n matrix, and A' is the transpose of A.  The matrix
*    A is stored in a special one dimensional row-pointer
*    format where the zeroes in each row are discarded and
*    the non-zeroes are stored contiguously.  A corresponding
*    integer array, pntr, holds pointers indicating the start of
*    each row in A.  Each element of the array indx contains
*    the row index of the corresponding element of A.
*
*\Usage:
*     call DCRSMM(  itrans, m, n, val, indx, pntr, x, ldx, y, ldy, nrhs )
*
*
*    itrans Integer (input)
*           On entry, itrans specifies the operation to be performed.
*           If itrans =  0, z := A*x + beta*y.
*           If itrans <> 1, z := A'*x + beta*y.
*           If itrans is any other value, then no operation is
*           performed.
*           The itrans argument is unchanged on exit.
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
*    pntr   Integer array (input)
*           On entry, pntr(j) contains the offset into val and indx 
*           for entries in the jth row pntr must have length > = n+1.
*           pntr is unchanged on exit.
*
*    x      real*8 array (input)
*           Real array of dimension at least n by nrhs when itrans = 0
*           and at least m by nrhs otherwise.  Before entry, the array x
*           must contain the vector operand x.
*           The argument x is unchanged on exit.
*
*    ldx    Integer (input)
*           Stride between row elements of x. Unchanged on exit.
*
*    y      real*8 array (output)
*           Real array of dimension at least m by nrhs when itrans = 0 and
*           at least n by nrhs otherwise.  On exit, y contains the result
*           vector.
*
*    ldy    Integer (input)
*           Stride between row elements of y. Unchanged on exit.
*
*    nrhs   Integer (input)
*           Number of vectors in the multivectors x and y. Unchanged on exit.
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
*    val = ( 11  12  14  21  22  23  32  33  45  51  52  54 )
*
*    with the corresponding pointer arrays
*
*   indx = (  0   1   3   0   1   2   1   2   4   0   1   3 ).
*
*                pntr = (0  3  6  8  9  12)
*
*    Thus, indx(j) indicates the column position of the jth element
*    of val and pntr(i) points to the first element of ith row.
*
*\Enddoc

      subroutine epetra_dcrsmm( itrans, m, n, val, indx, pntr, 
     &                          x, ldx, y, ldy, nrhs)

*     ----------------------------
*     Specifications for arguments
*     ----------------------------
      implicit none
      integer itrans, m, n, ldx, ldy
      integer nrhs, indx(0:*), pntr(0:*)
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
            call epetra_scrsmm5(m, n, val, indx, pntr,
     &                          x(ix), ldx, y(iy), ldy, nrhs1)
         else
            call epetra_sccsmm5(m, n, val, indx, pntr,
     &                          x(ix), ldx, y(iy), ldy, nrhs1)
         endif
         ix = ix + nrhs1*ldx
         iy = iy + nrhs1*ldy
         nrhs1 = 5
 10   continue

      return
      end

      subroutine epetra_sccsmm5( m, n, val, indx, pntr, 
     &                           x, ldx, y, ldy, nrhs)
*
*     ----------------------------
*     Specifications for arguments
*     ----------------------------
      implicit none
      integer m, n, nrhs, ldx, ldy
      integer indx(0:*), pntr(0:*)
      real*8 val(0:*), x(0:*), y(0:*)
*
*     ----------------------------------
*     Specifications for local variables
*     ----------------------------------
      integer i, j, k, jbgn, jend, incx, incy
      integer iy
      real*8 vali, xj1, xj2, xj3, xj4, xj5

*
*     --------------------------
*     First executable statement
*     --------------------------
*
c.....initialize soln
*
      if (ldy.eq.n) then
         do 11 k=0, n*nrhs-1
            y(k) = 0.0
 11      continue
      else
         do 15 k=1, nrhs
            iy = ldy*(k-1)
            do 20 i = 0, n-1
               y(i+iy) = 0.0D0
 20         continue
 15      continue
      endif
c
c.....do a series of SPAXPYs (sparse saxpys)
      if (nrhs.eq.1) then
         do 120 j = 0, m-1
            jbgn = pntr(j)
            jend = pntr(j+1)
            xj1 = x(j)
            do 130 i = jbgn, jend-1
               y(indx(i)) = y(indx(i)) + val(i) * xj1
 130        continue
 120     continue
         
      else if (nrhs.eq.2) then
         
      do 220 j = 0, m-1
         jbgn = pntr(j)
         jend = pntr(j+1)
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
         
      do 320 j = 0, m-1
         jbgn = pntr(j)
         jend = pntr(j+1)
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
         
      do 420 j = 0, m-1
         jbgn = pntr(j)
         jend = pntr(j+1)
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
         
      do 520 j = 0, m-1
         jbgn = pntr(j)
         jend = pntr(j+1)
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

      subroutine epetra_scrsmm5( m, n, val, indx, pntr, 
     &                           x, ldx, y, ldy, nrhs)
*
*     ----------------------------
*     Specifications for arguments
*     ----------------------------
      implicit none
      integer m, n, nrhs, ldx, ldy
      integer indx(0:*), pntr(0:*)
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
!$omp parallel default(private) shared (m, pntr, val, x, indx, y)
!$omp do 
         do 110 j = 0, m-1
            jbgn = pntr(j)
            jend = pntr(j+1)
            yj1 = 0.0
            do 120 i = jbgn, jend-1
               yj1 = yj1 + val(i) * x(indx(i))
 120        continue
            y(j) = yj1
 110     continue
!$omp end parallel
         

      else if (nrhs.eq.2) then
!$omp parallel default(private) shared (m, pntr, val, x, indx, y)
!$omp do
         do 210 j = 0, m-1
            jbgn = pntr(j)
            jend = pntr(j+1)
            yj1 = 0.0
            yj2 = 0.0
            do 220 i = jbgn, jend-1
               vali = val(i)
               incx = indx(i)
               yj1 = yj1 + vali * x(incx)
               incx = incx + ldx
               yj2 = yj2 + vali * x(incx)
 220        continue
            incy = j
            y(incy) = yj1
            incy = incy + ldy
            y(incy) = yj2
 210     continue
!$omp end parallel
         
      else if (nrhs.eq.3) then
!$omp parallel default(private) shared (m, pntr, val, x, indx, y)
!$omp do
         do 310 j = 0, m-1
            jbgn = pntr(j)
            jend = pntr(j+1)
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
            incy = j
            y(incy) = yj1
            incy = incy + ldy
            y(incy) = yj2
            incy = incy + ldy
            y(incy) = yj3
 310     continue
!$omp end parallel

      else if (nrhs.eq.4) then
!$omp parallel default(private) shared (m, pntr, val, x, indx, y)
!$omp do
         do 410 j = 0, m-1
            jbgn = pntr(j)
            jend = pntr(j+1)
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
            incy = j
            y(incy) = yj1
            incy = incy + ldy
            y(incy) = yj2
            incy = incy + ldy
            y(incy) = yj3
            incy = incy + ldy
            y(incy) = yj4
 410     continue
!$omp end parallel

      else if (nrhs.eq.5) then
!$omp parallel default(private) shared (m, pntr, val, x, indx, y)
!$omp do
         do 510 j = 0, m-1
            jbgn = pntr(j)
            jend = pntr(j+1)
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
 510     continue
!$omp end parallel


      endif
      return
      end
