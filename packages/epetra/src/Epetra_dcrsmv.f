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
*\Name: EPETRA_DCRSMV
*
*\Description:
*	EPETRA_DCRSMV performs one of the matrix-vector operations
*
*    y = A*x or y = A'*x.
*
*    where x, y are vectors, A is an
*    m-by-n matrix, and A' is the transpose of A.  The matrix
*    A is stored in a special one dimensional row-pointer
*    format where the zeroes in each row are discarded and
*    the non-zeroes are stored contiguously.  A corresponding
*    integer array, pntr, holds pointers indicating the start of
*    each row in A.  Each element of the array indx contains
*    the row index of the corresponding element of A.
*
*\Usage:
*     call DCRSMV( itrans, udiag, m, n, val, indx, pntr, x, y )
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
*    pntr   Integer array (input)
*           On entry, pntr(j) contains the offset into val and indx
*           for entries in the jth row pntr must have length > = n+1.
*           pntr is unchanged on exit.
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

      subroutine epetra_dcrsmv( itrans, m, n, 
     &                          val, indx, pntr, x, y)
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
      integer itrans, m, n, indx(0:*), pntr(0:*)
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
         jend = pntr(0)
         do 10 j = 0, m-1
            jbgn = jend
            jend = pntr(j+1)
            sum = 0.0
            do 20 i = jbgn, jend-1
               sum = sum + val(i) * x(indx(i))
 20         continue
            y(j) = sum
 10      continue
      else 
         do 111 i = 0, n-1
            y(i) = 0.0d0
 111     continue
c     
c.....do a series of SPAXPYs (sparse daxpys)
         jend = pntr(0)
         do 120 j = 0, m-1
            jbgn = jend
            jend = pntr(j+1)
            do 130 i = jbgn, jend-1
               y(indx(i)) = y(indx(i)) + val(i)*x(j)
 130        continue
 120     continue
      endif
      return
      end
