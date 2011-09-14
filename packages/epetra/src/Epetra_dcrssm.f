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
*\Name: EPETRA_DCRSSM
*
*\Description:
*	EPETRA_DCRSSM performs one of the triangular solve operations
*
*    y = T^{-1}*x or y = A^{-T}*x.
*
*    where x and y are multivectors, T is an m-by-n triangular.  The matrix
*    T is stored in a special one dimensional row-pointer
*    format where the zeroes in each row are discarded and
*    the non-zeroes are stored contiguously.  A corresponding
*    integer array, pntr, holds pointers indicating the start of
*    each row in T.  Each element of the array indx contains
*    the column index of the corresponding element of T.
*
*\Usage:
*     call EPETRA_DCRSSM( iupper, itrans, udiag, nodiag, m, n, val, indx, pntr, 
*                         x, ldx, y, ldy, xysame)
*
*
*    iupper Integer (input)
*           On entry, iupper specifies whether or not T is upper/lower triangular.
*           If iupper =  0, T is lower triangular.
*           If iupper <> 0, T is upper triangular..
*           If iupper is any other value, then no operation is
*           performed.
*           The iupper argument is unchanged on exit.
*
*    itrans Integer (input)
*           On entry, itrans specifies the operation to be performed.
*           If itrans =  0, y = T^{-1}*x
*           If itrans <> 1, y = A^{-T}*x.
*           If itrans is any other value, then no operation is
*           performed.
*           The itrans argument is unchanged on exit.
*
*    udiag  Integer (input)
*           On entry, udiag specifies whether or not the matrix should
*           be assumed to have a unit diagonal.
*           If udiag <> 0, add x to the result as though unit 
*           diagonal were present.
*           If udiag is any other value, then no unit diagonal is assumed.
*           The udiag argument is unchanged on exit.
*
*   nodiag  Integer (input)
*           On entry, nodiag specifies whether or not the matrix diagonal
*           is stored in this data structure.  If nodiag <> 0, then
*           it is assume that the diagonal is not present (which only
*           makes sense if udiag = 1).
*           The nodiag argument is unchanged on exit.
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
*    pntr Integer array (input)
*           On entry, pntr(j) contains the the offset into val and indx
*           for entries in the jth row pntr must have length > = n+1.
*           pntr is unchanged on exit.
*
*    x      real*8 array (input)
*           Real array of dimension at least n by nrhs.
*           Before entry, the array x
*           must contain the vector operand x.
*           Unchanged on exit.
*
*    ldx    Integer (input)
*           Stride between row elements of x. Unchanged on exit.
*
*    y      real*8 array (output)
*           Real array of dimension at least m by nrhs.
*           On exit it will contain the solution.
*           Note: It is possible (and desirable) for x and y to be
*           the same vector.
*
*    ldy    Integer (input)
*           Stride between row elements of y. Unchanged on exit.
*
*    xysame Integer (input)
*           On entry, if xysame = 1, it is assumed that x and y
*           are the same vector.  Setting this to 1 will save
*           some copying overhead if itrans = 0.  Otherwise
*           it is not used.
*           xysame is unchanged on exit.
*
*\Remarks:
*    1.  Although the example below stores the elements of each
*        column in natural order, this routine makes only one assumption
*        about the order of the non-zeroes within a column.  It expects
*        that for a lower (upper) triangular matrix, if the diagonal is stored,
*        it is stored at the end (beginning) of each row.
*
*\Examples:
*    If the original matrix is
*
*                   | 11   0   0   0   0 |
*                   | 21  22   0   0   0 |
*                   |  0  32  33   0   0 |
*                   |  0   0   0  44   0 |
*                   | 51  52   0  54  55 |
*
*    then the matrix is assumed to be store as
*
*    val = ( 11  21  22  32  33  44  51  52  54  55 )
*
*    with the corresponding pointer arrays
*
*   indx = (  0   0   1   1   2   3   0   1   3   4 ).
*
*                pntr = (0  3  6  8  9  12)
*
*    Thus, indx(j) indicates the (zero-based) column position of the jth element
*    of val and pntr(i) points to the first element of ith row.
*
*\Enddoc

      subroutine epetra_dcrssm( iupper, itrans, udiag, nodiag, m, n, 
     &                          val, indx, pntr, 
     &                          x, ldx, y, ldy, xysame, nrhs)
*
*     ----------------------------
*     Specifications for arguments
*     ----------------------------
      implicit none
      integer iupper, itrans, udiag, nodiag, m, n
      integer ldx, ldy, xysame, nrhs
      integer indx(0:*), pntr(0:*)
      real*8 val(0:*), x(0:*), y(0:*)

*     ----------------------------------
*     Specifications for local variables
*     ----------------------------------
      integer irhs, nrhs1, ix, iy

*     Strip mine nrhs by fives

*      ix = 0
*      iy = 0
*      do 5 irhs = 1, nrhs
*      call epetra_dcrssv(iupper, itrans, udiag, nodiag, m, n,
*     &    val, indx, pntr, x(ix), y(iy), xysame)
*         ix = ix + ldx
*         iy = iy + ldy
* 5    continue
*      return

      nrhs1 = mod(nrhs, 5)
      if (nrhs1.eq.0) nrhs1 = 5
      ix = 0
      iy = 0
      do 10 irhs = 1, nrhs, 5
         if (itrans.eq.0) then
            call epetra_scrssm5( iupper, udiag, nodiag, m, n, 
     &                          val, indx, pntr,
     &                          x(ix), ldx, y(iy), ldy, xysame, nrhs1)
         else
            call epetra_sccssm5( iupper, udiag, nodiag, m, n,
     &                          val, indx, pntr,
     &                          x(ix), ldx, y(iy), ldy, xysame, nrhs1)
         endif
         ix = ix + nrhs1*ldx
         iy = iy + nrhs1*ldy
         nrhs1 = 5
 10   continue

      return
      end

      subroutine epetra_scrssm5( iupper, udiag, nodiag, m, n,
     &                          val, indx, pntr, 
     &                          x, ldx, y, ldy, xysame, nrhs)
*
*
*     ----------------------------
*     Specifications for arguments
*     ----------------------------
      implicit none
      integer iupper, udiag, nodiag, m, n
      integer ldx, ldy, xysame, nrhs
      integer indx(0:*), pntr(0:*)
      real*8 val(0:*), x(0:*), y(0:*)

*     ----------------------------------
*     Specifications for local variables
*     ----------------------------------
      integer i,j, ibgn, iend, ibgnoff, iendoff
      integer jstrt, jstop, jinc
      integer incx, incy
      real*8 yj1, yj2, yj3, yj4, yj5
      real*8 vali

*
*     --------------------------
*     First executable statement
*     --------------------------
*
      if (iupper.ne.0) then
         jstrt = m - 1
         jstop = 0
         jinc = -1
         if (nodiag.ne.0) then
            ibgnoff = 0
            iendoff = 0
         else
            ibgnoff = 1
            iendoff = 0
         endif
      else
         jstrt = 0
         jstop = m - 1
         jinc = 1
         if (nodiag.ne.0) then
            ibgnoff = 0
            iendoff = 0
         else
            ibgnoff = 0
            iendoff = 1
         endif
      endif
      if (nrhs.eq.1) then
         do 10 j = jstrt, jstop, jinc
            ibgn = pntr(j) + ibgnoff
            iend = pntr(j+1) - iendoff - 1
            yj1 = x(j)
            do 20 i = ibgn, iend
               yj1 = yj1 - val(i) * y(indx(i))
 20         continue
            if (udiag.ne.0) then
              y(j) = yj1
            else
              if (iupper.ne.0) then
                 y(j) = yj1/val(ibgn-1)
              else
                 y(j) = yj1/val(iend+1)
               endif
            endif
 10      continue

      else if (nrhs.eq.2) then
         do 210 j = jstrt, jstop, jinc
            ibgn = pntr(j) + ibgnoff
            iend = pntr(j+1) - iendoff - 1
            incx = j
            yj1 = x(incx)
            incx = incx + ldx
            yj2 = x(incx)
            do 220 i = ibgn, iend
               vali = val(i)
               incy = indx(i)
               yj1 = yj1 - vali * y(incy)
               incy = incy + ldy
               yj2 = yj2 - vali * y(incy)
 220        continue
            if (udiag.ne.0) then
               incy = j
               y(incy) = yj1
               incy = incy + ldy
               y(incy) = yj2
            else
               if (iupper.ne.0) then
                  vali = 1.0/val(ibgn-1)
               else
                  vali = 1.0/val(iend+1)
               endif
               incy = j
               y(incy) = yj1*vali
               incy = incy + ldy
               y(incy) = yj2*vali
            endif
 210     continue

      else if (nrhs.eq.3) then
         do 310 j = jstrt, jstop, jinc
            ibgn = pntr(j) + ibgnoff
            iend = pntr(j+1) - iendoff - 1
            incx = j
            yj1 = x(incx)
            incx = incx + ldx
            yj2 = x(incx)
            incx = incx + ldx
            yj3 = x(incx)
            do 320 i = ibgn, iend
               vali = val(i)
               incy = indx(i)
               yj1 = yj1 - vali * y(incy)
               incy = incy + ldy
               yj2 = yj2 - vali * y(incy)
               incy = incy + ldy
               yj3 = yj3 - vali * y(incy)
 320        continue
            if (udiag.ne.0) then
               incy = j
               y(incy) = yj1
               incy = incy + ldy
               y(incy) = yj2
               incy = incy + ldy
               y(incy) = yj3
            else
               if (iupper.ne.0) then
                  vali = 1.0/val(ibgn-1)
               else
                  vali = 1.0/val(iend+1)
               endif
               incy = j
               y(incy) = yj1*vali
               incy = incy + ldy
               y(incy) = yj2*vali
               incy = incy + ldy
               y(incy) = yj3*vali
            endif
 310     continue
      else if (nrhs.eq.4) then
         do 410 j = jstrt, jstop, jinc
            ibgn = pntr(j) + ibgnoff
            iend = pntr(j+1) - iendoff - 1
            incx = j
            yj1 = x(incx)
            incx = incx + ldx
            yj2 = x(incx)
            incx = incx + ldx
            yj3 = x(incx)
            incx = incx + ldx
            yj4 = x(incx)
            do 420 i = ibgn, iend
               vali = val(i)
               incy = indx(i)
               yj1 = yj1 - vali * y(incy)
               incy = incy + ldy
               yj2 = yj2 - vali * y(incy)
               incy = incy + ldy
               yj3 = yj3 - vali * y(incy)
               incy = incy + ldy
               yj4 = yj4 - vali * y(incy)
 420        continue
            if (udiag.ne.0) then
               incy = j
               y(incy) = yj1
               incy = incy + ldy
               y(incy) = yj2
               incy = incy + ldy
               y(incy) = yj3
               incy = incy + ldy
               y(incy) = yj4
            else
               if (iupper.ne.0) then
                  vali = 1.0/val(ibgn-1)
               else
                  vali = 1.0/val(iend+1)
               endif
               incy = j
               y(incy) = yj1*vali
               incy = incy + ldy
               y(incy) = yj2*vali
               incy = incy + ldy
               y(incy) = yj3*vali
               incy = incy + ldy
               y(incy) = yj4*vali
            endif
 410     continue

      else if (nrhs.eq.5) then
         do 510 j = jstrt, jstop, jinc
            ibgn = pntr(j) + ibgnoff
            iend = pntr(j+1) - iendoff - 1
            incx = j
            yj1 = x(incx)
            incx = incx + ldx
            yj2 = x(incx)
            incx = incx + ldx
            yj3 = x(incx)
            incx = incx + ldx
            yj4 = x(incx)
            incx = incx + ldx
            yj5 = x(incx)
            do 520 i = ibgn, iend
               vali = val(i)
               incy = indx(i)
               yj1 = yj1 - vali * y(incy)
               incy = incy + ldy
               yj2 = yj2 - vali * y(incy)
               incy = incy + ldy
               yj3 = yj3 - vali * y(incy)
               incy = incy + ldy
               yj4 = yj4 - vali * y(incy)
               incy = incy + ldy
               yj5 = yj5 - vali * y(incy)
 520        continue
            if (udiag.ne.0) then
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
            else
               if (iupper.ne.0) then
                  vali = 1.0/val(ibgn-1)
               else
                  vali = 1.0/val(iend+1)
               endif
               incy = j
               y(incy) = yj1*vali
               incy = incy + ldy
               y(incy) = yj2*vali
               incy = incy + ldy
               y(incy) = yj3*vali
               incy = incy + ldy
               y(incy) = yj4*vali
               incy = incy + ldy
               y(incy) = yj5*vali
            endif
 510     continue
      endif
      return
      end

      subroutine epetra_sccssm5( iupper, udiag, nodiag, m, n,
     &                          val, indx, pntr, 
     &                          x, ldx, y, ldy, xysame, nrhs)
*
*
*     ----------------------------
*     Specifications for arguments
*     ----------------------------
      implicit none
      integer iupper, udiag, nodiag, m, n
      integer ldx, ldy, xysame, nrhs
      integer indx(0:*), pntr(0:*)
      real*8 val(0:*), x(0:*), y(0:*)

*     ----------------------------------
*     Specifications for local variables
*     ----------------------------------
      integer i, j, k, ibgn, iend, incy
      integer ibgnoff, iendoff
      integer jstrt, jstop, jinc
      integer iy, ix
      real*8 vali, yj1, yj2, yj3, yj4, yj5
      real*8 ytmp

*
*     --------------------------
*     First executable statement
*     --------------------------
*
c.....initialize soln
      if (xysame.eq.0) then
*
         if (ldy.eq.n.and.ldx.eq.n) then
            do 11 k=0, n*nrhs-1
               y(k) = x(k)
 11         continue
         else
            do 15 k=1, nrhs
               iy = ldy*(k-1)
               ix = ldx*(k-1)
               do 20 i = 0, min(m - 1, n - 1)
                  y(i+iy) = x(i+ix)
 20            continue
 15         continue
         endif
      endif

      if (iupper.ne.0) then
         jstrt = 0
         jstop = m - 1
         jinc = 1
         if (nodiag.ne.0) then
            ibgnoff = 0
            iendoff = 0
         else
            ibgnoff = 1
            iendoff = 0
         endif
      else
         jstrt = m - 1
         jstop = 0
         jinc = -1
         if (nodiag.ne.0) then
            ibgnoff = 0
            iendoff = 0
         else
            ibgnoff = 0
            iendoff = 1
         endif
      endif

      if (nrhs.eq.1) then
         do 120 j = jstrt, jstop, jinc
            ibgn = pntr(j) + ibgnoff
            iend = pntr(j+1) - iendoff -1
            if (udiag.eq.0) then
               if (iupper.ne.0) then
                  y(j) = y(j)/val(ibgn-1)
               else
                  y(j) = y(j)/val(iend+1)
               endif
            endif
            ytmp = y(j)
            do 130 i = ibgn, iend
               y(indx(i)) = y(indx(i)) - val(i)*ytmp
 130        continue
 120     continue
      else if (nrhs.eq.2) then
         do 220 j = jstrt, jstop, jinc
            ibgn = pntr(j) + ibgnoff
            iend = pntr(j+1) - iendoff -1
            if (udiag.eq.0) then
               if (iupper.ne.0) then
                  vali = 1.0/val(ibgn-1)
               else
                  vali = 1.0/val(iend+1)
               endif
               incy = j
               y(incy) = y(incy)*vali
               incy = incy + ldy
               y(incy) = y(incy)*vali
            endif
            incy = j
            yj1 = y(incy)
            incy = incy + ldy
            yj2 = y(incy)
            do 230 i = ibgn, iend
               incy = indx(i)
               vali = val(i)
               y(incy) = y(incy) - vali * yj1
               incy = incy + ldy
               y(incy) = y(incy) - vali * yj2
 230        continue
 220     continue
      else if (nrhs.eq.3) then
         do 320 j = jstrt, jstop, jinc
            ibgn = pntr(j) + ibgnoff
            iend = pntr(j+1) - iendoff - 1
            if (udiag.eq.0) then
               if (iupper.ne.0) then
                  vali = 1.0/val(ibgn-1)
               else
                  vali = 1.0/val(iend+1)
               endif
               incy = j
               y(incy) = y(incy)*vali
               incy = incy + ldy
               y(incy) = y(incy)*vali
               incy = incy + ldy
               y(incy) = y(incy)*vali
            endif
            incy = j
            yj1 = y(incy)
            incy = incy + ldy
            yj2 = y(incy)
            incy = incy + ldy
            yj3 = y(incy)
            do 330 i = ibgn, iend
               incy = indx(i)
               vali = val(i)
               y(incy) = y(incy) - vali * yj1
               incy = incy + ldy
               y(incy) = y(incy) - vali * yj2
               incy = incy + ldy
               y(incy) = y(incy) - vali * yj3
 330        continue
 320     continue
      else if (nrhs.eq.4) then
         do 420 j = jstrt, jstop, jinc
            ibgn = pntr(j) + ibgnoff
            iend = pntr(j+1) - iendoff -1
            if (udiag.eq.0) then
               if (iupper.ne.0) then
                  vali = 1.0/val(ibgn-1)
               else
                  vali = 1.0/val(iend+1)
               endif
               incy = j
               y(incy) = y(incy)*vali
               incy = incy + ldy
               y(incy) = y(incy)*vali
               incy = incy + ldy
               y(incy) = y(incy)*vali
               incy = incy + ldy
               y(incy) = y(incy)*vali
            endif
            incy = j
            yj1 = y(incy)
            incy = incy + ldy
            yj2 = y(incy)
            incy = incy + ldy
            yj3 = y(incy)
            incy = incy + ldy
            yj4 = y(incy)
            do 430 i = ibgn, iend
               incy = indx(i)
               vali = val(i)
               y(incy) = y(incy) - vali * yj1
               incy = incy + ldy
               y(incy) = y(incy) - vali * yj2
               incy = incy + ldy
               y(incy) = y(incy) - vali * yj3
               incy = incy + ldy
               y(incy) = y(incy) - vali * yj4
 430        continue
 420     continue
      else if (nrhs.eq.5) then
         do 520 j = jstrt, jstop, jinc
            ibgn = pntr(j) + ibgnoff
            iend = pntr(j+1) - iendoff -1
            if (udiag.eq.0) then
               if (iupper.ne.0) then
                  vali = 1.0/val(ibgn-1)
               else
                  vali = 1.0/val(iend+1)
               endif
               incy = j
               y(incy) = y(incy)*vali
               incy = incy + ldy
               y(incy) = y(incy)*vali
               incy = incy + ldy
               y(incy) = y(incy)*vali
               incy = incy + ldy
               y(incy) = y(incy)*vali
               incy = incy + ldy
               y(incy) = y(incy)*vali
            endif
            incy = j
            yj1 = y(incy)
            incy = incy + ldy
            yj2 = y(incy)
            incy = incy + ldy
            yj3 = y(incy)
            incy = incy + ldy
            yj4 = y(incy)
            incy = incy + ldy
            yj5 = y(incy)
            do 530 i = ibgn, iend
               incy = indx(i)
               vali = val(i)
               y(incy) = y(incy) - vali * yj1
               incy = incy + ldy
               y(incy) = y(incy) - vali * yj2
               incy = incy + ldy
               y(incy) = y(incy) - vali * yj3
               incy = incy + ldy
               y(incy) = y(incy) - vali * yj4
               incy = incy + ldy
               y(incy) = y(incy) - vali * yj5
 530        continue
 520     continue
      endif
      return
      end
