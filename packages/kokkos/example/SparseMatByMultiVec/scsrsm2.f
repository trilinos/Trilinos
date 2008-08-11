      subroutine scsrsv(n, valL, indxL, pntrL, valU, indxU, pntrU, b, x)
*
*     Performs the matrix-vector operation
*
*                               x = (LU)(inverse)b
*
*     where x and b are vectors and L and U are sparse matrices stored
*     in compress row format.
*
*     ----------------------------
*     Specifications for arguments
*     ----------------------------
      implicit double precision (a-h,o-z)
      integer n, indxL(*), pntrL(n+1), indxU(*), pntrU(n+1)
      real*8 valL(*), valU(*), b(n), x(n)
*
*     ----------------------------------
*     Specifications for local variables
*     ----------------------------------
      integer i,k
*
*     --------------------------
*     First executable statement
*     --------------------------
*
c     -----------------
c     Initialize x to b
c     -----------------
      do 10 i = 1, n
         x(i) = b(i)
   10 continue
c
c     ------------------------
c     Forward Solve: L * y = b
c     ------------------------
      do 80 k = 1, n
         klstrt = pntrL(k)
         klstop = pntrL(k+1) - 2
c
         temp = x(k)
cdir$ ivdep
         do 70 i = klstop, klstrt, -1
            temp = temp - valL(i)*x(indxL(i))
   70    continue
         x(k) = temp
  80  continue
c
c     -------------------------
c     Backward Solve: U * x = y
c     -------------------------
      do 100 k = n, 1, -1
         kustrt = pntrU(k)
         kustop = pntrU(k+1) - 1
c
         temp = x(k)
cdir$ ivdep
         do 90 i = kustrt+1, kustop
            temp = temp - valU(i)*x(indxU(i))
   90    continue
         x(k) = temp / valU(kustrt)
  100 continue
      return
      end
