      subroutine scscsv(n, valL, indxL, pntrL, valU, indxU, pntrU, b, x)
*
*     Performs the matrix-vector operation
*
*                               x = (LU)(inverse)b
*
*     where x and b are vectors and L and U are sparse matrices stored
*     in compress column format.
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
      do 30 k = 1, n
         klstrt = pntrL(k) + 1
         klstop = pntrL(k+1) - 1
c
         temp = x(k)
cdir$ ivdep
         do 20 i = klstrt, klstop
            x(indxL(i)) = x(indxL(i)) - temp*valL(i)
   20    continue
   30 continue
c
c     -------------------------
c     Backward Solve: U * x = y
c     -------------------------
      do 50 k = n, 1, -1
         kustrt = pntrU(k)
         kustop = pntrU(k+1) - 1
c
         x(k) = x(k) / valU(kustop)
         temp = x(k)
cdir$ ivdep
         do 40 i = kustrt, kustop-1
            x(indxU(i)) = x(indxU(i)) - temp*valU(i)
   40    continue
   50 continue
      return
      end
