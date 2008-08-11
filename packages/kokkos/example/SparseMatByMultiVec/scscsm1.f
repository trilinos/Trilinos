C     Unfinished!!!!!!!!

      subroutine scscmm1( m, n, valL, indxL, pntrL, valU, indxU, pntrU,
     $                    x, y, nrhs)
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
      integer m, n, nrhs, indxL(*), pntrL(*), indxU(*), pntrU(*)
      real*8 valL(*), valU(*), x(*), y(*)
*
*     ----------------------------------
*     Specifications for local variables
*     ----------------------------------
      integer i,j,jbgn, jend, indxi, rem, chunk, k
      real*8 vali, xj, xj1, xj2, xj3
*
*     --------------------------
*     First executable statement
*     --------------------------
*
c.....initialize soln 
      do 10 i = 1, m*nrhs
         y(i) = x(i)
 10   continue

      rem = mod(nrhs,chunk)
c
c     ------------------------
c     Forward Solve: L * y = x
c     ------------------------
c
c.....do a series of SPAXPYs (sparse saxpys)
      do 20 j = 1, n
        jbgn = pntrL(j)+1
        jend = pntrL(j+1) - 1
        do 30 i = jbgn, jend
           indxi = indxL(i)
           vali = valL(i)
           if (rem.eq.1) then
              y(indxi) = y(indxi) - vali * y(j)
          else if (rem.eq.2) then
              y(indxi) = y(indxi) - vali * y(j)
              y(indxi+1) = y(indxi+1) - vali * y(j+1)
          else if (rem.eq.3) then
              y(indxi) = y(indxi) - vali * y(j)
              y(indxi+1) = y(indxi+1) - vali * y(j+1)
              y(indxi+2) = y(indxi+2) - vali * y(j+2)
           endif
           do 40 k = rem+1, nrhs, chunk
              y(indxi) = y(indxi) - vali * y(j)
              y(indxi+1) = y(indxi+1) - vali * y(j+1)
              y(indxi+2) = y(indxi+2) - vali * y(j+2)
              y(indxi+3) = y(indxi+3) - vali * y(j+3)
 40        continue
 30     continue
 20   continue

c
c     -------------------------
c     Backward Solve: U * x = y
c     -------------------------
      do 50 j = n, 1, -1
        jbgn = pntrU(j)+1
        jend = pntrU(j+1) - 1
        do 30 i = jbgn, jend
           indxi = indxU(i)
           vali = valU(i)
           if (rem.eq.1) then
              y(indxi) = y(indxi) - vali * y(j)
          else if (rem.eq.2) then
              y(indxi) = y(indxi) - vali * y(j)
              y(indxi+1) = y(indxi+1) - vali * y(j+1)
          else if (rem.eq.3) then
              y(indxi) = y(indxi) - vali * y(j)
              y(indxi+1) = y(indxi+1) - vali * y(j+1)
              y(indxi+2) = y(indxi+2) - vali * y(j+2)
           endif
           do 40 k = rem+1, nrhs, chunk
              y(indxi) = y(indxi) - vali * y(j)
              y(indxi+1) = y(indxi+1) - vali * y(j+1)
              y(indxi+2) = y(indxi+2) - vali * y(j+2)
              y(indxi+3) = y(indxi+3) - vali * y(j+3)
 40        continue
 30     continue
 20   continue

      return
      end
