C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SMOGS3 (xscr, yscr, zscr, isbnd, x, y, z, numel,
     $     link, nlink, numnp, nscr)
C=======================================================================

C***********************************************************************

C  SUBROUTINE SMOGS  =  MESH SMOOTHING BY LAPLACE-S USING GAUSS-SEIDEL

C***********************************************************************

      real x(*), y(*), z(*), xscr(*), yscr(*), zscr(*)
      integer numel, nlink, link(nlink, *)
      logical isbnd(*)
      integer nscr(*)
      integer map(3,8)

      data map /2,4,5,  3,1,6,  4,2,7,  1,3,8,
     $          1,8,6,  2,5,7,  3,6,8,  4,7,5/

      if (nlink .ne. 8) return
      do 20 iel = 1, numel
         do 10 ino = 1, nlink

            ni = link(ino, iel)
            if (.not. isbnd(ni)) then
               n1 = link( map(1, ino), iel)
               n2 = link( map(2, ino), iel)
               n3 = link( map(3, ino), iel)
               xscr(ni) = xscr(ni) + x(n1) + x(n2) + x(n3)
               yscr(ni) = yscr(ni) + y(n1) + y(n2) + y(n3)
               zscr(ni) = zscr(ni) + z(n1) + z(n2) + z(n3)
               nscr(ni) = nscr(ni) + 3
            end if
 10      continue
 20   continue
      return
      end
