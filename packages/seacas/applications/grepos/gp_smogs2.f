C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SMOGS2 (xscr, yscr, isbnd, x, y, numel,
     $     link, nlink, numnp, nscr)
C=======================================================================

C***********************************************************************

C  SUBROUTINE SMOGS  =  MESH SMOOTHING BY LAPLACE-S USING GAUSS-SEIDEL

C***********************************************************************

      real x(*), y(*), xscr(*), yscr(*)
      integer numel, nlink, link(nlink, *)
      logical isbnd(*)
      integer nscr(*)

      if (nlink .ne. 4) return
      do 20 iel = 1, numel
         n1 = link(1, iel)
         n2 = link(2, iel)
         n3 = link(3, iel)
         n4 = link(4, iel)

         if (.not. isbnd(n1)) then
            xscr(n1) = xscr(n1) + x(n2) + x(n4)
            yscr(n1) = yscr(n1) + y(n2) + y(n4)
            nscr(n1) = nscr(n1) + 2
         end if

         if (.not. isbnd(n2)) then
            xscr(n2) = xscr(n2) + x(n1) + x(n3)
            yscr(n2) = yscr(n2) + y(n1) + y(n3)
            nscr(n2) = nscr(n2) + 2
         end if

         if (.not. isbnd(n3)) then
            xscr(n3) = xscr(n3) + x(n2) + x(n4)
            yscr(n3) = yscr(n3) + y(n2) + y(n4)
            nscr(n3) = nscr(n3) + 2
         end if

         if (.not. isbnd(n4)) then
            xscr(n4) = xscr(n4) + x(n1) + x(n3)
            yscr(n4) = yscr(n4) + y(n1) + y(n3)
            nscr(n4) = nscr(n4) + 2
         end if

 20   continue
      return
      end
