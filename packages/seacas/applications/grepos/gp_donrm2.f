C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DONRM2 (NODES, ISBND, COSIN, NUMEL, LINK, X, Y,
     $     NLINK, NUMNP)
C=======================================================================
      integer NODES(*)
      logical isbnd(*)
      REAL X(*), Y(*), COSIN(2,*)
      INTEGER LINK(NLINK, *)
      parameter (tol = 1.0e-6)

      do 10 inod = 1, numnp
         nodes(inod) = 0
 10   continue

      do 30 ilnk = 1, nlink
         do 20 iel = 1, numel
            nodes(link(ilnk, iel)) = 1
 20      continue
 30   continue

      do 60 iel = 1, numel
         do 50 iseg=1,nlink
            if (iseg .eq. nlink) then
               isegp1 = 1
            else
               isegp1 = iseg + 1
            end if

            XI = x( link(ISEG,iel) )
            YI = y( link(ISEG,iel) )

            XJ = x( link(ISEGp1,iel) )
            YJ = y( link(isegp1,iel) )

            DX = XI - XJ
            DY = YI - YJ
            RMAG = SQRT ( DX**2 + DY**2)

            cosin(1, link(iseg, iel)) = cosin(1, link(iseg, iel)) -
     $           dy / rmag
            cosin(2, link(iseg, iel)) = cosin(2, link(iseg, iel)) +
     $           dx / rmag

            cosin(1, link(isegp1, iel)) = cosin(1, link(isegp1, iel)) -
     $           dy / rmag
            cosin(2, link(isegp1, iel)) = cosin(2, link(isegp1, iel)) +
     $           dx / rmag
 50      continue
 60   continue

      do 70 inod = 1, numnp
         if (nodes(inod) .ne. 0) then
            if (abs(cosin(1, inod)) .gt. tol .or.
     $           abs(cosin(2, inod)) .gt. tol) then
               nodes(inod) = -1
               isbnd(inod) = .TRUE.
            end if
         end if
 70   continue

      RETURN
      END
