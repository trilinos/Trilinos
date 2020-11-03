C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE ARELAX (MXND, XN, YN, LXK, KXL, NXL, LLL, ARFACT)
C***********************************************************************

C  SUBROUTINE ARELAX = CALCULATES UNDER - RELAXATION FACTOR FOR AREA PULL
C                      AND LAPLACIAN SMOOTHER

C***********************************************************************

C  NOTE:
C     THE AREA PULL AND LAPLACIAN SMOOTHER WILL OVER - CORRECT
C     AND BECOME UNSTABLE WHEN TYPICAL MESH ELEMENTS ARE MUCH
C     LONGER THAN THEY ARE WIDE, SAY BY A FACTOR OF SIX OR MORE.
C     THIS ROUTINE COMPUTES AN APPROPRIATE UNDER - RELAXATION
C     FACTOR TO BE USED TO HELP CORRECT THIS PROBLEM.  ON REGIONS
C     WHICH HAVE GENERALLY NEAR SQUARE ELEMENTS WITH A SMALL
C     PERCENTAGE OF VERY LONG THIN ELEMENTS THIS FACTOR WILL
C     PROBABLY NOT ADEQUATELY HANDLE THE DIFFICULTY.  IN SUCH
C     SITUATIONS AN ALTERNATE SMOOTHER  (SUCH AS THE CENTROID - AREA -
C     PULL) SHOULD BE USED.
C     THE FACTOR RETURNED BY THIS ROUTINE MAY BE LARGER THAN ONE,
C     WHICH MEANS THAT OVER - RELAXATION IS APPROPRIATE.

C***********************************************************************

      DIMENSION NODES (4), LXK (4, MXND), KXL (2, 3 * MXND)
      DIMENSION NXL (2, 3 * MXND)
      DIMENSION XN (MXND), YN (MXND)

      LOGICAL CCW

      ARFACT = 1.0
      RATSUM = 0.
      NUM = 0

      DO 100 MYL = 1, LLL

C  SKIP BOUNDARY LINES

         IF (KXL(1, myL) .gt. 0 .and. KXL (2, myL) .GT. 0) THEN
            CCW = .TRUE.
            CALL GNXKA (MXND, XN, YN, KXL (1, MYL), NODES, AREA1, LXK,
     &         NXL, CCW)
            CALL GNXKA (MXND, XN, YN, KXL (2, MYL), NODES, AREA2, LXK,
     &         NXL, CCW)
            N1 = NXL (1, MYL)
            N2 = NXL (2, MYL)
            DXDY = (XN (N2) - XN (N1)) **2 + (YN (N2) - YN (N1)) **2
            IF (AREA1 + AREA2 .GT. 0) THEN
               RATIO = 2.0 * DXDY /  (AREA1 + AREA2)
               IF (RATIO .GE. 0.99) THEN
                  NUM = NUM + 1
                  RATSUM = RATSUM + RATIO
               ENDIF
            ENDIF
         ENDIF
  100 CONTINUE

      IF (NUM .LE. 0) RETURN
      ASPECT = RATSUM / DBLE(NUM)
      ARFACT = AMIN1 (2.0 / ASPECT, 1.5)

      RETURN

      END
