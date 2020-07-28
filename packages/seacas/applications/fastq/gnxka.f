C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE GNXKA (MXND, XN, YN, K, NODES, AREA, LXK, NXL, CCW)
C***********************************************************************
C  SUBROUTINE GNXKA = GENERATES A LIST OF THE FOUR NODES ASSOCIATED WITH
C                     ELEMENT K

C***********************************************************************
C  VARIABLES USED:
C     CCW    = .TRUE. IF LIST IS TO BE IN CCW ORDER AND AREA DEFINED
C    (Changed to always put in order and calculate area)

C***********************************************************************

      REAL XN (MXND), YN (MXND)
      INTEGER NODES(4)
      INTEGER LXK(4, MXND), NXL(2, 3 * MXND)

      LOGICAL CCW

      AREA = 0.0
      DO 10 I = 1, 4
        NODES (I) = 0
 10   CONTINUE

C... Let line 1 be the base line
      L = LXK(1, K)

C... Null List
      IF (L .LE. 0) THEN
        RETURN
      ENDIF

      NODES(1) = NXL(1, L)
      NODES(2) = NXL(2, L)

C... Find other ends of the two sides

      DO 110 I = 2, 4
         L  = LXK(I, K)
         M1 = NXL(1, L)
         M2 = NXL(2, L)

         IF (M1 .EQ. NODES(1)) THEN
           NODES(4) = M2
         ELSE IF (M2 .EQ. NODES(1)) THEN
           NODES(4) = M1
         END IF

         IF (M1 .EQ. NODES(2)) THEN
           NODES(3) = M2
         ELSE IF (M2 .EQ. NODES(2)) THEN
           NODES(3) = M1
         END IF

  110 CONTINUE

C... Compute signed area
      AREA = 0.5 *
     *  ((XN(NODES(3)) - XN(NODES(1))) *
     *   (YN(NODES(4)) - YN(NODES(2))) -
     &   (YN(NODES(3)) - YN(NODES(1))) *
     *   (XN(NODES(4)) - XN(NODES(2))))

      IF (AREA .LT. 0.) THEN
C ... Clockwise case  -  reverse the order
        NTMP = NODES(2)
        NODES(2) = NODES(4)
        NODES(4) = NTMP
        AREA = -AREA
      ENDIF

      RETURN
      END
