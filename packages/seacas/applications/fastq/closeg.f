C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE CLOSEG (MSNAP, SNAPDX, NSNAP, X, Y, II, INDEX, XBOT,
     &   XTOP, YBOT, YTOP)
C***********************************************************************

C  SUBROUTINE CLOSEG = SUBROUTINE TO RETURN CLOSEST GRID LINE

C***********************************************************************

C  SUBROUTINE CALLED BY:
C     DIGIT = A SUBROUTINE TO INPUT GEOMETRY

C***********************************************************************

C  VARIABLES USED:
C     X      = THE X LOCATION IN USER COORDINATES
C     Y      = THE Y LOCATION IN USER COORDINATES

C***********************************************************************

      DIMENSION SNAPDX(2, MSNAP), NSNAP(2)

C  FIND CLOSEST GRID CROSSING IN X OR Y

      XHOLD = X
      YHOLD = Y
      CALL SNAPPT (MSNAP, SNAPDX, NSNAP, XHOLD, YHOLD)
      IF (ABS(XHOLD - X) .LT. ABS(YHOLD - Y)) THEN
         INDEX = 1
      ELSE
         INDEX = 2
         XHOLD = YHOLD
      END IF

C  FIND INDEX TO GRID LINE

      DO 100 I = 1, NSNAP(INDEX)
         IF (SNAPDX(INDEX, I) .GE. XHOLD) THEN
            II = I
            GO TO 110
         END IF
  100 CONTINUE
      II = NSNAP(INDEX)
  110 CONTINUE

C  SET GRID LINE LIMITS

      IF (INDEX .EQ. 1) THEN
         XBOT = SNAPDX(1, II)
         XTOP = XBOT
         YBOT = SNAPDX(2, 1)
         YTOP = SNAPDX(2, NSNAP(2))
      ELSE
         XBOT = SNAPDX(1, 1)
         XTOP = SNAPDX(1, NSNAP(1))
         YBOT = SNAPDX(2, II)
         YTOP = YBOT
      END IF

      RETURN

      END
