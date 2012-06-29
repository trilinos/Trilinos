C $Id: closeg.f,v 1.1 1990/11/30 11:04:57 gdsjaar Exp $
C $Log: closeg.f,v $
C Revision 1.1  1990/11/30 11:04:57  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]CLOSEG.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE CLOSEG (MSNAP, SNAPDX, NSNAP, X, Y, II, INDEX, XBOT,
     &   XTOP, YBOT, YTOP)
C***********************************************************************
C
C  SUBROUTINE CLOSEG = SUBROUTINE TO RETURN CLOSEST GRID LINE
C
C***********************************************************************
C
C  SUBROUTINE CALLED BY:
C     DIGIT = A SUBROUTINE TO INPUT GEOMETRY
C
C***********************************************************************
C
C  VARIABLES USED:
C     X      = THE X LOCATION IN USER COORDINATES
C     Y      = THE Y LOCATION IN USER COORDINATES
C
C***********************************************************************
C
      DIMENSION SNAPDX(2, MSNAP), NSNAP(2)
C
C  FIND CLOSEST GRID CROSSING IN X OR Y
C
      XHOLD = X
      YHOLD = Y
      CALL SNAPPT (MSNAP, SNAPDX, NSNAP, XHOLD, YHOLD)
      IF (ABS(XHOLD - X) .LT. ABS(YHOLD - Y)) THEN
         INDEX = 1
      ELSE
         INDEX = 2
         XHOLD = YHOLD
      END IF
C
C  FIND INDEX TO GRID LINE
C
      DO 100 I = 1, NSNAP(INDEX)
         IF (SNAPDX(INDEX, I) .GE. XHOLD) THEN
            II = I
            GO TO 110
         END IF
  100 CONTINUE
      II = NSNAP(INDEX)
  110 CONTINUE
C
C  SET GRID LINE LIMITS
C
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
C
      RETURN
C
      END
