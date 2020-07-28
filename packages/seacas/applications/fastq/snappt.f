C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE SNAPPT (MSNAP, SNAPDX, NSNAP, X, Y)
C***********************************************************************

C  SUBROUTINE SNAPPT = GETS THE X, Y TO THE CLOSEST GRID POINT

C***********************************************************************

C  VARIABLES USED:
C     SNAPDX = SNAP GRID LINE ARRAY
C     NSNAP  = ARRAY OF GRID LINE COUNTERS
C     MSNAP  = THE DIMENSION OF THE GRID LINE ARRAY
C     XMID   = .TRUE. IF THE X VALUE FALLS BETWEEN TWO X GRID LINES
C     YMID   = .TRUE. IF THE Y VALUE FALLS BETWEEN TWO Y GRID LINES

C***********************************************************************

      DIMENSION SNAPDX (2, MSNAP), NSNAP (2)

      LOGICAL XMID,  YMID

      XMID = .FALSE.
      YMID = .FALSE.

C  GET THE BOUNDING X SNAP LINES

      IF (X .LE. SNAPDX (1, 1)) THEN
         XP = SNAPDX (1, 1)
      ELSEIF (X .GE. SNAPDX (1, NSNAP (1))) THEN
         XP = SNAPDX (1, NSNAP (1))
      ELSE
         XMID = .TRUE.
         DO 100 I = 2, NSNAP (1)
            IF (X .LE. SNAPDX (1, I)) THEN
               XP1 = SNAPDX (1, I - 1)
               XP2 = SNAPDX (1, I)
               GOTO 110
            ENDIF
  100    CONTINUE
  110    CONTINUE
      ENDIF

C  GET THE BOUNDING Y SNAP LINES

      IF (Y .LE. SNAPDX (2, 1)) THEN
         YP = SNAPDX (2, 1)
      ELSEIF (Y .GE. SNAPDX (2, NSNAP (2))) THEN
         YP = SNAPDX (2, NSNAP (2))
      ELSE
         YMID = .TRUE.
         DO 120 I = 2, NSNAP (2)
            IF (Y .LE. SNAPDX (2, I)) THEN
               YP1 = SNAPDX (2, I - 1)
               YP2 = SNAPDX (2, I)
               GOTO 130
            ENDIF
  120    CONTINUE
  130    CONTINUE
      ENDIF

C  NOW GET THE APPROPRIATE COMBINATION OF XLOW, XHIGH, XMID, YLOW,  ETC.

C  FIRST THE MOST COMMON CASE OF FITTING BETWEEN X AND Y GRIDS

      IF ( (YMID) .AND. (XMID)) THEN

C  GET THE SHORTEST DISTANCE TO THIS COMBINATION

         DIST1 = SQRT ( ( (XP1 - X) ** 2) +  ( (YP1 - Y) ** 2))
         DIST2 = SQRT ( ( (XP1 - X) ** 2) +  ( (YP2 - Y) ** 2))
         DIST3 = SQRT ( ( (XP2 - X) ** 2) +  ( (YP1 - Y) ** 2))
         DIST4 = SQRT ( ( (XP2 - X) ** 2) +  ( (YP2 - Y) ** 2))

         IF (DIST1 .LE. AMIN1 (DIST2, DIST3, DIST4)) THEN
            X = XP1
            Y = YP1
         ELSEIF (DIST2 .LE. AMIN1 (DIST1, DIST3, DIST4)) THEN
            X = XP1
            Y = YP2
         ELSEIF (DIST3 .LE. AMIN1 (DIST1, DIST2, DIST4)) THEN
            X = XP2
            Y = YP1
         ELSE
            X = XP2
            Y = YP2
         ENDIF

C  NOW THE CORNER CASES OF XLOW,  XHIGH,  YLOW,  AND YHIGH COMBINATIONS

      ELSEIF ( (.NOT.XMID) .AND. (.NOT.YMID)) THEN
         X = XP
         Y = YP

C  NOW THE EDGE CASES OF XLOW OR XHIGH AND YMID

      ELSEIF (.NOT.XMID) THEN
         X = XP
         DIST1 = SQRT ( ( (XP - X) ** 2) +  ( (YP1 - Y) ** 2))
         DIST2 = SQRT ( ( (XP - X) ** 2) +  ( (YP2 - Y) ** 2))
         IF (DIST1 .LT. DIST2) THEN
            Y = YP1
         ELSE
            Y = YP2
         ENDIF

C  NOW THE EDGE CASES OF XMID AND YHIGH OR YLOW

      ELSEIF (.NOT.YMID) THEN
         Y = YP
         DIST1 = SQRT ( ( (XP1 - X) ** 2) +  ( (YP - Y) ** 2))
         DIST2 = SQRT ( ( (XP2 - X) ** 2) +  ( (YP - Y) ** 2))
         IF (DIST1 .LT. DIST2) THEN
            X = XP1
         ELSE
            X = XP2
         ENDIF

C  NOW A CHECK TO MAKE SURE THAT SOMETHING DIDN'T FALL THROUGH

      ELSE
         CALL MESAGE (' **  ERROR  -  IMPOSSIBLE CASE IN SNAPPT  ** ')
         CALL PLTBEL
      ENDIF

      RETURN

      END
