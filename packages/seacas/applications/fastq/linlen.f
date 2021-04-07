C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE LINLEN (MP, COOR, LINKP, KNUM, LNUM, KT, I3, J1, J2,
     &   J3, DIST, ERR)
C***********************************************************************

C  SUBROUTINE LINLEN = CALCULATES THE LENGTH OF A GIVEN LINE

C***********************************************************************

C  VARIABLES USED:
C     NID   = AN ARRAY OF UNIQUE NODE IDENTIFIERS.
C     REAL  = .TRUE. FOR AN ACTUAL GENERATION
C           = .FALSE. FOR A TRIAL GENERATION
C     ERR   = .TRUE. IF AN ERROR WAS ENCOUNTERED
C     J1    = POINTER FOR THE FIRST POINT
C     J2    = POINTER FOR THE SECOND POINT
C     J3    = POINTER FOR THE THIRD POINT
C     MAXNP = MAXIMUM NUMBER OF NODES ON THE PERIMETER
C             NOTE: MAXNP MUST BE ADJUSTED FOR THE CURRENT
C                   LOCATION IN X, Y, & NID
C     KT    = THE LINE TYPE:
C           = 1 FOR STRAIGHT LINES
C           = 2 FOR CORNER LINES
C           = 3 FOR ARC WITH CENTER GIVEN
C           = 4 FOR ARC WITH THIRD POINT ON THE ARC
C           = 5 FOR PARABOLA
C           = 6 FOR ARC WITH RADIUS GIVEN

C***********************************************************************

      DIMENSION COOR (2, MP), LINKP (2, MP)

      LOGICAL ERR

      PI = ATAN2(0.0, -1.0)

      DIST = 0.
      ERR = .TRUE.

C  STRAIGHT LINE GENERATION

      IF (KT.EQ.1) THEN
         YDIFF = COOR (2, J2) -COOR (2, J1)
         XDIFF = COOR (1, J2) -COOR (1, J1)
         DIST = SQRT (YDIFF ** 2 + XDIFF ** 2)
         IF (DIST.EQ.0.) THEN
            WRITE (*, 10000) KNUM
            RETURN
         ENDIF

C  CORNER GENERATION

      ELSEIF (KT.EQ.2) THEN
         XDA = COOR (1, J3) -COOR (1, J1)
         YDA = COOR (2, J3) -COOR (2, J1)
         XDB = COOR (1, J2) -COOR (1, J3)
         YDB = COOR (2, J2) -COOR (2, J3)
         DA = SQRT (XDA ** 2 + YDA ** 2)
         DB = SQRT (XDB ** 2 + YDB ** 2)
         IF ((DA.EQ.0.) .OR. (DB.EQ.0.) )THEN
            WRITE (*, 10000) KNUM
            RETURN
         ENDIF
         DIST = DA+DB

C  CIRCULAR ARC

      ELSEIF ((KT.EQ.3) .OR. (KT.EQ.4) .OR. (KT.EQ.6) )THEN
         XSTART = COOR (1, J1)
         YSTART = COOR (2, J1)
         CALL ARCPAR (MP, KT, KNUM, COOR, LINKP, J1, J2, J3, I3,
     &      XCEN, YCEN, THETA1, THETA2, TANG, R1, R2, ERR, ICCW, ICW,
     &      XK, XA)

C  GENERATE THE CIRCLE

         ANG = THETA1
         DEL = TANG/30
         DO 100 I = 2, 29
            ANG = ANG+DEL
            RADIUS = XA * EXP (XK * ANG)
            XEND = XCEN+COS (ANG) * RADIUS
            YEND = YCEN+SIN (ANG) * RADIUS
            DIST = DIST+SQRT ((XEND-XSTART) ** 2 + (YEND-YSTART) ** 2)
            XSTART = XEND
            YSTART = YEND
  100    CONTINUE
         XEND = COOR (1, J2)
         YEND = COOR (2, J2)
         DIST = DIST+SQRT ((XEND-XSTART) ** 2 + (YEND-YSTART) ** 2)

C  ELIPSE

      ELSEIF (KT .EQ. 7) THEN
         XSTART = COOR (1, J1)
         YSTART = COOR (2, J1)
         CALL ELPSPR (MP, KT, KNUM, COOR, LINKP, J1, J2, J3,
     &      I3, XCEN, YCEN, THETA1, THETA2, TANG, IDUM1, IDUM2,
     &      AVALUE, BVALUE, ERR)

C  GENERATE THE ELIPSE

         ANG = THETA1
         DEL = TANG/30
         DO 110 I = 2, 29
            ANG = ANG+DEL
            RADIUS = SQRT ( (AVALUE **2 * BVALUE **2) /
     &         ( (BVALUE **2 * COS (ANG) **2) +
     &         (AVALUE **2 * SIN (ANG) **2) ) )
            XEND = XCEN+COS (ANG) * RADIUS
            YEND = YCEN+SIN (ANG) * RADIUS
            DIST = DIST+SQRT ((XEND-XSTART) ** 2 + (YEND-YSTART) ** 2)
            XSTART = XEND
            YSTART = YEND
  110    CONTINUE
         XEND = COOR (1, J2)
         YEND = COOR (2, J2)
         DIST = DIST+SQRT ((XEND-XSTART) ** 2 + (YEND-YSTART) ** 2)

C     PARABOLA

      ELSEIF (KT.EQ.5) THEN

C  CHECK LEGITIMACY OF DATA

         XMID =  (COOR (1, J1) +COOR (1, J2) ) * 0.5
         YMID =  (COOR (2, J1) +COOR (2, J2) ) * 0.5
         DOT =  (COOR (1, J2) -COOR (1, J1) ) * (COOR (1, J3) -XMID)
     &      + (COOR (2, J2) -COOR (2, J1) ) * (COOR (2, J3) -YMID)
         PERP = SQRT ((COOR (1, J2) -COOR (1, J1) ) ** 2 +
     &      (COOR (2, J2) - COOR (2, J1) ) ** 2) *
     &      SQRT ((COOR (1, J3) -XMID) ** 2 + (COOR (2, J3)
     &      -YMID) ** 2)
         IF (DOT.GE.0.05 * PERP) THEN
            WRITE (*, 10030) KNUM
            RETURN
         ENDIF

C  GETARC LENGTH

         HALFW = SQRT ((COOR (1, J2) -COOR (1, J1) ) ** 2 +
     &      (COOR (2, J2) - COOR (2, J1) ) ** 2 )  * 0.5
         IF (HALFW.EQ.0.) THEN
            WRITE (*, 10000) KNUM
            RETURN
         ENDIF
         HEIGHT = SQRT ((XMID-COOR (1, J3) ) ** 2
     &      + (YMID-COOR (2, J3) ) ** 2)
         COEF = HEIGHT/HALFW ** 2
         TCOEF = 2.0 * COEF

C  PARC IS A STATEMENT FUNCTION

         PLEFT = PARC (-TCOEF * HALFW, TCOEF)
         ARCTOT = 2.0 * PARC (TCOEF * HALFW, TCOEF)
         ARCDEL = ARCTOT/30
         ARCNXT = ARCDEL
         ARCNOW = 0.0
         THETA = ATAN2 (COOR (2, J2) -COOR (2, J1) , COOR (1, J2)
     &      - COOR (1, J1) )

C  CORRECT FOR ORIENTATION

         CROSS =  (COOR (1, J3) -XMID) *  (COOR (2, J2) -COOR (2, J1) )-
     &      (COOR (2, J3) -YMID) *  (COOR (1, J2) -COOR (1, J1) )
         IF (CROSS.LT.0.0) THETA = THETA+PI
         SINT = SIN (THETA)
         COST = COS (THETA)

C  FIND POINTS APPROXIMATELY BY INTEGRATION

         XL = -HALFW
         FL = SQRT (1.0+ (TCOEF * XL) ** 2)
         KOUNT = 1
         DELX = 2.0 * HALFW/200.0
         XSTART = COOR (1, J1)
         YSTART = COOR (2, J1)
         DO 120 I = 1, 100
            FM = SQRT (1.0+ (TCOEF * (XL+DELX) ) ** 2)
            XR = - HALFW + DBLE(I) * 2.0 * DELX
            FR = SQRT (1.0+ (TCOEF * XR) ** 2)
            ARCOLD = ARCNOW
            ARCNOW = ARCNOW+DELX * (FL+4.0 * FM+FR) / 3.0
            IF (ARCNOW.GE.ARCNXT) THEN

C  COMPUTE POSITION IN LOCAL COORDINATE SYSTEM

               FRAC =  (ARCNXT-ARCOLD) / (ARCNOW-ARCOLD)
               XK = XL+FRAC * 2.0 * DELX
               YK = COEF * XK ** 2

C  CORRECT FOR ORIENTATION PROBLEM

               IF (CROSS.LT.0.0) XK = -XK

C  ROTATE IN LINE WITH GLOBAL COORDINATE SYSTEM

               ROTX = XK * COST - YK * SINT
               ROTY = YK * COST + XK * SINT

C  RESTORE XK

               IF (CROSS.LT.0.0) XK = -XK

C  TRANSLATE

               XEND = ROTX+COOR (1, J3)
               YEND = ROTY+COOR (2, J3)
               DIST = DIST+SQRT ((XEND-XSTART) ** 2 + (YEND-YSTART) **2)
               KOUNT = KOUNT+1
               XSTART = XEND
               YSTART = YEND

C  PREPARE FOR NEXT POINT

               IF (KOUNT.GE.29) GOTO 130
               ARCNXT = ARCNXT+ARCDEL

C  RESTART INTEGRATION

               XR = XK
               FR = SQRT (1.0+ (TCOEF * XR) ** 2)

C  CORRECT FOR INTEGRATION ERROR

               ARCNOW = PARC (TCOEF * XR, TCOEF) -PLEFT
            ENDIF
            XL = XR
            FL = FR
  120    CONTINUE
  130    CONTINUE
         XEND = COOR (1, J2)
         YEND = COOR (2, J2)
         DIST = DIST+SQRT ((XEND-XSTART) ** 2+ (YEND-YSTART) ** 2)
      ENDIF

C     NORMAL EXIT

      ERR = .FALSE.
      RETURN

10000 FORMAT (' ZERO LINE LENGTH ENCOUNTERED FOR LINE', I5)
10030 FORMAT (' POINTS GIVEN FOR LINE', I5, ' DO NOT DEFINE A PARABOLA')

      END
