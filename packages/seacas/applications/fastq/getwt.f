C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE GETWT (MP, ML, MXLPS, NIX, ILIST, XLIST, ILINE, LCON,
     &   LTYPE, COOR, LINKP, FRACT, ADDLNK, ERR)
C***********************************************************************

C  SUBROUTINE GETWT = GETS A WEIGHT BASED ON A PERCENTAGE DISTANCE ALONG
C                     THE GIVEN SIDE LINE LIST.

C***********************************************************************

C  SUBROUTINE CALLED BY:
C     ADDWT = ADDS THE WEIGHTING FACTORS TO ANY NODES WITH
C             FLAGS CONTAINING WEIGHTS

C***********************************************************************

C  VARIABLES USED:
C     FRACT = THE FRACTION OF TOTAL DISTANCE ALONG THE X AXIS
C               (TAKES BACK THE WEIGHT OR Y VALUE AT THAT % DISTANCE)

C***********************************************************************

      DIMENSION ILIST (MXLPS), XLIST (MXLPS)
      DIMENSION ILINE (ML), LCON (3, ML), LTYPE (ML)
      DIMENSION COOR (2, MP), LINKP (2, MP)

      LOGICAL ADDLNK, ERR

      ADDLNK = .FALSE.
      ERR = .FALSE.
      PI = ATAN2(0.0, -1.0)
      TWOPI = 2.*PI

C  GET THE X VALUE

      X = FRACT* (XLIST (NIX)-XLIST (1))+XLIST (1)

C  FIND THE LINE THIS BELONGS TO

      DO 100 I = 1, NIX-1
         IF ((X.LE.XLIST (I+1)).AND. (X.GE.XLIST (I))) THEN
            IL = ILIST (I)
            GOTO 110
         ENDIF
  100 CONTINUE
      CALL MESAGE ('PROBLEMS IN GETWT - NO X SPAN FOUND')
      ERR = .TRUE.
      RETURN

  110 CONTINUE

C  NOW GET THE Y VALUE FOR THE X AND THE LINE  (AND TYPE OF LINE) GIVEN

      KT = LTYPE (IL)
      CALL LTSORT (MP, LINKP, LCON (1, IL), IP1, ADDLNK)
      CALL LTSORT (MP, LINKP, LCON (2, IL), IP2, ADDLNK)
      CALL LTSORT (MP, LINKP, IABS (LCON (3, IL)), IP3, ADDLNK)
      IF (LCON (3, IL).LT.0)IP3 = -IP3

C  CHECK FOR EXACT LINE END PLACEMENT

      EPS = ABS (XLIST (NIX)-XLIST (1))*.00001
      IF (ABS (X-COOR (1, IP1)).LT.EPS) THEN
         FRACT = COOR (2, IP1)
         RETURN
      ELSEIF (ABS (X-COOR (1, IP2)).LT.EPS) THEN
         FRACT = COOR (2, IP2)
         RETURN
      ENDIF

C  GET INTERMEDIATE Y VALUE BASED ON THE LINE TYPE

C  FIRST - STRAIGHT LINES

      IF (KT.EQ.1) THEN
         IF (COOR (1, IP1).GT.COOR (1, IP2)) THEN
            IHOLD = IP1
            IP1 = IP2
            IP2 = IHOLD
         ENDIF
         XFRACT =  (X-COOR (1, IP1))/ (COOR (1, IP2)-COOR (1, IP1))
         FRACT =  (XFRACT* (COOR (2, IP2)-COOR (2, IP1)))+COOR (2, IP1)

C  NEXT - CORNER LINES

      ELSEIF (KT.EQ.2) THEN
         IF (COOR (1, IP1).GT.COOR (1, IP2)) THEN
            IHOLD = IP1
            IP1 = IP2
            IP2 = IHOLD
         ENDIF
         IF (COOR (1, IP3).GT.X) THEN
            IP2 = IP3
            XFRACT =  (X-COOR (1, IP1))/ (COOR (1, IP2)-COOR (1, IP1))
            FRACT = (XFRACT * (COOR (2, IP2) - COOR (2, IP1)))
     &         + COOR (2, IP1)
         ELSEIF (COOR (1, IP3).LT.X) THEN
            IP1 = IP3
            XFRACT =  (X-COOR (1, IP1))/ (COOR (1, IP2)-COOR (1, IP1))
            FRACT =  (XFRACT* (COOR (2, IP2) - COOR (2, IP1)))
     &         + COOR (2, IP1)
         ELSE
            FRACT = COOR (2, IP3)
         ENDIF

C  NEXT - ARCS

      ELSEIF ((KT.EQ.3).OR. (KT.EQ.4).OR. (KT.EQ.6)) THEN

C  ARCWITH CENTER GIVEN
C  ARCGOES FROM 1ST POINT TO 2ND IN *COUNTER-CLOCKWISE* DIRECTION.

         IF (KT.EQ.3) THEN
            XCEN = COOR (1, IABS (IP3))
            YCEN = COOR (2, IABS (IP3))

C  CIRCLE WITH THIRD POINT ON ARC.

         ELSEIF (KT.EQ.4) THEN
            THETA1 = ATAN2 (COOR (2, IP3)-COOR (2, IP1), COOR (1, IP3)-
     &         COOR (1, IP1))+PI/2.0
            THETA2 = ATAN2 (COOR (2, IP3)-COOR (2, IP2), COOR (1, IP3)-
     &         COOR (1, IP2))+PI/2.0
            DET = -COS (THETA1)*SIN (THETA2)+COS (THETA2)*SIN (THETA1)
            X1 = 0.5 * (COOR (1, IP1)+COOR (1, IP3))
            Y1 = 0.5 * (COOR (2, IP1)+COOR (2, IP3))
            X2 = 0.5 * (COOR (1, IP2)+COOR (1, IP3))
            Y2 = 0.5 * (COOR (2, IP2)+COOR (2, IP3))
            R = (-SIN (THETA2) * (X2-X1)+COS (THETA2) * (Y2-Y1))/DET
            XCEN = X1 + R * COS (THETA1)
            YCEN = Y1 + R * SIN (THETA1)

C     CIRCLE WITH RADIUS GIVEN

         ELSEIF (KT.EQ.6) THEN
            DX = 0.5 * (COOR (1, IP2)-COOR (1, IP1))
            DY = 0.5 * (COOR (2, IP2)-COOR (2, IP1))
            CHORD = SQRT (DX*DX+DY*DY)
            R = ABS (COOR (1, IABS (IP3)))
            IF (R.LE.CHORD) THEN
               XCEN = 0.5 * (COOR (1, IP1)+COOR (1, IP2))
               YCEN = 0.5 * (COOR (2, IP1)+COOR (2, IP2))
            ELSE
               ARM = SQRT (R * R-CHORD * CHORD)
               IF (IP3.LT.0) THEN
                  XCEN = COOR (1, IP1)+DX+ARM * DY/CHORD
                  YCEN = COOR (2, IP1)+DY-ARM * DX/CHORD
               ELSE
                  XCEN = COOR (1, IP1)+DX-ARM * DY/CHORD
                  YCEN = COOR (2, IP1)+DY+ARM * DX/CHORD
               ENDIF
            ENDIF
         ENDIF
         R1 = SQRT ((COOR (1, IP1)-XCEN) **2 + (COOR (2, IP1)-YCEN) **2)
         R2 = SQRT ((COOR (1, IP2)-XCEN) **2 + (COOR (2, IP2)-YCEN) **2)
         IF ((R1.EQ.0.).OR. (R2.EQ.0.)) THEN
            ERR = .TRUE.
            WRITE (*, 10000)ILINE (IL)
            RETURN
         ENDIF
         THETA1 = ATAN2 (COOR (2, IP1)-YCEN, COOR (1, IP1)-XCEN)
         THETA2 = ATAN2 (COOR (2, IP2)-YCEN, COOR (1, IP2)-XCEN)

C  ARCWITH THE CENTER GIVEN

         IF (KT.EQ.3) THEN
            IF ((IP3.GE.0).AND. (THETA2.LE.THETA1))THETA2 = THETA2+TWOPI
            IF ((IP3.LT.0).AND. (THETA1.LE.THETA2))THETA1 = THETA1+TWOPI
            TANG = THETA2-THETA1

C  CIRCULAR ARC WITH 3RD POINT ON ARC - CLOCKWISE OR COUNTER-CLOCKWISE

         ELSEIF (KT.EQ.4) THEN
            THETA3 = ATAN2 (COOR (2, IP3)-YCEN, COOR (1, IP3)-XCEN)
            IF (THETA2.LE.THETA1)THETA2 = THETA2+TWOPI
            IF (THETA3.LE.THETA1)THETA3 = THETA3+TWOPI
            TANG = THETA2-THETA1
            IF (THETA3.GT.THETA2)TANG = - (TWOPI-TANG)

C  CIRCULAR ARC WITH RADIUS GIVEN - CLOCKWISE OR COUNTER-CLOCKWISE

         ELSEIF (KT.EQ.6) THEN
            IF ((IP3.GE.0).AND. (THETA2.LE.THETA1))THETA2 = THETA2+TWOPI
            IF ((IP3.LT.0).AND. (THETA1.LE.THETA2))THETA1 = THETA1+TWOPI
            TANG = THETA2-THETA1
         ENDIF

C  NOW ITERATE UNTIL THE X VALUE IS WITHIN SOME EPSILON

         AA =  (LOG (R2/R1))/ (THETA2-THETA1)
         BB = R2/EXP (AA * THETA2)
         ANG = THETA1
         EPS = ABS (COOR (1, IP1)-COOR (1, IP2)) * .000001
         DO 140 I = 1, 10
            DEL = TANG * .1
            DO 120 J = 1, 10
               ANG = ANG+DEL
               RADIUS = BB * EXP (AA * ANG)
               XTEST = XCEN+COS (ANG) * RADIUS
               IF (EPS.GE.ABS (XTEST-X)) THEN
                  FRACT = YCEN+SIN (ANG) * RADIUS
                  GOTO 150
               ELSEIF ((COOR (1, IP1) .LT. COOR (1, IP2))
     &            .AND. (XTEST .GT. X)) THEN
                  ANG = ANG-DEL
                  TANG = DEL
                  GOTO 130
               ELSEIF ((COOR (1, IP1) .GT. COOR (1, IP2))
     &            .AND. (XTEST .LT. X)) THEN
                  ANG = ANG-DEL
                  TANG = DEL
                  GOTO 130
               ENDIF
  120       CONTINUE
  130       CONTINUE
  140    CONTINUE
         ERR = .TRUE.
         WRITE (*, 10010)ILINE (IL)
         RETURN
  150    CONTINUE

C  FINALLY PARABOLAS

      ELSEIF (KT.EQ.5) THEN

C  CHECK LEGITIMACY OF DATA

         IF (COOR (1, IP1).GT.COOR (1, IP2)) THEN
            IJK = IP1
            IP1 = IP2
            IP2 = IJK
         ENDIF
         XMID =  (COOR (1, IP1)+COOR (1, IP2)) * 0.5
         YMID =  (COOR (2, IP1)+COOR (2, IP2)) * 0.5
         DOT =  (COOR (1, IP2)-COOR (1, IP1)) * (COOR (1, IP3)-XMID)
     &      + (COOR (2, IP2)-COOR (2, IP1)) * (COOR (2, IP3)-YMID)
         PERP = SQRT ((COOR (1, IP2)-COOR (1, IP1)) **2+ (COOR (2, IP2)-
     &      COOR (2, IP1)) **2) * SQRT ((COOR (1, IP3)-XMID) **2
     &      + (COOR (2, IP3) - YMID) **2)
         IF (DOT.GE.0.05 * PERP) THEN
            WRITE (*, 10020)ILINE (IL)
            ERR = .TRUE.
            RETURN
         ENDIF

C  GET TRANSFORMATION TO PARABOLA COORDINATE SYSTEM  (Y = 4AX **2)

         HALFW = SQRT ((COOR (1, IP2)-COOR (1, IP1)) **2 +
     &      (COOR (2, IP2) - COOR (2, IP1)) **2) *0.5
         HEIGHT = SQRT ((XMID-COOR (1, IP3)) **2 +
     &      (YMID - COOR (2, IP3)) **2)
         IF ((HEIGHT.EQ.0).OR. (HALFW.EQ.0.)) THEN
            WRITE (*, 10030)ILINE (IL)
            ERR = .TRUE.
            RETURN
         ENDIF
         A = HEIGHT/ (HALFW **2)
         XTOP = COOR (1, IP3)
         YTOP = COOR (2, IP3)
         THETA = ATAN2 (YMID-YTOP, XMID-XTOP)
         SINT = SIN (THETA)
         COST = COS (THETA)
         IF (SINT.EQ.0.0) THEN
            WRITE (*, 10040)ILINE (IL)
            ERR = .TRUE.
            RETURN
         ENDIF
         COTT = COST/SINT

C  FIND THE EQUATION OF THE LINE FOR X  =  CONSTANT IN NEW COORDINATES

         X0 = X-XTOP
         B = - (SINT * X0)- (COTT * COST * X0)

C  IF THE LINE HAS A ZERO SLOPE,  THEN FIND THE SIMPLE SOLUTION

         IF (COTT.EQ.0.0) THEN
            YNEW = B
         ELSE
            DIVIS = 1.- (4. * COTT * A * B)
            IF (DIVIS.LT.0.0) THEN
               WRITE (*, 10050)ILINE (IL)
               ERR = .TRUE.
               RETURN
            ENDIF
            XDIVIS = SQRT (DIVIS)
            Y1 =  (1.+XDIVIS)/ (2. * COTT * A)
            Y2 =  (1.-XDIVIS)/ (2. * COTT * A)
            IF ((ABS (Y1).LE.HALFW).AND. (ABS (Y2).GT.HALFW)) THEN
               YNEW = Y1
            ELSEIF ((ABS (Y2).LE.HALFW).AND. (ABS (Y1).GT.HALFW)) THEN
               YNEW = Y2
            ELSE
               WRITE (*, 10060)ILINE (IL)
            ENDIF
         ENDIF

C  TRANSLATE THIS XNEW TO A Y VALUE

         XNEW = A * YNEW * YNEW
         FRACT =  (XNEW * SINT)+ (YNEW * COST)+YTOP
      ENDIF

      RETURN

10000 FORMAT (' POINTS GIVEN FOR LINE', I5, ' DO NOT DEFINE AN ARC')
10010 FORMAT (' NO X ON ARC LINE', I5, ' FOUND IN GETWT')
10020 FORMAT (' POINTS FOR LINE', I5, ' DOES NOT DEFINE A PARABOLA')
10030 FORMAT (' ZERO LINE LENGTH FOR PARABOLA LINE', I5, ' IN GETWT')
10040 FORMAT (' PARABOLA ALIGNMENT PROBLEMS FOR LINE', I5, ' IN GETWT')
10050 FORMAT (' PARABOLA INTERSECTION PROBLEMS FOR LINE', I5,
     &   ' IN GETWT')
10060 FORMAT (' PARABOLA SOLUTION PROBLEMS FOR LINE', I5, ' IN GETWT')

      END
