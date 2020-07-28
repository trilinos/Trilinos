C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE CLOSEL (MP, ML, N, COOR, ILINE, LTYPE, LCON, LINKP,
     &   LINKL, X, Y, BIFIND, IFIND, ADDCEN, XCHOLD, YCHOLD)
C***********************************************************************

C  SUBROUTINE CLOSEL = FINDS CLOSEST PERPENDICULAR BISECTED LINE

C***********************************************************************

C  SUBROUTINE CALLED BY:
C     BISECT   =

C***********************************************************************

C  SUBROUTINES CALLED:
C     DLPARA = DETERMINES LINE PARAMETERS FROM TWO POINTS

C***********************************************************************

      DIMENSION COOR(2, MP), ILINE(ML), LCON(3, ML), LTYPE(ML), N(29)
      DIMENSION LINKL(2, ML), LINKP(2, MP)

      LOGICAL BIFIND, BAD, ADDLNK, ADDCEN, ERR

      PI = ATAN2(0.0, -1.0)
      TWOPI = PI + PI

C  FIND THE CLOSEST LINE ABOVE THE POINT INPUT

      BIFIND = .FALSE.
      ADDLNK = .FALSE.
      IFIND = 0
      DIST = 100000.
      DO 100 I = 1, N(19)
         CALL LTSORT (ML, LINKL, I, II, ADDLNK)
         IF (II .GT. 0) THEN
            KT = LTYPE(II)
            CALL LINEPR (ML, MP, LINKP, LCON, II, I1, I2, I3, J1, J2,
     &         J3)
            IF ((J1 .GT. 0) .AND. (J2 .GT. 0)) THEN
               IF (KT .EQ. 1) THEN

C  GET THE PARAMETERS FOR THE LINE

                  CALL DLPARA (COOR(1, J1), COOR(2, J1), COOR(1, J2),
     &               COOR(2, J2), XM1, B1, BAD)

C  GET DISTANCE FOR VERTICAL LINE

                  IF (BAD) THEN
                     DTRY = ABS(COOR(1, J1) - X)
                     XTRY = COOR(1, J1)
                     YTRY = Y

C  GET DISTANCE FOR HORIZONTAL LINE

                  ELSE IF (ABS(XM1) .LT. .000001) THEN
                     DTRY = ABS(COOR(2, J1) - Y)
                     XTRY = X
                     YTRY = COOR(2, J1)

C  GET PERPENDICULAR DISTANCE TO ARBITRARY LINE

                  ELSE
                     XM2 = -1./XM1
                     B2 = Y - (XM2*X)
                     XTRY = (B2 - B1)/(XM1 - XM2)
                     YTRY = (XM1*XTRY) + B1
                     DTRY = SQRT((X - XTRY)**2 + (Y - YTRY)**2)
                  END IF
                  IF (DTRY .LE. DIST) THEN
                     X1 = MIN(COOR(1, J1), COOR(1, J2))
                     X2 = MAX(COOR(1, J1), COOR(1, J2))
                     Y1 = MIN(COOR(2, J1), COOR(2, J2))
                     Y2 = MAX(COOR(2, J1), COOR(2, J2))
                     IF ((XTRY .GE. X1) .AND. (XTRY .LE. X2) .AND.
     &                  (YTRY .GE. Y1) .AND. (YTRY .LE. Y2)) THEN
                        DIST = DTRY
                        XHOLD = XTRY
                        YHOLD = YTRY
                        IFIND = I
                        BIFIND = .TRUE.
                     END IF
                  END IF

C  CHECK DISTANCES TO CIRCULAR ARCS

               ELSE IF ((KT .EQ. 3).OR.(KT .EQ. 4).OR.(KT .EQ. 6)) THEN

C  FIRST GET THETA1, THETA2, THETAT, R1, R2, AND RTRY

                  CALL ARCPAR (MP, KT, ILINE(II), COOR, LINKP, J1, J2,
     &               J3, I3, XCEN, YCEN, THETA1, THETA2, TANG, R1, R2,
     &               ERR, ICCW, ICW, XK, XA)

                  IF ((Y .EQ. YCEN) .AND. (X .EQ. XCEN)) THEN
                     RTRY = 0.
                     THETAT = (THETA1 + THETA2)*.5
                  ELSE
                     THETAT = ATAN2(Y - YCEN, X - XCEN)
                     RTRY = SQRT( ((X - XCEN)**2)  +  ((Y - YCEN)**2))
                  END IF

C  SEE IF THE POINT ANGLE IS WITHIN THE BEGINNING AND ENDING ANGLES

                  IF ( ((THETAT .LE. THETA2) .AND. (THETAT .GE. THETA1))
     &               .OR. ((THETAT + TWOPI .LE. THETA2) .AND.
     &               (THETAT + TWOPI .GE. THETA1)) ) THEN

C  CALCULATE THE ARC RADIUS AT THAT ANGLE

                     RADIUS = XA*EXP(XK*THETAT)
                     DTRY = ABS(RADIUS - RTRY)

C  CHECK TO SEE IF THE ARC IS THE CLOSEST

                     IF (DTRY .LE. DIST) THEN
                        DIST = DTRY
                        XHOLD = XCEN + COS(THETAT)*RADIUS
                        YHOLD = YCEN + SIN(THETAT)*RADIUS
                        IFIND = I
                        IF ((KT .EQ. 4).OR.(KT .EQ. 6)) THEN
                           ADDCEN = .TRUE.
                           XCHOLD = XCEN
                           YCHOLD = YCEN
                        ELSE
                           ADDCEN = .FALSE.
                        END IF
                        BIFIND = .TRUE.
                     END IF
                  END IF
               END IF
            END IF
         END IF
  100 CONTINUE
      X = XHOLD
      Y = YHOLD

      RETURN
      END
