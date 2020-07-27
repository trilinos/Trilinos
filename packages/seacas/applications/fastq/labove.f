C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE LABOVE (MP,ML,N,IPOINT,COOR,ILINE,LTYPE,LCON,
     +   LINKP,LINKL,X,Y,Y1,Y2,IFIND,INDEX,IFIRST,INEXT)
C***********************************************************************

C  SUBROUTINE LABOVE = GETS THE CLOSEST LINE ABOVE A GIVEN X,Y

C***********************************************************************

C  VARIABLES USED:
C     IFIND  = THE CLOSEST LINE FOUND ABOVE THE X,Y LOCATION
C     IFIRST = THE CCW END OF THE LINE (LEFT END)
C     INEXT  = THE CW END OF THE LINE (RIGHT END)

C***********************************************************************

      DIMENSION IPOINT(MP),COOR(2,MP),ILINE(ML),LTYPE(ML),LCON(3,ML)
      DIMENSION LINKP(2,MP),LINKL(2,ML)
      DIMENSION N(29)

      LOGICAL BAD,ADDLNK,ERR,UP

      PI = ATAN2(0.0, -1.0)

      TWOPI=PI+PI
      ADDLNK=.FALSE.
      DIST=Y2-Y1

      DO 100 I=1,N(19)
         CALL LTSORT(ML,LINKL,I,II,ADDLNK)
         IF(II.GT.0)THEN
            I1=LCON(1,II)
            I2=LCON(2,II)
            I3=LCON(3,II)
            CALL LTSORT(MP,LINKP,I1,J1,ADDLNK)
            CALL LTSORT(MP,LINKP,I2,J2,ADDLNK)
            IF(I3.NE.0)THEN
               CALL LTSORT(MP,LINKP,IABS(I3),J3,ADDLNK)
            ELSE
               J3=0
            ENDIF

C  SEE IF LINE EXISTS

            IF((J1.GT.0).AND.(J2.GT.0))THEN

C  CHECK A STRAIGHT LINE TO SEE IF IT SPANS THE POINT

               IF((LTYPE(II).EQ.1).AND.
     +            (((COOR(1,J1).LE.X).AND.(COOR(1,J2).GE.X)).OR.
     +            ((COOR(1,J2).LE.X).AND.(COOR(1,J1).GE.X))))THEN

C  SEE IF LINE IS ABOVE

                  CALL DLPARA(COOR(1,J1),COOR(2,J1),COOR(1,J2),
     +               COOR(2,J2),XM1,B1,BAD)
                  IF(.NOT.BAD)THEN
                     YTRY=(XM1*X)+B1
                     IF(YTRY.GT.Y)THEN

C  CHECK DISTANCE TO THE LINE - RECORD LINE NO. IF CLOSEST

                        DTRY=YTRY-Y
                        IF((DTRY.LE.DIST).OR.(ABS(DIST-DTRY).LT.
     +                     .001*DIST))THEN
                           DIST=DTRY
                           IFIND=I
                           INDEX=II
                           IF(COOR(1,J1).GT.X)THEN
                              IFIRST=I2
                              INEXT=I1
                           ELSE
                              IFIRST=I1
                              INEXT=I2
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF

C  CHECK AN ARC LINE

               ELSEIF(((LTYPE(II).EQ.3).OR.(LTYPE(II).EQ.4).OR.
     +            (LTYPE(II).EQ.6)).AND.(J3.GT.0))THEN
                  CALL ARCPAR(MP,LTYPE(II),ILINE(II),COOR,LINKP,J1,J2,
     +               J3,I3,XCEN,YCEN,THETA1,THETA2,TANG,R1,R2,ERR,
     +               ICCW,ICW,XK,XA)
                  IF(.NOT.ERR)THEN

C  SET THE ARC AS EITHER AN UPWARD ARCH OR DOWNWARD ARCH
C  (A CLOSED CIRCLE IS ALWAYS AN UPWARD ARCH)

                     IF(J1.EQ.J2)THEN
                        UP=.TRUE.
                     ELSEIF(COOR(1,ICCW).GE.COOR(1,ICW))THEN
                        UP=.TRUE.
                     ELSE
                        UP=.FALSE.
                     ENDIF

C  THE INPUT POINT IS RIGHT AT THE CENTER OF THE CIRCLE

                     IF((Y.EQ.YCEN).AND.(X.EQ.XCEN))THEN
                        RP=0.
                        THETAP=(THETA1+THETA2)*.5
                     ELSE
                        THETAP=ATAN2(Y-YCEN,X-XCEN)
                        RP=SQRT( ((X-XCEN)**2) + ((Y-YCEN)**2) )
                     ENDIF

C  SEE IF THE POINT ANGLE IS WITHIN THE BEGINNING AND ENDING ANGLES

                     IF( ( (THETAP.LE.THETA2) .AND. (THETAP.GE.THETA1) )
     &                  .OR.
     &                  ( (THETAP+TWOPI.LE.THETA2) .AND.
     &                  (THETAP+TWOPI.GE.THETA1) )
     &                  .OR.
     &                  ( (THETAP.LE.THETA1) .AND. (THETAP.GE.THETA2) )
     &                  .OR.
     &                  ( (THETAP+TWOPI.LE.THETA1) .AND.
     &                  (THETAP+TWOPI.GE.THETA2) ) ) THEN

C  SEE IF THE POINT TO CENTER DISTANCE IS LESS THAN THE
C  ARC RADIUS AT THAT ANGLE FOR AN UPWARD ARC OR GREATER FOR
C  A DOWNWARD ARC (BELOW THE ARC)

                        RTEST = XA * EXP ( XK * THETAP )
                        IF( ( (UP) .AND. (RTEST.GE.RP) ) .OR.
     +                     ( (.NOT.UP) .AND. (RTEST.LE.RP) ) )THEN

C  CHECK Y DISTANCE TO THE LINE - RECORD LINE NO. IF CLOSEST

                           CALL ARCY(XCEN,YCEN,THETA1,THETA2,XK,XA,
     +                        X,YTRY,ERR)
                           DTRY=YTRY-Y
                           IF( (.NOT.ERR) .AND.
     +                        (YTRY .GT. Y) .AND.
     +                        ( (DTRY.LE.DIST) .OR.
     +                        (ABS(DIST-DTRY).LT..001*DIST)))THEN
                              DIST=DTRY
                              IFIND=I
                              INDEX=II
                              IF(UP)THEN
                                 IFIRST=IPOINT(ICW)
                                 INEXT=IPOINT(ICCW)
                              ELSE
                                 IFIRST=IPOINT(ICCW)
                                 INEXT=IPOINT(ICW)
                              ENDIF
                           ENDIF
                        ENDIF

C  THE ONLY OTHER ARC POSSIBILITY IS IF THE X FALLS IN THE SPAN
C  BETWEEN THE TWO ENDPOINTS OR BETWEEN END POINTS AND CENTER POINT

                     ELSEIF ( ((X - COOR(1,J1) ) *
     +                  (X - COOR(1,J2) ) .LT.0)
     +                  .OR.
     +                  ((X - COOR(1,J1) ) *
     +                  (X - COOR(1,J3) ) .LT.0)
     +                  .OR.
     +                  ((X - COOR(1,J2) ) *
     +                  (X - COOR(1,J3) ) .LT.0) ) THEN

C  CHECK Y DISTANCE TO THE LINE - RECORD LINE NO. IF CLOSEST

                        CALL ARCY(XCEN,YCEN,THETA1,THETA2,XK,XA,
     +                     X,YTRY,ERR)
                        DTRY=YTRY-Y
                        IF( (.NOT.ERR) .AND.
     +                     (YTRY .GT. Y) .AND.
     +                     ( (DTRY.LE.DIST) .OR.
     +                     (ABS(DIST-DTRY).LT..001*DIST)))THEN
                           DIST=DTRY
                           IFIND=I
                           INDEX=II

C  TREAT THE BETWEEN END POINTS DIFFERENTLY THAN BETWEEN
C  ENDPOINT AND CENTER POINT

                           IF( ((X - COOR(1,J1) ) *
     +                        (X - COOR(1,J2) )) .LT.0) THEN
                              IF(UP)THEN
                                 IFIRST=IPOINT(ICW)
                                 INEXT=IPOINT(ICCW)
                              ELSE
                                 IFIRST=IPOINT(ICCW)
                                 INEXT=IPOINT(ICW)
                              ENDIF
                           ELSE
                              IF( COOR(2,ICCW) .GT.
     +                           COOR(2,ICW) ) THEN
                                 IFIRST=IPOINT(ICCW)
                                 INEXT=IPOINT(ICW)
                              ELSE
                                 IFIRST=IPOINT(ICW)
                                 INEXT=IPOINT(ICCW)
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
  100 CONTINUE

      RETURN

      END
