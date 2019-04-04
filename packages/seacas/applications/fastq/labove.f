C    Copyright(C) 2014-2017 National Technology & Engineering Solutions of
C    Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C
C    * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C
C    * Redistributions in binary form must reproduce the above
C      copyright notice, this list of conditions and the following
C      disclaimer in the documentation and/or other materials provided
C      with the distribution.
C
C    * Neither the name of NTESS nor the names of its
C      contributors may be used to endorse or promote products derived
C      from this software without specific prior written permission.
C
C    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C

C $Id: labove.f,v 1.2 1991/03/21 15:44:51 gdsjaar Exp $
C $Log: labove.f,v $
C Revision 1.2  1991/03/21 15:44:51  gdsjaar
C Changed all 3.14159... to atan2(0.0, -1.0)
C
c Revision 1.1.1.1  1990/11/30  11:10:57  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:10:56  gdsjaar
c Initial revision
c
C
CC* FILE: [.MAIN]LABOVE.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE LABOVE (MP,ML,N,IPOINT,COOR,ILINE,LTYPE,LCON,
     +   LINKP,LINKL,X,Y,Y1,Y2,IFIND,INDEX,IFIRST,INEXT)
C***********************************************************************
C
C  SUBROUTINE LABOVE = GETS THE CLOSEST LINE ABOVE A GIVEN X,Y
C
C***********************************************************************
C
C  VARIABLES USED:
C     IFIND  = THE CLOSEST LINE FOUND ABOVE THE X,Y LOCATION
C     IFIRST = THE CCW END OF THE LINE (LEFT END)
C     INEXT  = THE CW END OF THE LINE (RIGHT END)
C
C***********************************************************************
C
      DIMENSION IPOINT(MP),COOR(2,MP),ILINE(ML),LTYPE(ML),LCON(3,ML)
      DIMENSION LINKP(2,MP),LINKL(2,ML)
      DIMENSION N(29)
C
      LOGICAL BAD,ADDLNK,ERR,UP
C
      PI = ATAN2(0.0, -1.0)
C
      TWOPI=PI+PI
      ADDLNK=.FALSE.
      DIST=Y2-Y1
C
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
C
C  SEE IF LINE EXISTS
C
            IF((J1.GT.0).AND.(J2.GT.0))THEN
C
C  CHECK A STRAIGHT LINE TO SEE IF IT SPANS THE POINT
C
               IF((LTYPE(II).EQ.1).AND.
     +            (((COOR(1,J1).LE.X).AND.(COOR(1,J2).GE.X)).OR.
     +            ((COOR(1,J2).LE.X).AND.(COOR(1,J1).GE.X))))THEN
C
C  SEE IF LINE IS ABOVE
C
                  CALL DLPARA(COOR(1,J1),COOR(2,J1),COOR(1,J2),
     +               COOR(2,J2),XM1,B1,BAD)
                  IF(.NOT.BAD)THEN
                     YTRY=(XM1*X)+B1
                     IF(YTRY.GT.Y)THEN
C
C  CHECK DISTANCE TO THE LINE - RECORD LINE NO. IF CLOSEST
C
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
C
C  CHECK AN ARC LINE
C
               ELSEIF(((LTYPE(II).EQ.3).OR.(LTYPE(II).EQ.4).OR.
     +            (LTYPE(II).EQ.6)).AND.(J3.GT.0))THEN
                  CALL ARCPAR(MP,LTYPE(II),ILINE(II),COOR,LINKP,J1,J2,
     +               J3,I3,XCEN,YCEN,THETA1,THETA2,TANG,R1,R2,ERR,
     +               ICCW,ICW,XK,XA)
                  IF(.NOT.ERR)THEN
C
C  SET THE ARC AS EITHER AN UPWARD ARCH OR DOWNWARD ARCH
C  (A CLOSED CIRCLE IS ALWAYS AN UPWARD ARCH)
C
                     IF(J1.EQ.J2)THEN
                        UP=.TRUE.
                     ELSEIF(COOR(1,ICCW).GE.COOR(1,ICW))THEN
                        UP=.TRUE.
                     ELSE
                        UP=.FALSE.
                     ENDIF
C
C  THE INPUT POINT IS RIGHT AT THE CENTER OF THE CIRCLE
C
                     IF((Y.EQ.YCEN).AND.(X.EQ.XCEN))THEN
                        RP=0.
                        THETAP=(THETA1+THETA2)*.5
                     ELSE
                        THETAP=ATAN2(Y-YCEN,X-XCEN)
                        RP=SQRT( ((X-XCEN)**2) + ((Y-YCEN)**2) )
                     ENDIF
C
C  SEE IF THE POINT ANGLE IS WITHIN THE BEGINNING AND ENDING ANGLES
C
                     IF( ( (THETAP.LE.THETA2) .AND. (THETAP.GE.THETA1) )
     &                  .OR.
     &                  ( (THETAP+TWOPI.LE.THETA2) .AND.
     &                  (THETAP+TWOPI.GE.THETA1) )
     &                  .OR.
     &                  ( (THETAP.LE.THETA1) .AND. (THETAP.GE.THETA2) )
     &                  .OR.
     &                  ( (THETAP+TWOPI.LE.THETA1) .AND.
     &                  (THETAP+TWOPI.GE.THETA2) ) ) THEN
C
C  SEE IF THE POINT TO CENTER DISTANCE IS LESS THAN THE
C  ARC RADIUS AT THAT ANGLE FOR AN UPWARD ARC OR GREATER FOR
C  A DOWNWARD ARC (BELOW THE ARC)
C
                        RTEST = XA * EXP ( XK * THETAP )
                        IF( ( (UP) .AND. (RTEST.GE.RP) ) .OR.
     +                     ( (.NOT.UP) .AND. (RTEST.LE.RP) ) )THEN
C
C  CHECK Y DISTANCE TO THE LINE - RECORD LINE NO. IF CLOSEST
C
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
C
C  THE ONLY OTHER ARC POSSIBILITY IS IF THE X FALLS IN THE SPAN
C  BETWEEN THE TWO ENDPOINTS OR BETWEEN END POINTS AND CENTER POINT
C
                     ELSEIF ( ((X - COOR(1,J1) ) *
     +                  (X - COOR(1,J2) ) .LT.0)
     +                  .OR.
     +                  ((X - COOR(1,J1) ) *
     +                  (X - COOR(1,J3) ) .LT.0)
     +                  .OR.
     +                  ((X - COOR(1,J2) ) *
     +                  (X - COOR(1,J3) ) .LT.0) ) THEN
C
C  CHECK Y DISTANCE TO THE LINE - RECORD LINE NO. IF CLOSEST
C
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
C
C  TREAT THE BETWEEN END POINTS DIFFERENTLY THAN BETWEEN
C  ENDPOINT AND CENTER POINT
C
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
C
      RETURN
C
      END
