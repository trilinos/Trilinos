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

C $Id: closel.f,v 1.2 1991/03/21 15:44:25 gdsjaar Exp $
C $Log: closel.f,v $
C Revision 1.2  1991/03/21 15:44:25  gdsjaar
C Changed all 3.14159... to atan2(0.0, -1.0)
C
c Revision 1.1.1.1  1990/11/30  11:05:01  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:04:59  gdsjaar
c Initial revision
c 
C
CC* FILE: [.MAIN]CLOSEL.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE CLOSEL (MP, ML, N, COOR, ILINE, LTYPE, LCON, LINKP,
     &   LINKL, X, Y, BIFIND, IFIND, ADDCEN, XCHOLD, YCHOLD)
C***********************************************************************
C
C  SUBROUTINE CLOSEL = FINDS CLOSEST PERPENDICULAR BISECTED LINE
C
C***********************************************************************
C
C  SUBROUTINE CALLED BY:
C     BISECT   =
C
C***********************************************************************
C
C  SUBROUTINES CALLED:
C     DLPARA = DETERMINES LINE PARAMETERS FROM TWO POINTS
C
C***********************************************************************
C
      DIMENSION COOR(2, MP), ILINE(ML), LCON(3, ML), LTYPE(ML), N(29)
      DIMENSION LINKL(2, ML), LINKP(2, MP)
C
      LOGICAL BIFIND, BAD, ADDLNK, ADDCEN, ERR
C
      PI = ATAN2(0.0, -1.0)
      TWOPI = PI + PI
C
C  FIND THE CLOSEST LINE ABOVE THE POINT INPUT
C
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
C
C  GET THE PARAMETERS FOR THE LINE
C
                  CALL DLPARA (COOR(1, J1), COOR(2, J1), COOR(1, J2),
     &               COOR(2, J2), XM1, B1, BAD)
C
C  GET DISTANCE FOR VERTICAL LINE
C
                  IF (BAD) THEN
                     DTRY = ABS(COOR(1, J1) - X)
                     XTRY = COOR(1, J1)
                     YTRY = Y
C
C  GET DISTANCE FOR HORIZONTAL LINE
C
                  ELSE IF (ABS(XM1) .LT. .000001) THEN
                     DTRY = ABS(COOR(2, J1) - Y)
                     XTRY = X
                     YTRY = COOR(2, J1)
C
C  GET PERPENDICULAR DISTANCE TO ARBITRARY LINE
C
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
C
C  CHECK DISTANCES TO CIRCULAR ARCS
C
               ELSE IF ((KT .EQ. 3).OR.(KT .EQ. 4).OR.(KT .EQ. 6)) THEN
C
C  FIRST GET THETA1, THETA2, THETAT, R1, R2, AND RTRY
C
                  CALL ARCPAR (MP, KT, ILINE(II), COOR, LINKP, J1, J2,
     &               J3, I3, XCEN, YCEN, THETA1, THETA2, TANG, R1, R2,
     &               ERR, ICCW, ICW, XK, XA)
C
                  IF ((Y .EQ. YCEN) .AND. (X .EQ. XCEN)) THEN
                     RTRY = 0.
                     THETAT = (THETA1 + THETA2)*.5
                  ELSE
                     THETAT = ATAN2(Y - YCEN, X - XCEN)
                     RTRY = SQRT( ((X - XCEN)**2)  +  ((Y - YCEN)**2))
                  END IF
C
C  SEE IF THE POINT ANGLE IS WITHIN THE BEGINNING AND ENDING ANGLES
C
                  IF ( ((THETAT .LE. THETA2) .AND. (THETAT .GE. THETA1))
     &               .OR. ((THETAT + TWOPI .LE. THETA2) .AND.
     &               (THETAT + TWOPI .GE. THETA1)) ) THEN
C
C  CALCULATE THE ARC RADIUS AT THAT ANGLE
C
                     RADIUS = XA*EXP(XK*THETAT)
                     DTRY = ABS(RADIUS - RTRY)
C
C  CHECK TO SEE IF THE ARC IS THE CLOSEST
C
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
C
      RETURN
      END
