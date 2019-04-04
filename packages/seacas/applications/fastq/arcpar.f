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

C $Id: arcpar.f,v 1.2 1991/03/21 15:44:18 gdsjaar Exp $
C $Log: arcpar.f,v $
C Revision 1.2  1991/03/21 15:44:18  gdsjaar
C Changed all 3.14159... to atan2(0.0, -1.0)
C
c Revision 1.1.1.1  1990/11/30  11:03:41  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:03:40  gdsjaar
c Initial revision
c
C
CC* FILE: [.MAIN]ARCPAR.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE ARCPAR (MP, KT, KNUM, COOR, LINKP, IPNTR1, IPNTR2,
     &   IPNTR3, IP3, XCEN, YCEN, THETA1, THETA2, TANG, R1, R2, ERR,
     &   ICCW, ICW, XK, XA)
C***********************************************************************
C
C  SUBROUTINE ARCPAR = THIS ROUTINE CALCULATES THE ARC PARAMETERS
C
C***********************************************************************
C
C  VARIABLES USED:
C     TANG   = TOTAL ANGLE SCRIBED BY THE ARC
C     THETA1 = FIRST CCW ANGLE OF THE ARC
C     THETA2 = SECOND CCW ANGLE OF THE ARC
C     IPNTR1 = POINTER TO FIRST COORDINATE VALUE
C     IPNTR2 = POINTER TO SECOND COORDINATE VALUE
C     IPNTR3 = POINTER TO THIRD COORDINATE VALUE
C     IP3    = THE THIRD POINT NUMBER (CAN BE NEGATED)
C
C***********************************************************************
C
      DIMENSION COOR (2, MP), LINKP (2, MP)
C
      LOGICAL ERR
C
      PI = ATAN2(0.0, -1.0)
      TWOPI = PI + PI
      ERR = .FALSE.
C
C  ARC WITH CENTER GIVEN
C  ARC GOES FROM 1ST POINT TO 2ND IN *COUNTER-CLOCKWISE* DIRECTION.
C
      IF (KT .EQ. 3) THEN
         XCEN = COOR (1, IPNTR3)
         YCEN = COOR (2, IPNTR3)
C
C  CIRCLE WITH THIRD POINT ON ARC
C
      ELSEIF (KT .EQ. 4) THEN
         THETA1 = ATAN2 (COOR (2, IPNTR3) - COOR (2, IPNTR1),
     &      COOR (1, IPNTR3) - COOR (1, IPNTR1)) + PI / 2.0
         THETA2 = ATAN2 (COOR (2, IPNTR3) - COOR (2, IPNTR2),
     &      COOR (1, IPNTR3) - COOR (1, IPNTR2)) + PI / 2.0
         DET = - COS (THETA1) * SIN (THETA2) +
     &      COS (THETA2) * SIN (THETA1)
         X11 = 0.5 * (COOR (1, IPNTR1) + COOR (1, IPNTR3))
         Y11 = 0.5 * (COOR (2, IPNTR1) + COOR (2, IPNTR3))
         X21 = 0.5 * (COOR (1, IPNTR2) + COOR (1, IPNTR3))
         Y21 = 0.5 * (COOR (2, IPNTR2) + COOR (2, IPNTR3))
         R =  ( - SIN (THETA2) * (X21 - X11) +
     &      COS (THETA2) * (Y21 - Y11)) / DET
         XCEN = X11 + R * COS (THETA1)
         YCEN = Y11 + R * SIN (THETA1)
C
C     CIRCLE WITH RADIUS GIVEN
C
      ELSEIF (KT .EQ. 6) THEN
         DX = 0.5 * (COOR (1, IPNTR2) - COOR (1, IPNTR1))
         DY = 0.5 * (COOR (2, IPNTR2) - COOR (2, IPNTR1))
         CHORD = SQRT (DX * DX + DY * DY)
         R = ABS (COOR (1, IPNTR3))
         IF (R .LE. CHORD) THEN
            XCEN = 0.5 * (COOR (1, IPNTR1) + COOR (1, IPNTR2))
            YCEN = 0.5 * (COOR (2, IPNTR1) + COOR (2, IPNTR2))
         ELSE
            ARM = SQRT (R * R - CHORD * CHORD)
            IF (IP3.LT.0) THEN
               XCEN = COOR (1, IPNTR1) + DX + ARM * DY / CHORD
               YCEN = COOR (2, IPNTR1) + DY - ARM * DX / CHORD
            ELSE
               XCEN = COOR (1, IPNTR1) + DX - ARM * DY / CHORD
               YCEN = COOR (2, IPNTR1) + DY + ARM * DX / CHORD
            ENDIF
         ENDIF
      ENDIF
C
C  CHECK TO MAKE SURE THAT THE BEGINNING AND ENDING RADIUS EXIST
C
      IF ((( COOR (1, IPNTR1) .EQ. XCEN) .AND.
     &   (COOR (2, IPNTR1) .EQ. YCEN)) .OR.
     &   ((COOR (1, IPNTR2) .EQ. XCEN) .AND.
     &   (COOR (2, IPNTR2) .EQ. YCEN))) THEN
         CALL PLTFLU
         WRITE (*, 10000) ABS (KNUM)
         ERR = .TRUE.
         GOTO 100
      ENDIF
      R1 = SQRT ( (COOR (1, IPNTR1) - XCEN) **2 +  (COOR (2, IPNTR1) -
     &   YCEN) **2)
      R2 = SQRT ( (COOR (1, IPNTR2) - XCEN) **2 +  (COOR (2, IPNTR2) -
     &   YCEN) **2)
      THETA1 = ATAN2 (COOR (2, IPNTR1) - YCEN, COOR (1, IPNTR1) - XCEN)
      THETA2 = ATAN2 (COOR (2, IPNTR2) - YCEN, COOR (1, IPNTR2) - XCEN)
C
C  ARC WITH THE CENTER GIVEN
C
      IF (KT .EQ. 3) THEN
         IF  (IPNTR1 .EQ. IPNTR2) THEN
            THETA2  =  THETA1  +  TWOPI
         ELSEIF ((IP3 .GE. 0) .AND. (THETA2 .LE. THETA1)) THEN
            THETA2 = THETA2 + TWOPI
         ELSEIF ((IP3 .LT. 0) .AND. (THETA1 .LE. THETA2)) THEN
            THETA1 = THETA1 + TWOPI
         ENDIF
         TANG = THETA2 - THETA1
         IF (IP3 .LT. 0) THEN
            ICCW = IPNTR2
            ICW = IPNTR1
         ELSE
            ICCW = IPNTR1
            ICW = IPNTR2
         ENDIF
C
C  CIRCULAR ARC WITH 3RD POINT ON ARC - CLOCKWISE OR COUNTER-CLOCKWISE
C
      ELSEIF (KT .EQ. 4) THEN
         THETA3 = ATAN2 (COOR (2, IPNTR3) - YCEN, COOR (1, IPNTR3) -
     &      XCEN)
         IF (THETA2 .LE. THETA1) THETA2 = THETA2 + TWOPI
         IF (THETA3 .LE. THETA1) THETA3 = THETA3 + TWOPI
         TANG = THETA2 - THETA1
         IF (THETA3 .GT. THETA2) THEN
            TANG = - (TWOPI - TANG)
            ICCW = IPNTR2
            ICW = IPNTR1
         ELSE
            ICCW = IPNTR1
            ICW = IPNTR2
         ENDIF
C
C     CIRRCULAR ARC WITH RADIUS GIVEN - CLOCKWISE OR COUNTER-CLOCKWISE
C
      ELSEIF (KT .EQ. 6) THEN
         IF ( (IP3 .GE. 0) .AND. (THETA2 .LE. THETA1))
     &      THETA2 = THETA2 + TWOPI
         IF ( (IP3 .LT. 0) .AND. (THETA1 .LE. THETA2))
     &      THETA1 = THETA1 + TWOPI
         TANG = THETA2 - THETA1
         IF (IP3 .GE. 0) THEN
            ICCW = IPNTR1
            ICW = IPNTR2
         ELSE
            ICCW = IPNTR2
            ICW = IPNTR1
         ENDIF
      ENDIF
C
C  DEFINE THE SPIRAL PARAMETERS  (R = XA * EXP  (XK * THETA))
C
      XK = (LOG (R2 / R1)) / (THETA2 - THETA1)
      DIVID = EXP (XK * THETA2)
      IF (DIVID .GT. 0.) THEN
         XA = R2 / DIVID
      ELSE
         WRITE (*, 10010) IABS (KNUM)
         ERR = .TRUE.
         GOTO 100
      ENDIF
C
  100 CONTINUE
C
      RETURN
C
10000 FORMAT (' CENTER POINT FOR LINE', I5, ' LIES ON ONE OF',
     &   ' THE ENDPOINTS')
10010 FORMAT (' DEFINITION FOR ARC LINE', I5, ' IS INVLAID')
C
      END
