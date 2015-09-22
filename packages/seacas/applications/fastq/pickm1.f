C    Copyright (c) 2014, Sandia Corporation.
C    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C    the U.S. Government retains certain rights in this software.
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C        * Redistributions of source code must retain the above copyright
C          notice, this list of conditions and the following disclaimer.
C    
C        * Redistributions in binary form must reproduce the above
C          copyright notice, this list of conditions and the following
C          disclaimer in the documentation and/or other materials provided
C          with the distribution.
C    
C        * Neither the name of Sandia Corporation nor the names of its
C          contributors may be used to endorse or promote products derived
C          from this software without specific prior written permission.
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

C $Id: pickm1.f,v 1.2 1991/03/21 15:44:59 gdsjaar Exp $
C $Log: pickm1.f,v $
C Revision 1.2  1991/03/21 15:44:59  gdsjaar
C Changed all 3.14159... to atan2(0.0, -1.0)
C
c Revision 1.1.1.1  1990/11/30  11:13:24  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:13:23  gdsjaar
c Initial revision
c 
C
CC* FILE: [.QMESH]PICKM1.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE PICKM1 (N, X, Y, ANGLE, M1, IFIRST, REAL)
C***********************************************************************
C
C  SUBROUTINE PICKM1 = DETERMINES A REASONABLE SHAPE FOR A LOGICAL
C                      RECTANGLE WITH PERIMETER GIVEN IN X AND Y
C
C***********************************************************************
C
      DIMENSION X (N), Y (N), ANGLE (N)
      DIMENSION SMANG (7), INDEX (7)
C
      LOGICAL REAL
C
C  FORM THE LIST OF SMALLEST ANGLES
C
      NSA = 6
      DO 100 I = 1, NSA
         SMANG (I) = 10.
         INDEX (I) = 0
  100 CONTINUE
C
      PI = ATAN2(0.0, -1.0)
      TWOPI = PI + PI
      AGOLD = ATAN2 (Y (1) - Y (N), X (1) - X (N))
C
      DO 130 J = 1, N
C
C  GET THE ANGLE FORMED BY THIS SET OF POINTS
C
         NEXT = J + 1
         IF  (NEXT .GT. N) NEXT = 1
         AGNEW = ATAN2 (Y (NEXT) - Y (J) ,  X (NEXT) - X (J))
         DIFF = AGNEW - AGOLD
         IF (DIFF .GT. PI)DIFF = DIFF - TWOPI
         IF (DIFF .LT.  - PI)DIFF = DIFF + TWOPI
         ANGLE (J) = PI - DIFF
         AGOLD = AGNEW
C
C  SORT THIS ANGLE AGAINST PREVIOUS ANGLES TO SEE IF IT IS ONE OF
C  THE SMALLEST
C
         SMANG (NSA + 1) = ANGLE (J)
         INDEX (NSA + 1) = J
         DO 110 II = 1, NSA
            I = NSA + 1 - II
            IF  (SMANG (I + 1) .GE. SMANG (I)) GO TO 120
            TEMP = SMANG (I)
            ITEMP = INDEX (I)
            SMANG (I) = SMANG (I + 1)
            INDEX (I) = INDEX (I + 1)
            SMANG (I + 1) = TEMP
  110       INDEX (I + 1) = ITEMP
  120    CONTINUE
C
  130 CONTINUE
C
C  DETERMINE OPTIMUM ORIGIN / SHAPE COMBINATION
C
      ATOL = PI * 150. / 180.
      IFIRST = 1
      M1 = N / 4
      M2 = N / 2 - M1
      I2 = 1 + M1
      I3 = I2 + M2
      I4 = I3 + M1
      GBEST = ANGLE (1) + ANGLE (I2) + ANGLE (I3) + ANGLE (I4)
      BADANG = AMAX1 (ANGLE (1), ANGLE (I2), ANGLE (I3), ANGLE (I4))
C
      MMAX = N / 2 - 1
      AMAXEL = FLOAT (N / 4) * FLOAT ( (N + 2) / 4)
      DO 150 ISA = 1, NSA
         IF (SMANG (ISA) .LE. ATOL) THEN
            I1 = INDEX (ISA)
            DO 140 M = 1, MMAX
               M2 = N / 2 - M
               I2 = I1 + M
               IF  (I2 .GT. N) I2 = I2 - N
               I3 = I2 + M2
               IF  (I3 .GT. N) I3 = I3 - N
               I4 = I3 + M
               IF  (I4 .GT. N) I4 = I4 - N
               AFAC = ANGLE (I1) + ANGLE (I2) + ANGLE (I3) + ANGLE (I4)
               ERAT = AMIN1 (AMAXEL / FLOAT (M * M2) ,  5.)
               EFAC =  (ERAT + 15.) / 16.
               GVAL = AFAC * EFAC
               IF (GVAL .LT. GBEST) THEN
                  BADANG = AMAX1 (ANGLE (I1), ANGLE (I2), ANGLE (I3),
     &               ANGLE (I4))
                  IFIRST = I1
                  M1 = M
                  GBEST = GVAL
               ENDIF
  140       CONTINUE
         ENDIF
  150 CONTINUE
      IF ( (REAL) .AND. (BADANG .GT. 2.62)) THEN
         CALL MESAGE (' **  WARNING: CORNER (S) OF THE REGION HAVE  **')
         CALL MESAGE (' **           LARGE ANGLES  (> 150 DEGREES.) **')
         CALL MESAGE (' **           POORLY FORMED MESH MAY RESULT  **')
      ENDIF
C
      RETURN
      END
