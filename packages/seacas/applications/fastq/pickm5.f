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

C $Id: pickm5.f,v 1.2 1991/03/21 15:45:03 gdsjaar Exp $
C $Log: pickm5.f,v $
C Revision 1.2  1991/03/21 15:45:03  gdsjaar
C Changed all 3.14159... to atan2(0.0, -1.0)
C
c Revision 1.1.1.1  1990/11/30  11:13:30  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:13:28  gdsjaar
c Initial revision
c 
C
CC* FILE: [.QMESH]PICKM5.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE PICKM5 (N, X, Y, ANGLE, IST2, IST3, IST4, IST5, INDST,
     &   M1, N1, N2, N3, N4)
C***********************************************************************
C
C  SUBROUTINE PICKM5 = DETERMINES A REASONABLE SHAPE FOR A LOGICAL
C                      PENTAGON WITH PERIMETER GIVEN IN X AND Y
C
C***********************************************************************
C
C  VARIABLES IN: N....... NUMBER OF POINTS
C                X, Y.... ARRAY OF COORDINATES
C                ANGLE... ARRAY OF ANGLES
C           OUT: Ni...... POINTERS TO THE CHOSEN CORNERS
C
C
C  WRITTEN BY: HORACIO RECALDE                       DATE: JAN, 1988
C
C**********************************************************************
C
      PARAMETER (NSA = 10, NSA2 = 12)
      PARAMETER (BIGNUM = 99999.0, MAXTRY = 10)
C
      DIMENSION X(N), Y(N), ANGLE(N), IST1(NSA2), IST2(N), IST3(N)
      DIMENSION IST4(N), IST5(N), INDEX(NSA + 1), INDST(N)
      DIMENSION SMANG(NSA + 1)
      LOGICAL CHECK
C
C--- SET THE MAXIMUM SIDE AND TOLERANCE
C
      PI = ATAN2(0.0, -1.0)
      MMAX = N/2 - 1
      TOL = 150./180.*PI
C
C---GET THE 10 SMALLEST ANGLES
C
      CALL CSMANG (N, X, Y, ANGLE, NSA, SMANG, INDEX)
      IF (ANGLE (INDEX (1)) .GT. TOL) THEN
         M1 = -1
         RETURN
      END IF
C
C---CHECK CONDITIONS ON THE 5 SMALLEST ANGLES
C
      CALL CHCOND (N, NSA, SMANG, INDEX, M1, N1, N2, N3, N4, CHECK)
      IF (CHECK) RETURN
C
      SUMANG = BIGNUM
C
C---PUSH ANGLE INDEX ON FIRST STACK
C
      CALL SKINIT (IST1, NSA2, NSA, IERROR)
      DO 100 I = NSA, 1, -1
         CALL SKPUSH (IST1, NSA2, INDEX (I), IERROR)
  100 CONTINUE
C
C---POP TOP OF FIRST STACK
C
  110 CONTINUE
      CALL SKPOP (IST1, NSA2, M1, IERROR)
      IF (IERROR .EQ. 2) GO TO 200
      IF (ANGLE (M1) .GT. SUMANG) GO TO 200
      SUM1 = ANGLE (M1)
C
C---SELECT SECOND TRIAL
C
      MM1 = M1 + 2
      IF (MM1 .GT. N) MM1 = MM1 - N
      NN1 = MIN (M1 + MMAX, M1 + N - 8)
      IF (NN1 .GT. N) NN1 = NN1 - N
      CALL SORTST (N, ANGLE, MM1, NN1, IRANGE, INDST)
C
C---PUSH ANGLE INDEX ON SECOND STACK IN DECENDING ORDER
C
      CALL SKINIT (IST2, N, N - 2, IERROR)
      DO 120 I = MIN (MAXTRY, IRANGE), 1, -1
         CALL SKPUSH (IST2, N, INDST (I), IERROR)
  120 CONTINUE
C
C---POP TOP OF SECOND STACK
C
  130 CONTINUE
      CALL SKPOP (IST2, N, M2, IERROR)
      IF (IERROR .EQ. 2) GO TO 110
      IF (SUM1 + ANGLE (M2) .GE. SUMANG) GO TO 110
      SUM2 = SUM1 + ANGLE (M2)
      N1 = M2 - M1
      IF (N1 .LT. 0) N1 = N + N1
C
C---SELECT THIRD TRIAL NODE AND SORT IN ASCENDING ORDER
C
      MM2 = M2 + 2
      IF (MM2 .GT. N) MM2 = MM2 - N
      NN2 = MIN (M1 + N1 + MMAX, M1 + N - 6)
      IF (NN2 .GT. N) NN2 = NN2 - N
      CALL SORTST (N, ANGLE, MM2, NN2, IRANGE, INDST)
C
C---PUSH ANGLE INDEX ON THIRD STACK IN DECENDING ORDER
C
      CALL SKINIT (IST3, N, N - 2, IERROR)
      DO 140 I = MIN (MAXTRY, IRANGE), 1, -1
         CALL SKPUSH (IST3, N, INDST (I), IERROR)
  140 CONTINUE
C
C---POP TOP OF THIRD STACK
C
  150 CONTINUE
      CALL SKPOP (IST3, N, M3, IERROR)
      IF (IERROR .EQ. 2) GO TO 130
      IF (SUM2 + ANGLE (M3) .GE. SUMANG) GO TO 130
      SUM3 = SUM2 + ANGLE (M3)
      N2 = M3 - M2
      IF (N2 .LT. 0) N2 = N + N2
      IF (N1 + N2 .GT. N/2 -1) GO TO 150
C
C---SELECT FOURTH TRIAL AND SORT IN ASCENDING ORDER
C
      MM3 = M3 + 2
      IF (MM3 .GT. N) MM3 = MM3 - N
      NN3 = MIN (M1 + N1 + N2 + MMAX, M1 + N - 4)
      IF (NN3 .GT. N) NN3 = NN3 - N
      CALL SORTST (N, ANGLE, MM3, NN3, IRANGE, INDST)
C
C---PUSH ANGLE INDEX ON FOURTH STACK IN DECENDING ORDER
C
      CALL SKINIT (IST4, N, N - 2, IERROR)
      DO 160 I = MIN (MAXTRY, IRANGE), 1, -1
         CALL SKPUSH (IST4, N, INDST (I), IERROR)
  160 CONTINUE
C
C---POP TOP OF FOURTH STACK
C
  170 CONTINUE
      CALL SKPOP (IST4, N, M4, IERROR)
      IF (IERROR .EQ. 2) GO TO 150
      IF (SUM3 + ANGLE (M4) .GE. SUMANG) GO TO 150
      SUM4 = SUM3 + ANGLE (M4)
      N3 = M4 - M3
      IF (N3 .LT. 0) N3 = N + N3
      IF (N2 + N3 .GT. N/2 -1) GO TO 170
C
C---SELECT FIFTH TRIAL AND SORT IN ASCENDING ORDER
C
      MM4 = M4 + 2
      IF (MM4 .GT. N) MM4 = MM4 - N
      NN4 = MIN (M1 + N1 + N2 + N3 + MMAX, M1 + N - 2)
      IF (NN4 .GT. N) NN4 = NN4 - N
      CALL SORTST (N, ANGLE, MM4, NN4, IRANGE, INDST)
C
C---PUSH ANGLE INDEX ON FIFTH STACK IN DECENDING ORDER
C
      CALL SKINIT (IST5, N, N - 2, IERROR)
      DO 180 I = MIN (MAXTRY, IRANGE), 1, -1
         CALL SKPUSH (IST5, N, INDST (I), IERROR)
  180 CONTINUE
C
C---POP TOP OF FIFTH STACK
C
  190 CONTINUE
      CALL SKPOP (IST5, N, M5, IERROR)
      IF (IERROR .EQ. 2) GO TO 170
      IF (SUM4 + ANGLE (M5) .GT. SUMANG) GO TO 170

      N4 = M5 - M4
      IF (N4 .LT. 0) N4 = N + N4
      N5 = M1 - M5
      IF (N5 .LT. 0) N5 = N + N5
C
C--- CHECK COMPATIBILITY EQUATIONS
C
      IF ( (N1 + N2 + N3 .LT. N4 + N5 + 2) .OR.
     &   (N2 + N3 + N4 .LT. N5 + N1 + 2) .OR.
     &   (N3 + N4 + N5 .LT. N1 + N2 + 2) .OR.
     &   (N4 + N5 + N1 .LT. N2 + N3 + 2) .OR.
     &   (N5 + N1 + N2 .LT. N3 + N4 + 2)) THEN
      ELSE
         SUMANG = SUM4 + ANGLE (M5)
         M1HOLD = M1
         N1HOLD = N1
         N2HOLD = N2
         N3HOLD = N3
         N4HOLD = N4
      END IF
      GO TO 190
C
  200 CONTINUE
      IF (SUMANG .EQ. BIGNUM) THEN
         M1 = 0
      ELSE
         M1 = M1HOLD
         N1 = N1HOLD
         N2 = N2HOLD
         N3 = N3HOLD
         N4 = N4HOLD
      END IF
C
      RETURN
      END
