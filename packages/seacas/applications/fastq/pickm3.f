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

C $Id: pickm3.f,v 1.3 1991/03/22 19:57:40 gdsjaar Exp $
C $Log: pickm3.f,v $
C Revision 1.3  1991/03/22 19:57:40  gdsjaar
C Guess: substitute MM1 for MM
C
c Revision 1.2  1991/03/21  15:45:01  gdsjaar
c Changed all 3.14159... to atan2(0.0, -1.0)
c
c Revision 1.1.1.1  1990/11/30  11:13:27  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:13:25  gdsjaar
c Initial revision
c 
C
CC* FILE: [.QMESH]PICKM3.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE PICKM3 (N, X, Y, ANGLE, M1, M2, IFIRST)
C***********************************************************************
C
C  SUBROUTINE PICKM3 = DETERMINES A REASONABLE SHAPE FOR A LOGICAL
C                      TRIANGLE WITH PERIMETER GIVEN IN X AND Y
C
C***********************************************************************
C
      PARAMETER  (RLARGE = 100000.)
      DIMENSION X (N), Y (N), ANGLE (N)
      DIMENSION SMANG (7), INDEX (7)
      DIMENSION ISORT (3)
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
C  DETERMINE OPTIMUM ORIGIN / SHAPE COMBINATION FOR A TRIANGLE
C
      ATOL = PI * 150. / 180.
C
C  FIND SIDE DIVISION USING 5 SMALLEST ANGLES AND CHECK CONDITION
C
      DO 140 I = 1,  3
         ISORT (I) = INDEX (I)
  140 CONTINUE
      DO 160 I = 1,  2
         DO 150 J = I + 1,  3
            IF (ISORT (I) .GT. ISORT (J)) THEN
               ITMP = ISORT (I)
               ISORT (I) = ISORT (J)
               ISORT (J) = ITMP
            ENDIF
  150    CONTINUE
  160 CONTINUE
      I1 = ISORT (1)
      I2 = ISORT (2)
      I3 = ISORT (3)
      MM1 = I2  -  I1
      IF (MM1 .LT. 0) MM1 = N  +  MM1
      MM2 = I3  -  I2
      IF (MM2 .LT. 0) MM2 = N  +  MM2
      MM3 = N  -  MM1  -  MM2
      MAX = MAX0 (MM1,  MM2,  MM3)
      IF (MAX .LE. N - MAX - 2) THEN
C
C  ADD UP ASSIGNED ANGLES
C
         IFIRST = I1
         M1 = MM1
         M2 = MM2
         GBEST = ANGLE (I1) + ANGLE (I2) + ANGLE (I3)
      ELSE
         IFIRST = 1
         GBEST = RLARGE
      END IF
C
C  LIMIT THE SIZE OF ANY SIDE
C
      MMAX = (N - 2) / 2
C
C  GO AROUND THE PERIMETER USING THE 10 SMALLEST ANGLES AS POSSIBLE
C  STARTING POINTS,  AND THEN FIND THE BEST COMBINATION OF SIDE LENGTHS
C
      DO 190 ISA = 1, NSA
         IF (SMANG (ISA) .LE. ATOL) THEN
            I1 = INDEX (ISA)
            SUM1 = ANGLE (I1)
            IF (SUM1  .GE.  GBEST) GO TO 190
C
C  ASSIGN A TRIAL SECOND NODE
C
            DO 180 N1 = 2, MMAX
               I2 = I1 + N1
               IF (I2 .GT. N)I2 = I2 - N
               SUM2 = SUM1 + ANGLE (I2)
               IF (SUM2  .GE.  GBEST) GO TO 180
C
C  ASSIGN A TRIAL THIRD NODE
C
               DO 170 N2 = 2, N - N1 - 2
                  I3 = I2 + N2
                  IF (I3 .GT. N)I3 = I3 - N
                  GVAL = SUM2 + ANGLE (I3)
                  IF (GVAL  .GE.  GBEST) GO TO 170
C
C  FIND SIDE DIVISION AND CHECK CONDITION
C
                  MM1 = I2  -  I1
                  IF  (MM1 .LT. 0) MM1 = N + MM1
                  MM2 = I3  -  I2
                  IF  (MM2 .LT. 0) MM2 = N + MM2
C ... Guess by GDS, MM1 substituted for MM?
                  MM3 = N - MM1  - MM2
                  MAX = MAX0 (MM1, MM2, MM3)
                  IF  (MAX .LE. N - MAX - 2) THEN
C
C  ADD UP ASSIGNED ANGLES AND COMPARE TO PREVIOUS TRIALS
C
                     IFIRST = I1
                     M1 = MM1
                     M2 = MM2
                     GBEST = GVAL
                  ENDIF
  170          CONTINUE
  180       CONTINUE
         ENDIF
  190 CONTINUE
C
      RETURN
C
      END
