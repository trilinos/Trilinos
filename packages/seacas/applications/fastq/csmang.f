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

C $Id: csmang.f,v 1.2 1991/03/21 15:44:31 gdsjaar Exp $
C $Log: csmang.f,v $
C Revision 1.2  1991/03/21 15:44:31  gdsjaar
C Changed all 3.14159... to atan2(0.0, -1.0)
C
c Revision 1.1.1.1  1990/11/30  11:05:35  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:05:34  gdsjaar
c Initial revision
c 
C
CC* FILE: [.QMESH]CSMANG.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE CSMANG (N, X, Y, ANGLE, NSANG, SMANG, INDEX)
C***********************************************************************
C
C  SUBROUTINE CSMANG = CALCULATES THE "NSANG" SMALLEST ANGLES
C                      AND PLACES THEM IN THE SMANG ARRAY WITH
C                      THE INDICES BEING PLACED IN THE INDEX ARRAY
C
C***********************************************************************
C
C     OBSERVATION: - IT DOES NOT MATTER THE ANGLE ORIENTATION.
C                  - THE ANGLES ARE STORE IN THE ANGLE ARRAY AS THEY
C                    APPEAR.
C                  - THE SMALLEST ANGLES ARE IN ASCENDING ORDER.
C                - THE INDEX ARRAY RETURNS THE SMALLEST ANGLE POSITION
C                    IN ASCENDING ORDER.
C
C
C     MODIFIED BY : HORACIO RECALDE             DATE:JAN 1988
C***********************************************************************
C
      DIMENSION X(N), Y(N), ANGLE(N)
      DIMENSION SMANG (NSANG + 1), INDEX (NSANG + 1)

      PI = ATAN2(0.0, -1.0)
      TWOPI = 2.0 * PI
C
C  FORM THE LIST OF SMALLEST ANGLES
C
      NSA = NSANG
      DO 100 I = 1,NSA
         SMANG(I) = 10.
         INDEX(I) = 0
  100 CONTINUE
C
      AGOLD = ATAN2 (Y (1) - Y(N), X (1) - X (N))
C
      DO 130 J = 1, N
C
C  GET THE ANGLE FORMED BY THIS SET OF POINTS
C
         NEXT = J + 1
         IF (NEXT .GT. N) NEXT = 1
         AGNEW = ATAN2 (Y (NEXT) - Y (J), X (NEXT) - X (J))
         DIFF = AGNEW - AGOLD
         IF (DIFF .GT. PI)  DIFF = DIFF - TWOPI
         IF (DIFF .LT. -PI) DIFF = DIFF + TWOPI
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
            IF (SMANG (I + 1) .GE. SMANG (I)) GO TO 120
            TEMP = SMANG(I)
            ITEMP = INDEX(I)
            SMANG (I) = SMANG (I + 1)
            INDEX (I) = INDEX (I + 1)
            SMANG (I + 1) = TEMP
  110       INDEX (I + 1) = ITEMP
  120    CONTINUE
C
  130 CONTINUE
C
      RETURN
C
      END
