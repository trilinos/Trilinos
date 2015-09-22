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

C $Id: arcy.f,v 1.2 1998/07/14 18:18:21 gdsjaar Exp $
C $Log: arcy.f,v $
C Revision 1.2  1998/07/14 18:18:21  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.1.1.1  1990/11/30 11:03:44  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:03:43  gdsjaar
c Initial revision
c 
C
CC* FILE: [.MAIN]ARCY.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE ARCY (XCEN, YCEN, THETA1, THETA2, XK, XA, X, YTRY, ERR)
C***********************************************************************
C
C  SUBROUTINE ARCY = ITERATIVELY SOLVES THE LOGARITHMIC SPIRAL PROBLEM
C                    TO DETERMINE A Y VALUE GIVEN AN X THAT INTERSECTS
C                    THE ARC
C
C***********************************************************************
C
      LOGICAL ERR
C
C  START WITH 10 INCREMENTS, EACH PASS INCREMENTS DECREASE TEN FOLD
C
      ANGINC = (THETA2 - THETA1) * .05
      ANG = THETA1
      F1 = SPIRAL (XA, XK, X, XCEN, YCEN, ANG)
      IF (F1 .EQ. 0.0) THEN
         THETA = ANG
         GO TO 110
      END IF
      ANG2 = ANG + ANGINC
  100 CONTINUE
      F2 = SPIRAL (XA, XK, X, XCEN, YCEN, ANG2)
      IF (F2 .EQ. 0.0) THEN
         THETA = ANG2
      ELSE IF (F1*F2 .LT. 0.0) THEN
         THETA = SOLVE(XA, XK, X, XCEN, YCEN, ANG, ANG2)
      ELSE
         ANG = ANG2
         ANG2 = ANG2 + ANGINC
         IF (ANG2 .LE. THETA2) GO TO 100
         ERR = .TRUE.
         GO TO 120
      END IF
C
  110 CONTINUE
      YTRY = (XA * EXP(XK * THETA)) * SIN(THETA) + YCEN
C
  120 CONTINUE
C
C  FIND THE SECOND ROOT IF THE FIRST ONE HAS BEEN LOCATED
C
      IF(.NOT.ERR)THEN
         ANG=THETA+ANGINC
         F1 = SPIRAL (XA, XK, X, XCEN, YCEN, ANG)
         F2 = SPIRAL (XA, XK, X, XCEN, YCEN, THETA2)
         IF (F1 .EQ. 0.0) THEN
            THETA = ANG
         ELSEIF (F2 .EQ. 0.0) THEN
            THETA = THETA2
         ELSE IF (F1*F2 .LT. 0.0) THEN
            THETA = SOLVE(XA, XK, X, XCEN, YCEN, ANG, THETA2)
         ELSE
            GO TO 130
         END IF
      END IF
C
      YTRY2 = (XA * EXP(XK * THETA)) * SIN(THETA) + YCEN
      YTRY = MAX(YTRY,YTRY2)
  130 CONTINUE
      RETURN
C
      END
