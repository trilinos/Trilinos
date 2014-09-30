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

C $Id: arelax.f,v 1.2 1999/06/17 19:16:50 gdsjaar Exp $
      SUBROUTINE ARELAX (MXND, XN, YN, LXK, KXL, NXL, LLL, ARFACT)
C***********************************************************************
C
C  SUBROUTINE ARELAX = CALCULATES UNDER - RELAXATION FACTOR FOR AREA PULL
C                      AND LAPLACIAN SMOOTHER
C
C***********************************************************************
C
C  NOTE:
C     THE AREA PULL AND LAPLACIAN SMOOTHER WILL OVER - CORRECT
C     AND BECOME UNSTABLE WHEN TYPICAL MESH ELEMENTS ARE MUCH
C     LONGER THAN THEY ARE WIDE, SAY BY A FACTOR OF SIX OR MORE.
C     THIS ROUTINE COMPUTES AN APPROPRIATE UNDER - RELAXATION
C     FACTOR TO BE USED TO HELP CORRECT THIS PROBLEM.  ON REGIONS
C     WHICH HAVE GENERALLY NEAR SQUARE ELEMENTS WITH A SMALL
C     PERCENTAGE OF VERY LONG THIN ELEMENTS THIS FACTOR WILL
C     PROBABLY NOT ADEQUATELY HANDLE THE DIFFICULTY.  IN SUCH
C     SITUATIONS AN ALTERNATE SMOOTHER  (SUCH AS THE CENTROID - AREA -
C     PULL) SHOULD BE USED.
C     THE FACTOR RETURNED BY THIS ROUTINE MAY BE LARGER THAN ONE,
C     WHICH MEANS THAT OVER - RELAXATION IS APPROPRIATE.
C
C***********************************************************************
C
      DIMENSION NODES (4), LXK (4, MXND), KXL (2, 3 * MXND)
      DIMENSION NXL (2, 3 * MXND)
      DIMENSION XN (MXND), YN (MXND)
C
      LOGICAL CCW
C
      ARFACT = 1.0
      RATSUM = 0.
      NUM = 0
C
      DO 100 MYL = 1, LLL
C
C  SKIP BOUNDARY LINES
C
         IF (KXL(1, myL) .gt. 0 .and. KXL (2, myL) .GT. 0) THEN
            CCW = .TRUE.
            CALL GNXKA (MXND, XN, YN, KXL (1, MYL), NODES, AREA1, LXK,
     &         NXL, CCW)
            CALL GNXKA (MXND, XN, YN, KXL (2, MYL), NODES, AREA2, LXK,
     &         NXL, CCW)
            N1 = NXL (1, MYL)
            N2 = NXL (2, MYL)
            DXDY = (XN (N2) - XN (N1)) **2 + (YN (N2) - YN (N1)) **2
            IF (AREA1 + AREA2 .GT. 0) THEN
               RATIO = 2.0 * DXDY /  (AREA1 + AREA2)
               IF (RATIO .GE. 0.99) THEN
                  NUM = NUM + 1
                  RATSUM = RATSUM + RATIO
               ENDIF
            ENDIF
         ENDIF
  100 CONTINUE
C
      IF (NUM .LE. 0) RETURN
      ASPECT = RATSUM / FLOAT (NUM)
      ARFACT = AMIN1 (2.0 / ASPECT, 1.5)
C
      RETURN
C
      END
