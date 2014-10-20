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

C $Id: isolap.f,v 1.1 1990/11/30 11:10:35 gdsjaar Exp $
C $Log: isolap.f,v $
C Revision 1.1  1990/11/30 11:10:35  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]ISOLAP.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE ISOLAP (MXND, XN, YN, LXK, KXL, NXL, LXN, NNN, NNNOLD,
     &   WFAC, NIT, EPS, RO)
C***********************************************************************
C
C  SUBROUTINE ISOLAP = MESH SMOOTHING BY LAPLACE-S USING GAUSS-SEIDEL
C
C***********************************************************************
C
C  VARIABLES USED:
C     WFAC = WEIGTH (0. = LAPLACIAN, 1. = ISOPARAMETRIC)
C     NIT  = THE MAX NUMBER OF ITERATIONS TO DO.
C     EPS  = MINIMUM DISTANCE NODES MUST MOVE TO CONTINUE ITERATIONS
C     RO   = AN UNDER- OR OVER-RELAXATION FACTOR (NORMALLY 1.0)
C
C***********************************************************************
C
      DIMENSION KLIST(20), NODES(4)
      DIMENSION XN(MXND), YN(MXND)
      DIMENSION LXN(4, MXND), NXL(2, 3*MXND)
      DIMENSION LXK(4, MXND), KXL(2, 3*MXND)
C
      LOGICAL BIG, CCW, ERR
C
      IF (RO .LT. 0.01) RO = 1.
      EPS2 = (EPS*RO)**2
C
C  ITERATION LOOP
C
      DO 140 IT = 1, NIT
         BIG = .FALSE.
C
C  NODE LOOP
C
         DO 130 NODE = NNNOLD + 1, NNN
C
C  SKIP CONTINUATION AND BOUNDARY LINES
C
            IF ((LXN(1, NODE) .GT. 0) .AND. (LXN(2, NODE) .GT. 0)) THEN
C
C  FIND ELEMENTS ATTACHED TO NODE
C
               CALL GKXN (MXND, KXL, LXN, NODE, KS, KLIST, ERR)
C
               SUMX = 0.0
               SUMY = 0.0
               DO 120 KL = 1, KS
                  CCW = .FALSE.
                  KK = KLIST(KL)
                  CALL GNXKA (MXND, XN, YN, KK, NODES, AREA, LXK, NXL,
     &               CCW)
C
                  DO 100 IN = 1, 4
                     IF (NODES(IN) .EQ. NODE) THEN
                        J1 = IN + 1
                        IF (J1 .GT. 4) J1 = 1
                        GO TO 110
                     END IF
  100             CONTINUE
  110             CONTINUE
                  J2 = J1 + 1
                  IF (J2 .GT. 4) J2 = 1
                  J3 = J2 + 1
                  IF (J3 .GT. 4) J3 = 1
C
                  SUMX = SUMX + XN(NODES(J1)) + XN(NODES(J3))
     &               - WFAC*XN(NODES(J2))
                  SUMY = SUMY + YN(NODES(J1)) + YN(NODES(J3))
     &               - WFAC*YN(NODES(J2))
  120          CONTINUE
C
C  REDEFINE THIS NODE-S COORDINATES
C
               SUMX = SUMX/(FLOAT(KS)*(2.0 - WFAC))
               SUMY = SUMY/(FLOAT(KS)*(2.0 - WFAC))
               XDEL = RO*(SUMX-XN(NODE))
               YDEL = RO*(SUMY-YN(NODE))
               XN(NODE) = XN(NODE) + XDEL
               YN(NODE) = YN(NODE) + YDEL
C
C  CHECK FOR CONVERGENCE
C
               IF ((XDEL*XDEL + YDEL*YDEL) .GT. EPS2) BIG = .TRUE.
            ENDIF
  130    CONTINUE
C
C  IF NO SIGNIFICANT MOVEMENTS OCCURRED,  RETURN
C
         IF (.NOT.BIG) RETURN
  140 CONTINUE
      RETURN
      END
