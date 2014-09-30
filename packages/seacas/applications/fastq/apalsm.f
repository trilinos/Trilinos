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

C $Id: apalsm.f,v 1.1 1990/11/30 11:03:37 gdsjaar Exp $
C $Log: apalsm.f,v $
C Revision 1.1  1990/11/30 11:03:37  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]APALSM.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE APALSM (MXND, XN, YN, LXK, KXL, NXL, LXN, NNN, NNNOLD,
     &   NIT, TOL, RO, ALPHA, ERR)
C***********************************************************************
C
C  SUBROUTINE APALSM = AREA PULL AND LAPLACIAN MESH SMOOTHER
C
C***********************************************************************
C
C  NOTE:
C     IN THIS SMOOTHER EACH NODE IS SUCCESSIVELY MOVED BY
C     AN AMOUNT GIVEN AS A WEIGHTED AVERAGE OF AN *AREA PULL*
C     VECTOR AND THE LAPLACIAN VECTOR (AVERAGE OF VECTORS POINTING
C     TO NEIGHBORS).  THE *AREA PULL* VECTOR IS OBTAINED BY LETTING
C     EACH ELEMENT PULL IN PERPENDICULARLY ON ITS SIDES WITH FORCE
C     PROPORTIONAL TO THE LENGTH OF THAT SIDE TIMES A FACTOR
C     INVOLVING THE AREAS OF THIS ELEMENT AND ITS NEIGHBOR SHARING
C     THAT SIDE.
C
C***********************************************************************
C
C  VARIABLES USED:
C     NIT   = MAX ITERATIONS TO DO
C     TOL   = PERCENT OF SMALLEST CONNECTING LINE TO USE AS NODE MOVEMENT
C             CONVERGENCE TOLERANCE.
C     RO    = UNDER OR OVER-RELAXATION FACTOR.
C     ALPHA = WEIGHT GIVEN TO AREA PULL VECTOR.  USUALLY = 0.5.
C             WEIGHT GIVEN TO LAPLACIAN VECTOR = 1.-ALPHA.
C
C***********************************************************************
C
      DIMENSION LXK(4, MXND), KXL(2, 3*MXND), NXL(2, 3*MXND)
      DIMENSION LXN(4, MXND), XN(MXND), YN(MXND)
      DIMENSION LINES(20), NS1(4), NS2(4)
      DIMENSION KLIB(8), NLIB(4, 8), ALIB(8), XCLIB(8), YCLIB(8)
C
      LOGICAL BIG, ERR
C
      ERR = .FALSE.
      TOL2 = TOL**2
      BETA = 1. - ALPHA
C
C  ITERATION LOOP
C
      DO 160 IT = 1, NIT
         BIG = .FALSE.
C
C  NODE LOOP
C
         NNT = 0
         DO 150 NODE = NNNOLD  +  1, NNN
C
C  CHECK FOR CONTINUATIONS,  BOUNDARY,  OR RELAXED NODE
C
            IF ((LXN(3, NODE) .GE.  0) .AND. (LXN(2, NODE) .GT. 0)
     &         .AND. (LXN(1, NODE) .GT. 0)) THEN
               NNT = NNT + 1
C
C  INITIALIZE
C
               KNUM = 0
               XA = 0.
               YA = 0.
               XL = 0.
               YL = 0.
C
C  PROCESS EACH LINE CONNECTED TO THIS NODE
C
               CALL GETLXN (MXND, LXN, NODE, LINES, KOUNT, ERR)
               IF (ERR) RETURN
               DO 100 IL = 1, KOUNT
                  L = LINES(IL)
                  N1 = NXL(1, L)
                  N2 = NXL(2, L)
C
C  FETCH ELEMENT DATA
C
                  IF (KXL(1, L) .GT. 0) CALL APALIB (MXND, XN, YN, LXK,
     &               NXL, KXL(1, L), NS1, AREA1, XCEN1, YCEN1, KNUM,
     &               KLIB, NLIB, ALIB, XCLIB, YCLIB)
                  IF (KXL(2, L) .GT. 0) CALL APALIB (MXND, XN, YN, LXK,
     &               NXL, KXL(2, L), NS2, AREA2, XCEN2, YCEN2, KNUM,
     &               KLIB, NLIB, ALIB, XCLIB, YCLIB)
C
C  GET FORCE DIRECTION MODULO PI RADIANS.
C  CORRECT FOR WRONG DIRECTION BY ALIGNING WITH THE VECTOR
C  FROM (XCEN1, YCEN1) TO (XCEN2, YCEN2).
C
                  XF = -(YN(N2) - YN(N1))
                  YF = XN(N2) - XN(N1)
                  DOT = XF*(XCEN2 - XCEN1) + YF*(YCEN2 - YCEN1)
                  IF (DOT  .LT.  0.) THEN
                     XF = -XF
                     YF = -YF
                  END IF
C
C  TAKE CARE OF ZERO AREAS
C
                  IF ((AREA1  .LE.  0) .OR. (AREA2  .LE.  0)) THEN
                     AREA1 = 1.0
                     AREA2 = 1.0
                  END IF
C
C  UPDATE AREA PULL VECTOR SUM
C
                  WEIGHT = (AREA2 - AREA1)/(AREA2 + AREA1)
                  XA = XA  +  WEIGHT*XF
                  YA = YA  +  WEIGHT*YF
C
C  UPDATE LAPLACIAN VECTOR SUM
C
                  NOE = N1 + N2 - NODE
                  DX = XN(NOE) - XN(NODE)
                  DY = YN(NOE) - YN(NODE)
                  XL = XL + DX
                  YL = YL + DY
C
C  UPDATE LEAST LENGTH
C
                  DIST2 = DX*DX  +  DY*DY
                  IF (IL .EQ. 1) DMIN2 = DIST2
                  DMIN2 = MIN(DMIN2, DIST2)
  100          CONTINUE
C
C  COMPUTE NET MOVEMENT VECTOR
C
               RK = 1.0/FLOAT(KOUNT)
               XNET = (ALPHA*XA  +  BETA*XL)*RK
               YNET = (ALPHA*YA  +  BETA*YL)*RK
C
C  MOVE THE NODE
C
               YN(NODE) = YN(NODE)  +  YNET * RO
               XN(NODE) = XN(NODE)  +  XNET * RO
C
C  CHECK FOR SIGNIFICANT MOVEMENT
C
               IF (XNET*XNET + YNET*YNET .GT. TOL2*DMIN2) THEN
C
C  SIGNIFICANT MOVEMENT - REMOVE RELAXATION FLAGS
C
C  FIRST FROM DIRECTLY CONNECTED NODES
C
                  DO 110 IL = 1, KOUNT
                     L = LINES(IL)
                     NOE = NXL(1, L) + NXL(2, L) - NODE
                     LXN(3, NOE) = ABS(LXN(3, NOE))
  110             CONTINUE
C
C  NEXT FROM DIAGONALLY OPPOSITE NODES (MAX 8)
C
                  DO 140 IK = 1, KNUM
                     DO 120 I = 1, 4
                        IF (NLIB(I, IK) .EQ. NODE) THEN
                           IDIAG = I + 2
                           IF (IDIAG .GE.  5) IDIAG = IDIAG - 4
                           NDIAG = NLIB(IDIAG, IK)
                           LXN(3, NDIAG) = ABS(LXN(3, NDIAG))
                           GO TO 130
                        END IF
  120                CONTINUE
                     CALL MESAGE ('ERROR IN APALSM')
                     ERR = .TRUE.
                     RETURN
  130                CONTINUE
  140             CONTINUE
C
C  INSIGNIFICANT MOVEMENT
C  INSERT RELAXATION FLAG
C
               ELSE
                  LXN(3, NODE) = -ABS(LXN(3, NODE))
               END IF
C
            END IF
  150    CONTINUE
         IF (.NOT.BIG) GO TO 170
  160 CONTINUE
      IT = NIT
C
C  REMOVE ALL FLAGS
C
  170 CONTINUE
      DO 180 NODE = NNNOLD  +  1, NNN
         LXN(3, NODE) = ABS(LXN(3, NODE))
  180 CONTINUE
C
      RETURN
      END
