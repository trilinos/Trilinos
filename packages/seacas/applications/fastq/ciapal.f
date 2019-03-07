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

C $Id: ciapal.f,v 1.1 1990/11/30 11:04:43 gdsjaar Exp $
C $Log: ciapal.f,v $
C Revision 1.1  1990/11/30 11:04:43  gdsjaar
C Initial revision
C
CC* FILE: [.QMESH]CIAPAL.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE CIAPAL (MXND, XN, YN, LXK, KXL, NXL, LXN, NNN, NNNOLD,
     &   NIT, EPS, RO, ALPHA)
C***********************************************************************
C
C  SUBROUTINE CIAPAL = CENTROID INVERSE AREA PUSH AND LAPLACIAN SMOOTH
C
C***********************************************************************
C
C  NOTE:
C     IN THIS METHOD EACH CENTROID OF AN ELEMENT PUSHES OUT
C     ON THE SURROUNDING NODES WITH A FORCE INVERSELY PROPORTIONAL
C     TO THE AREA OF THE ELEMENT WHILE IT SIMULTANEOUSLY PULLS ON
C     EACH NODE WITH A FORCE PROPORTIONAL TO THE LENGTH OF THE LINE
C     CONNECTING THE CENTROID WITH EACH NODE.
C
C***********************************************************************
C
C  VARIABLES USED:
C     NIT   = MAX NUMBER OF ITERATIONS TO DO
C     EPS   = NODE MOVEMENT TOLERANCE FOR CONVERGENCE
C     RO    = UNDER OR OVER-RELAXATION FACTOR.
C     ALPHA = WEIGHT GIVEN TO AREA-PUSH VECTOR.  USUALLY=0.5.
C             WEIGHT GIVEN TO LAPLACIAN VECTOR = 1.-ALPHA.
C
C***********************************************************************
C
      DIMENSION NODES (4)
      DIMENSION KLIST (20), AREA (20), XCEN (20), YCEN (20)
      DIMENSION LXK (4, MXND), KXL (2, 3*MXND)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)
      DIMENSION XN (MXND), YN (MXND)
C
      LOGICAL BIG, CCW, ERR
      EPS2 =  (EPS * RO) **2
      BETA = 1.0 - ALPHA
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
C  SKIP CONTINUATION LINES AND BOUNDARY LINES
C
            IF ((LXN (1, NODE).GT.0) .AND. (LXN (2, NODE).GT.0)) THEN
C
C  GET ELEMENT LIST  (IGNORE ERR IF IT IS BECAUSE TOO MANY WERE FOUND)
C
               CALL GKXN (MXND, KXL, LXN, NODE, NUMK, KLIST, ERR)
               IF ((ERR) .AND. (NUMK.LT.20))RETURN
C
C  GET AREAS AND CENTROIDS
C
               DO 110 IK = 1, NUMK
                  KK = KLIST (IK)
                  CCW = .TRUE.
                  CALL GNXKA (MXND, XN, YN, KK, NODES, AREA (IK), LXK,
     &               NXL, CCW)
                  XSUM = 0.
                  YSUM = 0.
                  DO 100 I = 1, 4
                     NN = NODES (I)
                     XSUM = XSUM + XN (NN)
                     YSUM = YSUM + YN (NN)
  100             CONTINUE
                  XCEN (IK) = 0.25 * XSUM
                  YCEN (IK) = 0.25 * YSUM
  110          CONTINUE
C
C  COMPUTE AND SUM THE FORCE VECTORS
C
               FX = 0.
               FY = 0.
               SUMW = 0.
               SDX = 0.
               SDY = 0.
               DO 120 IK = 1, NUMK
                  DX = XCEN (IK) - XN (NODE)
                  DY = YCEN (IK) - YN (NODE)
                  AL2 = DX * DX + DY * DY
                  ARL = AREA (IK) * SQRT (AL2)
                  WEIGHT = 1.0E20
                  IF (ARL.GT.0.)WEIGHT = 1.0/ARL
                  SUMW = SUMW + WEIGHT
                  FX = FX - WEIGHT * DX
                  FY = FY - WEIGHT * DY
                  SDX = SDX + DX
                  SDY = SDY + DY
  120          CONTINUE
C
C  MOVE THE NODE
C
               RSUMW = 1.0/SUMW
               RNUMK = 1.0/FLOAT (NUMK)
               FX = ALPHA * FX * RSUMW + BETA * SDX * RNUMK
               FY = ALPHA * FY * RSUMW + BETA * SDY * RNUMK
               DELX = RO * FX
               DELY = RO * FY
               XN (NODE) = XN (NODE) + DELX
               YN (NODE) = YN (NODE) + DELY
               IF (DELX ** 2 + DELY ** 2 .GT. EPS2)BIG = .TRUE.
            ENDIF
  130    CONTINUE
         IF (.NOT.BIG) RETURN
  140 CONTINUE
      RETURN
      END
