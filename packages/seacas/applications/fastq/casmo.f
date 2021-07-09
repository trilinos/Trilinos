C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE CASMO (MXND, XN, YN, LXK, KXL, NXL, LXN, NNN, NNNOLD,
     &   NIT, EPS, RO)
C***********************************************************************

C  SUBROUTIINE CASMO  =  CENTROID-AREA-PULL METHOD MESH SMOOTHING

C***********************************************************************

C  NOTE:
C     IN THIS METHOD EACH NODE IS PULLED TOWARD THE CENTROIDS OF
C     ADJACENT ELEMENTS BY FORCES PROPORTIONAL TO THE RESPECTIVE
C     ELEMENT AREAS.
C     IDEA BY STEVE PETTY AND RONDALL JONES

C***********************************************************************

C  VARIABLES USED:
C     NIT  =  MAX NUMBER OF ITERATIONS TO DO
C     EPS  =  NODE MOVEMENT TOLERANCE FOR CONVERGENCE
C     RO   =  UNDER OR OVER-RELAXATION FACTOR.

C***********************************************************************

      DIMENSION NODES(4)
      DIMENSION KLIST(20), AREA(20), XCEN(20), YCEN(20)
      DIMENSION LXK(4, MXND), KXL(2, 3*MXND), NXL(2, 3*MXND)
      DIMENSION LXN(4, MXND), XN(MXND), YN(MXND)

      LOGICAL ERR, BIG, CCW

C  ITERATION LOOP

      DO 140 IT = 1, NIT
         BIG = .FALSE.

C  NODE LOOP

         DO 130 NODE = NNNOLD + 1, NNN

C  SKIP CONTINUATIONS AND BOUNDARY NODES

            IF((LXN(1, NODE).GT.0).AND.(LXN(2, NODE).GT.0))THEN

C  GET ELEMENT LIST (IGNORE ERR IF CAUSED BY TOO MANY ELEMENTS)

               CALL GKXN (MXND, KXL, LXN, NODE, NUMK, KLIST, ERR)
               IF((ERR).AND.(NUMK.LT.20))RETURN

C  GET AREAS AND CENTROIDS

               DO 110 IK = 1, NUMK
                  KK = KLIST(IK)
                  CCW = .TRUE.
                  CALL GNXKA (MXND, XN, YN, KK, NODES, AREA(IK), LXK,
     &               NXL, CCW)
                  XSUM = 0.
                  YSUM = 0.
                  DO 100 I = 1, 4
                     NN = NODES(I)
                     XSUM = XSUM + XN(NN)
                     YSUM = YSUM + YN(NN)
  100             CONTINUE
                  XCEN(IK) = 0.25*XSUM
                  YCEN(IK) = 0.25*YSUM
  110          CONTINUE

C  COMPUTE AND SUM THE FORCE VECTORS

               FX = 0.
               FY = 0.
               SUMW = 0.
               DO 120 IK = 1, NUMK
                  DX = XCEN(IK) - XN(NODE)
                  DY = YCEN(IK) - YN(NODE)
                  WEIGHT = AREA(IK)
                  SUMW = SUMW + WEIGHT
                  FX = FX + WEIGHT*DX
                  FY = FY + WEIGHT*DY
  120          CONTINUE

C  NORMALIZE THE RESULTANT VECTOR

               RSUMW = 1.0/SUMW
               FX = FX*RSUMW
               FY = FY*RSUMW

C  MOVE THE NODE

               DELX = RO*FX
               DELY = RO*FY
               XN(NODE) = XN(NODE) + DELX
               YN(NODE) = YN(NODE) + DELY
            ENDIF
  130    CONTINUE
         IF(.NOT.BIG)RETURN
  140 CONTINUE
      RETURN
      END
