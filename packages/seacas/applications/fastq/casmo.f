C $Id: casmo.f,v 1.1 1990/11/30 11:04:15 gdsjaar Exp $
C $Log: casmo.f,v $
C Revision 1.1  1990/11/30 11:04:15  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]CASMO.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE CASMO (MXND, XN, YN, LXK, KXL, NXL, LXN, NNN, NNNOLD,
     &   NIT, EPS, RO)
C***********************************************************************
C
C  SUBROUTIINE CASMO  =  CENTROID-AREA-PULL METHOD MESH SMOOTHING
C
C***********************************************************************
C
C  NOTE:
C     IN THIS METHOD EACH NODE IS PULLED TOWARD THE CENTROIDS OF
C     ADJACENT ELEMENTS BY FORCES PROPORTIONAL TO THE RESPECTIVE
C     ELEMENT AREAS.
C     IDEA BY STEVE PETTY AND RONDALL JONES
C
C***********************************************************************
C
C  VARIABLES USED:
C     NIT  =  MAX NUMBER OF ITERATIONS TO DO
C     EPS  =  NODE MOVEMENT TOLERANCE FOR CONVERGENCE
C     RO   =  UNDER OR OVER-RELAXATION FACTOR.
C
C***********************************************************************
C
      DIMENSION NODES(4)
      DIMENSION KLIST(20), AREA(20), XCEN(20), YCEN(20)
      DIMENSION LXK(4, MXND), KXL(2, 3*MXND), NXL(2, 3*MXND)
      DIMENSION LXN(4, MXND), XN(MXND), YN(MXND)
C
      LOGICAL ERR, BIG, CCW
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
C  SKIP CONTINUATIONS AND BOUNDARY NODES
C
            IF((LXN(1, NODE).GT.0).AND.(LXN(2, NODE).GT.0))THEN
C
C  GET ELEMENT LIST (IGNORE ERR IF CAUSED BY TOO MANY ELEMENTS)
C
               CALL GKXN (MXND, KXL, LXN, NODE, NUMK, KLIST, ERR)
               IF((ERR).AND.(NUMK.LT.20))RETURN
C
C  GET AREAS AND CENTROIDS
C
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
C
C  COMPUTE AND SUM THE FORCE VECTORS
C
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
C
C  NORMALIZE THE RESULTANT VECTOR
C
               RSUMW = 1.0/SUMW
               FX = FX*RSUMW
               FY = FY*RSUMW
C
C  MOVE THE NODE
C
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
