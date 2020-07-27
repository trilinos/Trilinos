C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE CIAPAL (MXND, XN, YN, LXK, KXL, NXL, LXN, NNN, NNNOLD,
     &   NIT, EPS, RO, ALPHA)
C***********************************************************************

C  SUBROUTINE CIAPAL = CENTROID INVERSE AREA PUSH AND LAPLACIAN SMOOTH

C***********************************************************************

C  NOTE:
C     IN THIS METHOD EACH CENTROID OF AN ELEMENT PUSHES OUT
C     ON THE SURROUNDING NODES WITH A FORCE INVERSELY PROPORTIONAL
C     TO THE AREA OF THE ELEMENT WHILE IT SIMULTANEOUSLY PULLS ON
C     EACH NODE WITH A FORCE PROPORTIONAL TO THE LENGTH OF THE LINE
C     CONNECTING THE CENTROID WITH EACH NODE.

C***********************************************************************

C  VARIABLES USED:
C     NIT   = MAX NUMBER OF ITERATIONS TO DO
C     EPS   = NODE MOVEMENT TOLERANCE FOR CONVERGENCE
C     RO    = UNDER OR OVER-RELAXATION FACTOR.
C     ALPHA = WEIGHT GIVEN TO AREA-PUSH VECTOR.  USUALLY=0.5.
C             WEIGHT GIVEN TO LAPLACIAN VECTOR = 1.-ALPHA.

C***********************************************************************

      DIMENSION NODES (4)
      DIMENSION KLIST (20), AREA (20), XCEN (20), YCEN (20)
      DIMENSION LXK (4, MXND), KXL (2, 3*MXND)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)
      DIMENSION XN (MXND), YN (MXND)

      LOGICAL BIG, CCW, ERR
      EPS2 =  (EPS * RO) **2
      BETA = 1.0 - ALPHA

C  ITERATION LOOP

      DO 140 IT = 1, NIT
         BIG = .FALSE.

C  NODE LOOP

         DO 130 NODE = NNNOLD + 1, NNN

C  SKIP CONTINUATION LINES AND BOUNDARY LINES

            IF ((LXN (1, NODE).GT.0) .AND. (LXN (2, NODE).GT.0)) THEN

C  GET ELEMENT LIST  (IGNORE ERR IF IT IS BECAUSE TOO MANY WERE FOUND)

               CALL GKXN (MXND, KXL, LXN, NODE, NUMK, KLIST, ERR)
               IF ((ERR) .AND. (NUMK.LT.20))RETURN

C  GET AREAS AND CENTROIDS

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

C  COMPUTE AND SUM THE FORCE VECTORS

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

C  MOVE THE NODE

               RSUMW = 1.0/SUMW
               RNUMK = 1.0/DBLE(NUMK)
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
