C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE APALSM (MXND, XN, YN, LXK, KXL, NXL, LXN, NNN, NNNOLD,
     &   NIT, TOL, RO, ALPHA, ERR)
C***********************************************************************

C  SUBROUTINE APALSM = AREA PULL AND LAPLACIAN MESH SMOOTHER

C***********************************************************************

C  NOTE:
C     IN THIS SMOOTHER EACH NODE IS SUCCESSIVELY MOVED BY
C     AN AMOUNT GIVEN AS A WEIGHTED AVERAGE OF AN *AREA PULL*
C     VECTOR AND THE LAPLACIAN VECTOR (AVERAGE OF VECTORS POINTING
C     TO NEIGHBORS).  THE *AREA PULL* VECTOR IS OBTAINED BY LETTING
C     EACH ELEMENT PULL IN PERPENDICULARLY ON ITS SIDES WITH FORCE
C     PROPORTIONAL TO THE LENGTH OF THAT SIDE TIMES A FACTOR
C     INVOLVING THE AREAS OF THIS ELEMENT AND ITS NEIGHBOR SHARING
C     THAT SIDE.

C***********************************************************************

C  VARIABLES USED:
C     NIT   = MAX ITERATIONS TO DO
C     TOL   = PERCENT OF SMALLEST CONNECTING LINE TO USE AS NODE MOVEMENT
C             CONVERGENCE TOLERANCE.
C     RO    = UNDER OR OVER-RELAXATION FACTOR.
C     ALPHA = WEIGHT GIVEN TO AREA PULL VECTOR.  USUALLY = 0.5.
C             WEIGHT GIVEN TO LAPLACIAN VECTOR = 1.-ALPHA.

C***********************************************************************

      DIMENSION LXK(4, MXND), KXL(2, 3*MXND), NXL(2, 3*MXND)
      DIMENSION LXN(4, MXND), XN(MXND), YN(MXND)
      DIMENSION LINES(20), NS1(4), NS2(4)
      DIMENSION KLIB(8), NLIB(4, 8), ALIB(8), XCLIB(8), YCLIB(8)

      LOGICAL BIG, ERR

      ERR = .FALSE.
      TOL2 = TOL**2
      BETA = 1. - ALPHA

C  ITERATION LOOP

      DO 160 IT = 1, NIT
         BIG = .FALSE.

C  NODE LOOP

         NNT = 0
         DO 150 NODE = NNNOLD  +  1, NNN

C  CHECK FOR CONTINUATIONS,  BOUNDARY,  OR RELAXED NODE

            IF ((LXN(3, NODE) .GE.  0) .AND. (LXN(2, NODE) .GT. 0)
     &         .AND. (LXN(1, NODE) .GT. 0)) THEN
               NNT = NNT + 1

C  INITIALIZE

               KNUM = 0
               XA = 0.
               YA = 0.
               XL = 0.
               YL = 0.

C  PROCESS EACH LINE CONNECTED TO THIS NODE

               CALL GETLXN (MXND, LXN, NODE, LINES, KOUNT, ERR)
               IF (ERR) RETURN
               DO 100 IL = 1, KOUNT
                  L = LINES(IL)
                  N1 = NXL(1, L)
                  N2 = NXL(2, L)

C  FETCH ELEMENT DATA

                  IF (KXL(1, L) .GT. 0) CALL APALIB (MXND, XN, YN, LXK,
     &               NXL, KXL(1, L), NS1, AREA1, XCEN1, YCEN1, KNUM,
     &               KLIB, NLIB, ALIB, XCLIB, YCLIB)
                  IF (KXL(2, L) .GT. 0) CALL APALIB (MXND, XN, YN, LXK,
     &               NXL, KXL(2, L), NS2, AREA2, XCEN2, YCEN2, KNUM,
     &               KLIB, NLIB, ALIB, XCLIB, YCLIB)

C  GET FORCE DIRECTION MODULO PI RADIANS.
C  CORRECT FOR WRONG DIRECTION BY ALIGNING WITH THE VECTOR
C  FROM (XCEN1, YCEN1) TO (XCEN2, YCEN2).

                  XF = -(YN(N2) - YN(N1))
                  YF = XN(N2) - XN(N1)
                  DOT = XF*(XCEN2 - XCEN1) + YF*(YCEN2 - YCEN1)
                  IF (DOT  .LT.  0.) THEN
                     XF = -XF
                     YF = -YF
                  END IF

C  TAKE CARE OF ZERO AREAS

                  IF ((AREA1  .LE.  0) .OR. (AREA2  .LE.  0)) THEN
                     AREA1 = 1.0
                     AREA2 = 1.0
                  END IF

C  UPDATE AREA PULL VECTOR SUM

                  WEIGHT = (AREA2 - AREA1)/(AREA2 + AREA1)
                  XA = XA  +  WEIGHT*XF
                  YA = YA  +  WEIGHT*YF

C  UPDATE LAPLACIAN VECTOR SUM

                  NOE = N1 + N2 - NODE
                  DX = XN(NOE) - XN(NODE)
                  DY = YN(NOE) - YN(NODE)
                  XL = XL + DX
                  YL = YL + DY

C  UPDATE LEAST LENGTH

                  DIST2 = DX*DX  +  DY*DY
                  IF (IL .EQ. 1) DMIN2 = DIST2
                  DMIN2 = MIN(DMIN2, DIST2)
  100          CONTINUE

C  COMPUTE NET MOVEMENT VECTOR

               RK = 1.0/DBLE(KOUNT)
               XNET = (ALPHA*XA  +  BETA*XL)*RK
               YNET = (ALPHA*YA  +  BETA*YL)*RK

C  MOVE THE NODE

               YN(NODE) = YN(NODE)  +  YNET * RO
               XN(NODE) = XN(NODE)  +  XNET * RO

C  CHECK FOR SIGNIFICANT MOVEMENT

               IF (XNET*XNET + YNET*YNET .GT. TOL2*DMIN2) THEN

C  SIGNIFICANT MOVEMENT - REMOVE RELAXATION FLAGS

C  FIRST FROM DIRECTLY CONNECTED NODES

                  DO 110 IL = 1, KOUNT
                     L = LINES(IL)
                     NOE = NXL(1, L) + NXL(2, L) - NODE
                     LXN(3, NOE) = ABS(LXN(3, NOE))
  110             CONTINUE

C  NEXT FROM DIAGONALLY OPPOSITE NODES (MAX 8)

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

C  INSIGNIFICANT MOVEMENT
C  INSERT RELAXATION FLAG

               ELSE
                  LXN(3, NODE) = -ABS(LXN(3, NODE))
               END IF

            END IF
  150    CONTINUE
         IF (.NOT.BIG) GO TO 170
  160 CONTINUE
      IT = NIT

C  REMOVE ALL FLAGS

  170 CONTINUE
      DO 180 NODE = NNNOLD  +  1, NNN
         LXN(3, NODE) = ABS(LXN(3, NODE))
  180 CONTINUE

      RETURN
      END
