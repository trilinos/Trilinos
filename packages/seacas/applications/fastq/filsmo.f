C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE FILSMO (MXND, MLN, XN, YN, ZN, LXK, KXL, NXL, LXN,
     &   LLL, NNN, NNN2, LNODES, BNSIZE, NLOOP, XMIN, XMAX, YMIN, YMAX,
     &   ZMIN, ZMAX, DEV1, KREG)
C***********************************************************************

C  SUBROUTINE FILSMO = MESH SMOOTHING DONE BY ISOPARAMETRIC/EQUAL
C                      ANGULAR SMOOTHING OF THE ADDED INTERIOR (FREE)
C                      BOUNDARY ROW AND THEN A LENGTH-WEIGHTED/EQUAL
C                      ANGULAR BOUNDARY LAPLACIAN OF THE INTERIOR NODES.
C                      THE FREE BOUNDARY IS FINALLY SMOOTHED AGAIN.

C***********************************************************************

C  VARIABLES USED:
C     WFAC = WEIGHT (0. = LAPLACIAN, 1. = ISOPARAMETRIC)
C     NIT  = THE MAX NUMBER OF ITERATIONS TO DO.
C     EPS  = MINIMUM DISTANCE NODES MUST MOVE TO CONTINUE ITERATIONS
C     RO   = AN UNDER- OR OVER-RELAXATION FACTOR (NORMALLY 1.0)

C***********************************************************************

      COMMON /TIMING/ TIMEA, TIMEP, TIMEC, TIMEPC, TIMEAJ, TIMES

      DIMENSION XN(MXND), YN(MXND), ZN(MXND)
      DIMENSION LXN(4, MXND), NXL(2, 3*MXND)
      DIMENSION LXK(4, MXND), KXL(2, 3*MXND)
      DIMENSION LINES(20), LNODES (MLN, MXND), BNSIZE (2, MXND)

      LOGICAL BIG, ERR, GRAPH, DONE

      CHARACTER*3 DEV1

      CALL GETIME (TIME1)
      GRAPH = .FALSE.
      DONE = .FALSE.
      WT = 10.

      NIT = MAX0 (5 * NLOOP, 40)
      TOL = .03
      VRO = 1.
      RO = 1.
      WFAC = 1.0
      WFAC2 = .5
      CALL MNORM  (MXND,  XN,  YN,  NXL,  LLL,  STDLEN)
      EPS  =  TOL * STDLEN
      IF (RO .LT. 0.01) RO = 1.
      EPS2 = (EPS * RO)**2

C  FIRST SMOOTH THE ADDED ROW

      IF (NLOOP .GT. 0) THEN
         CALL ROWSMO (MXND, MLN, XN, YN, ZN, LXK, KXL, NXL, LXN, NNN,
     &      WFAC, WFAC2, NIT, EPS, RO, NNN2, LNODES, BNSIZE, LLL,
     &      GRAPH, XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, DEV1, KREG)
      ENDIF

C  NOW SMOOTH THE INTERIOR NODES

C  ITERATION LOOP

      DO 140 IT = 1, NIT
         BIG = .FALSE.

C  NODE LOOP

         DO 130 NODE = 1, NNN
            IF ( (LXN (1, NODE) .GT. 0) .AND.
     &         (LXN (2, NODE) .GT. 0) .AND.
     &         (LNODES (4, NODE) .EQ. - 2) ) THEN
               DONE = .TRUE.
               FX = 0.
               FY = 0.
               SL = 0.
               VL = 0.

C  LOOP THROUGH ALL LINES CONNECTED TO NODE

               CALL GETLXN (MXND, LXN, NODE, LINES, KOUNT, ERR)
               IF (ERR) GOTO 150
               DO 100 IL = 1, KOUNT
                  L = LINES (IL)
                  NEND = NXL (1, L) + NXL (2, L) - NODE
                  DX = XN (NEND) - XN (NODE)
                  DY = YN (NEND) - YN (NODE)
                  AL = SQRT (DX * DX + DY * DY)

C  CHECK FOR A BOUNDARY NODE AT THE OTHER END
C  OF THE LINE - TRY TO AVERAGE ANGULAR ERRORS WITH THE BOUNDARY WHERE
C  POSSIBLE - THIS MEANS ADDING IN AN EXTRA VECTOR TO PULL THE NODE
C  BACK TO WHERE IT OUGHT TO BE TO BE AT EQUAL ANGLES

                  IF (LXN (2, NEND) .LT. 0) THEN
                     CALL SETN02 (MXND, NXL, LXK, KXL, L, NEND, NODE,
     &                  N0, N2)
                     CALL EQLANG (MXND, XN, YN, LXN, NODE, N0, N2,
     &                  NEND, AL, VRO, VXDEL, VYDEL)
                     VL = SQRT (VXDEL * VXDEL + VYDEL * VYDEL)
                     FX = FX + (VXDEL * WT * VL)
                     FY = FY + (VYDEL * WT * VL)
                     SL = SL + VL * WT
                  ENDIF
                  FX = FX + DX * AL
                  FY = FY + DY * AL
                  SL = SL + AL
  100          CONTINUE

C  MOVE THE NODE

               DELX = RO * FX/SL
               DELY = RO * FY/SL

C  ERASE THE NODE'S LINES IF GRAPH IS ON

               IF (GRAPH) THEN
                  CALL LCOLOR('BLACK')
                  DO 110 II = 1, KOUNT
                     IDRAW = LINES(II)
                     NODE1 = NXL (1, IDRAW)
                     NODE2 = NXL (2, IDRAW)
                     CALL D2NODE (MXND, XN, YN, NODE1, NODE2)
  110             CONTINUE
                  CALL LCOLOR ('WHITE')
               ENDIF

               XN (NODE) = XN (NODE)+DELX
               YN (NODE) = YN (NODE)+DELY

C  REPLOT THE NODE'S LINES IF GRAPH IS ON

               IF (GRAPH) THEN
                  DO 120 II = 1, KOUNT
                     IDRAW = LINES(II)
                     NODE1 = NXL (1, IDRAW)
                     NODE2 = NXL (2, IDRAW)
                     CALL D2NODE (MXND, XN, YN, NODE1, NODE2)
  120             CONTINUE
                  CALL SFLUSH
               ENDIF
               IF (DELX ** 2 + DELY ** 2 .GT. EPS2) BIG = .TRUE.
            ENDIF
  130    CONTINUE
         IF (.NOT.BIG) GOTO 150
  140 CONTINUE
  150 CONTINUE

C  NOW RESMOOTH THE ADDED ROW IF THE MESH HAS CHANGED INTERNALLY

      IF ((NLOOP .GT. 0) .AND. (DONE)) THEN
         CALL ROWSMO (MXND, MLN, XN, YN, ZN, LXK, KXL, NXL, LXN, NNN,
     &      WFAC, WFAC2, NIT, EPS, RO, NNN2, LNODES, BNSIZE, LLL,
     &      GRAPH, XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, DEV1, KREG)
      ENDIF

C  NOW RESET ALL THE NODES AS BEING SMOOTHED

      DO 160 I = 1, NNN
         LNODES (4, I) = IABS (LNODES (4, I))
  160 CONTINUE

      CALL GETIME (TIME2)
      TIMES = TIMES + TIME2 - TIME1
      RETURN

      END
