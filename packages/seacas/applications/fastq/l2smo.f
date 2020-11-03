C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE L2SMO (MXND, XN, YN, NXL, LXN, NNN, NNNOLD, NIT, EPS,
     &   RO)
C***********************************************************************

C  SUBROUTINE L2SMO = LAPLACIAN SQUARED METHOD OF MESH SMOOTHING

C***********************************************************************

C  NOTE:
C     THIS METHOD IS LIKE LAPLACIAN SMOOTHING EXCEPT EACH VECTOR
C     POINTING TO A NEIGHBOR NODE HAS A LENGTH OF
C      (DISTANCE TO THAT NODE)**2.

C***********************************************************************

C  VARIABLES USED:
C     NIT = MAX NUMBER OF ITERATIONS TO DO
C     EPS = NODE MOVEMENT TOLERANCE FOR CONVERGENCE
C     RO  = UNDER OR OVER-RELAXATION FACTOR.

C***********************************************************************

      DIMENSION LINES (20)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND), XN (MXND), YN (MXND)

      LOGICAL BIG, ERR

      EPS2 = (EPS * RO) ** 2

C  ITERATION LOOP

      DO 120 IT = 1, NIT
         BIG = .FALSE.

C  NODE LOOP

         DO 110 NODE = NNNOLD + 1, NNN
            IF ((LXN (1, NODE).GT.0) .AND. (LXN (2, NODE).GT.0))THEN
               FX = 0.
               FY = 0.
               SL = 0.

C  LOOP THROUGH ALL LINES CONNECTED TO NODE

               CALL GETLXN (MXND, LXN, NODE, LINES, KOUNT, ERR)
               IF (ERR) RETURN
               DO 100 IL = 1, KOUNT
                  L = LINES (IL)
                  NEND = NXL (1, L) + NXL (2, L) - NODE
                  DX = XN (NEND) - XN (NODE)
                  DY = YN (NEND) - YN (NODE)
                  AL2 = DX * DX + DY * DY
                  AL = SQRT (AL2)
                  FX = FX + DX * AL
                  FY = FY + DY * AL
                  SL = SL + AL
  100          CONTINUE

C  MOVE THE NODE

               DELX = RO * FX/SL
               DELY = RO * FY/SL
               XN (NODE) = XN (NODE) + DELX
               YN (NODE) = YN (NODE) + DELY
               IF (DELX ** 2 + DELY ** 2 .GT. EPS2) BIG = .TRUE.
            ENDIF
  110    CONTINUE
         IF (.NOT.BIG)RETURN
  120 CONTINUE
      RETURN
      END
