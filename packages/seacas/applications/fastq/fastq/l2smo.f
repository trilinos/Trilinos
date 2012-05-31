C $Id: l2smo.f,v 1.1 1990/11/30 11:10:53 gdsjaar Exp $
C $Log: l2smo.f,v $
C Revision 1.1  1990/11/30 11:10:53  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]L2SMO.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE L2SMO (MXND, XN, YN, NXL, LXN, NNN, NNNOLD, NIT, EPS,
     &   RO)
C***********************************************************************
C
C  SUBROUTINE L2SMO = LAPLACIAN SQUARED METHOD OF MESH SMOOTHING
C
C***********************************************************************
C
C  NOTE:
C     THIS METHOD IS LIKE LAPLACIAN SMOOTHING EXCEPT EACH VECTOR
C     POINTING TO A NEIGHBOR NODE HAS A LENGTH OF
C      (DISTANCE TO THAT NODE)**2.
C
C***********************************************************************
C
C  VARIABLES USED:
C     NIT = MAX NUMBER OF ITERATIONS TO DO
C     EPS = NODE MOVEMENT TOLERANCE FOR CONVERGENCE
C     RO  = UNDER OR OVER-RELAXATION FACTOR.
C
C***********************************************************************
C
      DIMENSION LINES (20)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND), XN (MXND), YN (MXND)
C
      LOGICAL BIG, ERR
C
      EPS2 = (EPS * RO) ** 2
C
C  ITERATION LOOP
C
      DO 120 IT = 1, NIT
         BIG = .FALSE.
C
C  NODE LOOP
C
         DO 110 NODE = NNNOLD + 1, NNN
            IF ((LXN (1, NODE).GT.0) .AND. (LXN (2, NODE).GT.0))THEN
               FX = 0.
               FY = 0.
               SL = 0.
C
C  LOOP THRU ALL LINES CONNECTED TO NODE
C
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
C
C  MOVE THE NODE
C
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
