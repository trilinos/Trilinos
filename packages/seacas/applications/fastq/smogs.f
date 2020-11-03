C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE SMOGS (MXND, XN, YN, NXL, LXN, NNN, NNNOLD, NIT, EPS,
     &   RO)
C***********************************************************************

C  SUBROUTINE SMOGS  =  MESH SMOOTHING BY LAPLACE-S USING GAUSS-SEIDEL

C***********************************************************************

C  VARIABLES USED:
C     NIT   =  THE MAX NUMBER OF ITERATIONS TO DO.
C     EPS   =  MINIMUM DISTANCE NODES MUST MOVE TO CONTINUE ITERATIONS
C     RO    =  AN UNDER- OR OVER-RELAXATION FACTOR (NORMALLY 1.0)

C***********************************************************************

      DIMENSION LINES(20)
      DIMENSION LXN(4, MXND), NXL(2, 3*MXND), XN(MXND), YN(MXND)

      LOGICAL BIG, ERR

      IF (RO .LT. 0.01) RO = 1.
      EPS2 = (EPS*RO)**2

C  ITERATION LOOP

      DO 120 IT = 1, NIT
         BIG = .FALSE.

C  NODE LOOP

         DO 110 NODE = NNNOLD + 1, NNN

C  SKIP CONTINUATION AND BOUNDARY LINES

            IF ((LXN(1, NODE) .GT. 0) .AND. (LXN(2, NODE) .GT. 0)) THEN

C  SUM COORDINATES OF ALL NEIGHBORING NODES

               SUMX = 0.0
               SUMY = 0.0
               CALL GETLXN (MXND, LXN, NODE, LINES, KOUNT, ERR)

C  IGNORE ERR BECAUSE IT IS ALREADY TAKEN CARE OF IN THE SKIP

               DO 100 IL = 1, KOUNT
                  L = LINES(IL)
                  IM = NXL(1, L) + NXL(2, L) - NODE
                  SUMX = SUMX + XN(IM)
                  SUMY = SUMY + YN(IM)
  100          CONTINUE

C  REDEFINE THIS NODE - S COORDINATES

               SUMX = SUMX/DBLE(KOUNT)
               SUMY = SUMY/DBLE(KOUNT)
               XDEL = RO*(SUMX - XN(NODE))
               YDEL = RO*(SUMY - YN(NODE))
               XN(NODE) = XN(NODE) + XDEL
               YN(NODE) = YN(NODE) + YDEL

C  CHECK FOR CONVERGENCE

               IF ((XDEL*XDEL + YDEL*YDEL) .GT. EPS2) BIG = .TRUE.
            ENDIF
  110    CONTINUE

C  IF NO SIGNIFICANT MOVEMENTS OCCURRED,  RETURN

         IF (.NOT.BIG) RETURN
  120 CONTINUE
      RETURN
      END
