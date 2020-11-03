C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE CONDNO (MXND, NODES, QRAT, SRAT, COND, SIDES, XN, YN,
     &   LSIDE)
C***********************************************************************

C  SUBROUTINE CONDNO = COMPUTES EVALUATOR FUNCTIONS FOR RESTRUCTURING

C***********************************************************************

C  VARIABLES USED:
C     LSIDE = .TRUE. IF SIDES ARRAY IS TO BE FILLED
C     NODES = THE FOUR NODES OF THE ELEMENT
C     QRAT  = THE RATIO OF SMALLEST TO LARGEST ANGLE
C     SRAT  = THE RATIO OF SMALLEST TO LARGEST SIDE
C     COND  = SRAT*QRAT

C***********************************************************************

      DIMENSION NODES (4), SIDES (4), XN (MXND), YN (MXND)

      LOGICAL LSIDE

      N1 = NODES (1)
      N2 = NODES (2)
      N3 = NODES (3)
      N4 = NODES (4)
      SS1 = (XN (N1) - XN (N2)) **2 +  (YN (N1) - YN (N2)) **2
      SS2 = (XN (N2) - XN (N3)) **2 +  (YN (N2) - YN (N3)) **2
      SS3 = (XN (N3) - XN (N4)) **2 +  (YN (N3) - YN (N4)) **2
      SS4 = (XN (N4) - XN (N1)) **2 +  (YN (N4) - YN (N1)) **2
      AMAX = AMAX1 (SS1, SS3)
      AMIN = AMIN1 (SS1, SS3)
      BMAX = AMAX1 (SS2, SS4)
      BMIN = AMIN1 (SS2, SS4)
      IF (AMIN * BMIN .GT. 0.0) THEN
         SRAT = SQRT (SQRT (AMAX * BMAX /  (AMIN * BMIN)))
      ELSE
         SRAT = 1.0E10
      ENDIF
      COND = QRAT * SRAT
      IF (LSIDE) THEN
         SIDES (1) = SQRT (SS1)
         SIDES (2) = SQRT (SS2)
         SIDES (3) = SQRT (SS3)
         SIDES (4) = SQRT (SS4)
      ENDIF

      RETURN

      END
