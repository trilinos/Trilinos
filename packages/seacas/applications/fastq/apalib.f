C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE APALIB (MXND, XN, YN, LXK, NXL, K, NODES, AREA, XCEN,
     &   YCEN, KNUM, KLIB, NLIB, ALIB, XCLIB, YCLIB)
C***********************************************************************

C  SUBROUTINE APALIB = LIBRARY OF ELEMENT DATA USED TO AVOID DUPLICATE
C                      COMPUTATIONS

C***********************************************************************

      DIMENSION NODES (4), KLIB (8), NLIB (4, 8)
      DIMENSION ALIB (8), XCLIB (8), YCLIB (8)
      DIMENSION XN (MXND), YN (MXND), LXK (4, MXND), NXL (2, 3 * MXND)

      LOGICAL CCW

C  SEARCH LIBRARY

      IF (KNUM .GT. 0) THEN
         DO 110 I = 1, KNUM
            IF  (K - KLIB (I) .EQ. 0) THEN

C  FETCH FROM LIBRARY

               IK = I
               DO 100 J = 1, 4
                  NODES (J) = NLIB (J, IK)
  100          CONTINUE
               AREA = ALIB (IK)
               XCEN = XCLIB (IK)
               YCEN = YCLIB (IK)
               RETURN
            ENDIF
  110    CONTINUE
      ENDIF

C  COMPUTE NEW DATA

      CCW = .TRUE.
      CALL GNXKA (MXND, XN, YN, K, NODES, AREA, LXK, NXL, CCW)
      N1 = NODES (1)
      N2 = NODES (2)
      N3 = NODES (3)
      N4 = NODES (4)
      XCEN =  (XN (N1) + XN (N2) + XN (N3) + XN (N4)) * 0.25
      YCEN =  (YN (N1) + YN (N2) + YN (N3) + YN (N4)) * 0.25
      IF (KNUM .GE. 8) RETURN

C  FILE NEW DATA IN LIBRARY

      KNUM = KNUM + 1
      DO 120 I = 1, 4
         NLIB (I, KNUM) = NODES (I)
  120 CONTINUE
      KLIB (KNUM) = K
      ALIB (KNUM) = AREA
      XCLIB (KNUM) = XCEN
      YCLIB (KNUM) = YCEN
      RETURN

      END
