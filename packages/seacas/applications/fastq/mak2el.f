C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE MAK2EL (MP, MXNPER, MXND, NNN0, NNN, KKK, X, Y, NID,
     &   XN, YN, NUID, LXK, COOR, IP3)
C***********************************************************************

C  SUBROUTINE MAK2EL = GENERATES  (ADDS) ELEMENT CONNECTIVITY FOR 2 NODES

C***********************************************************************

      DIMENSION COOR (2, MP), X (MXNPER), Y (MXNPER), NID (MXNPER)
      DIMENSION XN (MXND), YN (MXND), NUID (MXND), LXK (4, MXND)

C  PUT NODES AND NUID'S INTO THE PROPER LOCATIONS

      KOUNT = 0
      DO 100 I = NNN0, NNN-1
         KOUNT = KOUNT + 1
         KKK = KKK + 1
         XN (I) = X (KOUNT)
         YN (I) = Y (KOUNT)
         NUID (I) = NID (KOUNT)
         LXK (1, KKK) = I
         LXK (2, KKK) = I + 1
         LXK (3, KKK) = 0
         LXK (4, KKK) = 0
         IF (IP3.GT.0)THEN
            X1 = X (I + 1) - X (I)
            Y1 = Y (I + 1) - Y (I)
            X2 = X (I + 1) - COOR (1, IP3)
            Y2 = Y (I + 1) - COOR (2, IP3)
            CROSSP =  (X1 * Y2) -  (Y1 * X2)
            IF (CROSSP.GT.0)THEN
               LXK (1, KKK) = I + 1
               LXK (2, KKK) = I
            ENDIF
            IF (CROSSP * CROSSP .LT.
     &         (.01 * ((X1 * X1) +  (Y1 * Y1)) *
     &         ((X2 * X2) + (Y2 * Y2)))) WRITE (*, 10000) KKK
         ENDIF
  100 CONTINUE

      XN (NNN) = X (KOUNT + 1)
      YN (NNN) = Y (KOUNT + 1)
      NUID (NNN) = NID (KOUNT + 1)

      RETURN

10000 FORMAT (' ** WARNING **  -  COLINEAR REFERENCE NODE FOR ELEMENT:',
     &   I5)
      END
