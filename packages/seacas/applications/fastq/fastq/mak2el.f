C $Id: mak2el.f,v 1.1 1990/11/30 11:11:48 gdsjaar Exp $
C $Log: mak2el.f,v $
C Revision 1.1  1990/11/30 11:11:48  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]MAK2EL.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE MAK2EL (MP, MXNPER, MXND, NNN0, NNN, KKK, X, Y, NID,
     &   XN, YN, NUID, LXK, COOR, IP3)
C***********************************************************************
C
C  SUBROUTINE MAK2EL = GENERATES  (ADDS) ELEMENT CONNECTIVITY FOR 2 NODES
C
C***********************************************************************
C
      DIMENSION COOR (2, MP), X (MXNPER), Y (MXNPER), NID (MXNPER)
      DIMENSION XN (MXND), YN (MXND), NUID (MXND), LXK (4, MXND)
C
C  PUT NODES AND NUID'S INTO THE PROPER LOCATIONS
C
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
C
      XN (NNN) = X (KOUNT + 1)
      YN (NNN) = Y (KOUNT + 1)
      NUID (NNN) = NID (KOUNT + 1)
C
      RETURN
C
10000 FORMAT (' ** WARNING **  -  COLINEAR REFERENCE NODE FOR ELEMENT:',
     &   I5)
      END
