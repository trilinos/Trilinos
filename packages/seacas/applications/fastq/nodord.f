C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE NODORD (NPNODE, XN, YN, LISTN, NUID, NNN)
C***********************************************************************

C  SUBROUTINE NODORD = ORDER THE NODE TABLE INTO INCREASING VALUES OF
C                      THE VARIABLE LISTN

C***********************************************************************

      DIMENSION LISTN (NPNODE)
      DIMENSION XN (NPNODE), YN (NPNODE), NUID (NPNODE)

      NN = NNN
      M = NN
  100 CONTINUE
      M =  (9 * M) / 16
      IF (M .LE. 0) RETURN
      M1 = M + 1
      DO 120 J = M1, NN
         L = J
         I = J - M
  110    CONTINUE
         IF (LISTN (L) .LT. LISTN (I)) THEN
            KLISTN = LISTN (I)
            KNUID = NUID (I)
            TXN = XN (I)
            TYN = YN (I)
            LISTN (I) = LISTN (L)
            NUID (I) = NUID (L)
            XN (I) = XN (L)
            YN (I) = YN (L)
            LISTN (L) = KLISTN
            NUID (L) = KNUID
            XN (L) = TXN
            YN (L) = TYN
            L = I
            I = I - M
            IF (I .GE. 1)GOTO 110
         ENDIF
  120 CONTINUE
      GOTO 100

      END
