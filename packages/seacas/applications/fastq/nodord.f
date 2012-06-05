C $Id: nodord.f,v 1.1 1990/11/30 11:12:46 gdsjaar Exp $
C $Log: nodord.f,v $
C Revision 1.1  1990/11/30 11:12:46  gdsjaar
C Initial revision
C
C
CC* FILE: [.RENUM]NODORD.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE NODORD (NPNODE, XN, YN, LISTN, NUID, NNN)
C***********************************************************************
C
C  SUBROUTINE NODORD = ORDER THE NODE TABLE INTO INCREASING VALUES OF
C                      THE VARIABLE LISTN
C
C***********************************************************************
C
      DIMENSION LISTN (NPNODE)
      DIMENSION XN (NPNODE), YN (NPNODE), NUID (NPNODE)
C
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
C
      END
