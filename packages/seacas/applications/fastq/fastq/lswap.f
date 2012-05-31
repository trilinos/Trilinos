C $Id: lswap.f,v 1.1 1990/11/30 11:11:33 gdsjaar Exp $
C $Log: lswap.f,v $
C Revision 1.1  1990/11/30 11:11:33  gdsjaar
C Initial revision
C
CC* FILE: [.QMESH]LSWAP.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE LSWAP (MXND, LXK, KXL, K1, L1, K2, L2, ERR)
C***********************************************************************
C
C  SUBROUTINE LSWAP = EXCHANGE LINE L1 IN ELEMENT K1 WITH LINE L2 IN
C                     ELEMENT K2
C
C***********************************************************************
C
      DIMENSION LXK (4, MXND), KXL (2, 3*MXND)
C
      LOGICAL ERR
      ERR = .TRUE.
C
C  INSERT L2 FOR L1
C
      DO 130 I = 1, 4
         IF (LXK (I, K1) .EQ. L1) THEN
            LXK (I, K1) = L2
C
C  INSERT L1 FOR L2
C
            DO 120 J = 1, 4
               IF (LXK (J, K2) .EQ. L2) THEN
                  LXK (J, K2) = L1
C
C  INSERT K2 FOR K1
C
                  DO 110 K = 1, 2
                     IF (KXL (K, L1) .EQ. K1) THEN
                        KXL (K, L1) = K2
C
C  INSERT K1 FOR K2
C
                        DO 100 L = 1, 2
                           IF (KXL (L, L2) .EQ. K2) THEN
                              KXL (L, L2) = K1
C
C  EVERYTHING INSERTED OK
C
                              ERR = .FALSE.
                              RETURN
                           ENDIF
  100                   CONTINUE
                        WRITE ( * , 10000)K1, L1, K2, L2
                        RETURN
                     ENDIF
  110             CONTINUE
                  WRITE ( * , 10000)K1, L1, K2, L2
                  RETURN
               ENDIF
  120       CONTINUE
            WRITE ( * , 10000)K1, L1, K2, L2
            RETURN
         ENDIF
  130 CONTINUE
      WRITE ( * , 10000)K1, L1, K2, L2
C
      RETURN
C
10000 FORMAT (' ERROR IN LSWAP.  K1, L1, K2, L2 :', 4I5)
C
      END
