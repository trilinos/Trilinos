C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE LSWAP (MXND, LXK, KXL, K1, L1, K2, L2, ERR)
C***********************************************************************

C  SUBROUTINE LSWAP = EXCHANGE LINE L1 IN ELEMENT K1 WITH LINE L2 IN
C                     ELEMENT K2

C***********************************************************************

      DIMENSION LXK (4, MXND), KXL (2, 3*MXND)

      LOGICAL ERR
      ERR = .TRUE.

C  INSERT L2 FOR L1

      DO 130 I = 1, 4
         IF (LXK (I, K1) .EQ. L1) THEN
            LXK (I, K1) = L2

C  INSERT L1 FOR L2

            DO 120 J = 1, 4
               IF (LXK (J, K2) .EQ. L2) THEN
                  LXK (J, K2) = L1

C  INSERT K2 FOR K1

                  DO 110 K = 1, 2
                     IF (KXL (K, L1) .EQ. K1) THEN
                        KXL (K, L1) = K2

C  INSERT K1 FOR K2

                        DO 100 L = 1, 2
                           IF (KXL (L, L2) .EQ. K2) THEN
                              KXL (L, L2) = K1

C  EVERYTHING INSERTED OK

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

      RETURN

10000 FORMAT (' ERROR IN LSWAP.  K1, L1, K2, L2 :', 4I5)

      END
