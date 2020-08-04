C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE CCROSS (J1, J2, I1, I2, JXI, IXJ, ISTART, ICLEAR,
     &   NOROOM, ERR)
C***********************************************************************

C  SUBROUTINE CROSS = CREATE OR ADD TO THE CROSS - REFERENCE ARRAY FOR
C                     JXI (J1, J2) IN IXJ (I1, I2)

C***********************************************************************

C  NOTE:
C     THE NEW ITEMS MUST BEGIN AT J1=1,  J2=ISTART.
C     THE CROSS REFERENCE ARRAY WILL BE CLEARED FROM I1=1,  I2=ICLEAR
C     TO THE END OF THE ARRAY.

C***********************************************************************

      DIMENSION JXI (J1, J2), IXJ (I1, I2)

      LOGICAL ERR, NOROOM

C  CLEAR

      ERR = .TRUE.
      NOROOM = .FALSE.
      DO 110 J = ICLEAR, I2
         DO 100 I = 1, I1
            IXJ (I, J)  =  0
  100    CONTINUE
  110 CONTINUE

C  REFILE EACH ITEM

      DO 150 J = ISTART, J2
         DO 140 I = 1, J1
            L = IABS (JXI (I, J))
            IF (L .NE. 0) THEN
               IF (L .GT. I2) THEN
                  WRITE ( * , 10000)L, I2
                  RETURN
               ENDIF

C  FIND EMPTY SPOT FOR THIS ITEM

               DO 120 K = 1, I1
                  KK = K
                  IF (IXJ (K, L) .EQ. 0)GO TO 130
  120          CONTINUE
               CALL MESAGE ('NO ROOM FOR REFERENCE - ERROR IN CROSS')
               NOROOM = .TRUE.
               RETURN
  130          CONTINUE

C  FILE THIS ITEM

               IXJ (KK, L)  =  J
            ENDIF
  140    CONTINUE
  150 CONTINUE
      ERR = .FALSE.
      RETURN

10000 FORMAT (' OUT-OF-BOUNDS REFERENCE IN CROSS (INDEX = ', I5,
     &   ', MAX = ', I5, ')')
      END
