C $Id: cross.f,v 1.1 1990/11/30 11:05:31 gdsjaar Exp $
C $Log: cross.f,v $
C Revision 1.1  1990/11/30 11:05:31  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]CROSS.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE CROSS (J1, J2, I1, I2, JXI, IXJ, ISTART, ICLEAR,
     &   NOROOM, ERR)
C***********************************************************************
C
C  SUBROUTINE CROSS = CREATE OR ADD TO THE CROSS - REFERENCE ARRAY FOR
C                     JXI (J1, J2) IN IXJ (I1, I2)
C
C***********************************************************************
C
C  NOTE:
C     THE NEW ITEMS MUST BEGIN AT J1=1,  J2=ISTART.
C     THE CROSS REFERENCE ARRAY WILL BE CLEARED FROM I1=1,  I2=ICLEAR
C     TO THE END OF THE ARRAY.
C
C***********************************************************************
C
      DIMENSION JXI (J1, J2), IXJ (I1, I2)
C
      LOGICAL ERR, NOROOM
C
C  CLEAR
C
      ERR = .TRUE.
      NOROOM = .FALSE.
      DO 110 J = ICLEAR, I2
         DO 100 I = 1, I1
            IXJ (I, J)  =  0
  100    CONTINUE
  110 CONTINUE
C
C  REFILE EACH ITEM
C
      DO 150 J = ISTART, J2
         DO 140 I = 1, J1
            L = IABS (JXI (I, J))
            IF (L .NE. 0) THEN
               IF (L .GT. I2) THEN
                  WRITE ( * , 10000)L, I2
                  RETURN
               ENDIF
C
C  FIND EMPTY SPOT FOR THIS ITEM
C
               DO 120 K = 1, I1
                  KK = K
                  IF (IXJ (K, L) .EQ. 0)GO TO 130
  120          CONTINUE
               CALL MESAGE ('NO ROOM FOR REFERENCE - ERROR IN CROSS')
               NOROOM = .TRUE.
               RETURN
  130          CONTINUE
C
C  FILE THIS ITEM
C
               IXJ (KK, L)  =  J
            ENDIF
  140    CONTINUE
  150 CONTINUE
      ERR = .FALSE.
      RETURN
C
10000 FORMAT (' OUT-OF-BOUNDS REFERENCE IN CROSS (INDEX = ', I5,
     &   ', MAX = ', I5, ')')
      END
