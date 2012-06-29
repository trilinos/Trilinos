C $Id: chkkxl.f,v 1.1 1990/11/30 11:04:34 gdsjaar Exp $
C $Log: chkkxl.f,v $
C Revision 1.1  1990/11/30 11:04:34  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]CHKKXL.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE CHKKXL (MXND, LXK, KXL, LLL, ERR)
C***********************************************************************
C
C  SUBROUTINE CHKKXL = CHECKS TO SEE IF THE KXL COMPARES CORRECTLY TO
C                      THE LXK ARRAY
C
C***********************************************************************
C
      DIMENSION LXK (4, MXND), KXL (2, 3 * MXND)
C
      LOGICAL ERR
C
      ERR = .TRUE.
C
      DO 130 L = 1, LLL
         DO 120 IK = 1, 2
            K = KXL (IK, L)
            IF (K .NE. 0) THEN
               DO 100 I = 1, 4
                  IF (LXK (I, K) .EQ. L)GOTO 110
  100          CONTINUE
               WRITE ( * , 10000)IK, L, K
               RETURN
            ENDIF
  110       CONTINUE
  120    CONTINUE
  130 CONTINUE
      ERR = .FALSE.
C
      RETURN
C
10000 FORMAT ('KXL(', I4, ',', I4,') = ', I4,
     &   ' IS NOT IN LXK ARRAY  -  CHKKXL')
C
      END
