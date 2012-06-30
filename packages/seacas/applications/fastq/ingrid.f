C $Id: ingrid.f,v 1.1 1990/11/30 11:09:38 gdsjaar Exp $
C $Log: ingrid.f,v $
C Revision 1.1  1990/11/30 11:09:38  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]INGRID.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE INGRID (MSNAP, SNAPDX, NSNAP, II, RIN, IFOUND, ERR)
C***********************************************************************
C
C  SUBROUTINE INGRID = INPUTS A X OR Y GRID CARD
C
C***********************************************************************
C
      DIMENSION SNAPDX(2, MSNAP), NSNAP(2), RIN(IFOUND)
      LOGICAL ERR
C
      ERR = .FALSE.
      IF (II .LT. 1 .OR. II .GT. 2) THEN
         ERR = .TRUE.
         WRITE (*, 10000) II
         GO TO 110
      END IF
C
      DO 100 I = 1, IFOUND
         CALL ADDSNP (MSNAP, SNAPDX, NSNAP, II, RIN(I), ERR)
         IF (ERR) GO TO 110
  100 CONTINUE
C
  110 CONTINUE
      RETURN
C
10000 FORMAT (' GRID INDEX OUT-OF-RANGE [1,2] IN INGRID: ', I3)
      END
