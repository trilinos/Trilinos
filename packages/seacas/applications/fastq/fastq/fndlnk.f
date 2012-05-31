C $Id: fndlnk.f,v 1.1 1990/11/30 11:07:43 gdsjaar Exp $
C $Log: fndlnk.f,v $
C Revision 1.1  1990/11/30 11:07:43  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]FNDLNK.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE FNDLNK (MXND, LXK, NXL, K, N1, N2, L, ERR)
C***********************************************************************
C
C  SUBROUTINE FNDLNK = FIND THE LINE IN ELEMENT K WITH NODES N1 AND N2
C
C***********************************************************************
C
      DIMENSION LXK (4, MXND), NXL (2, 3 * MXND)
C
      LOGICAL ERR
C
      ERR = .FALSE.
      DO 100 I = 1, 4
         LL = LXK (I, K)
         M1 = NXL (1, LL)
         M2 = NXL (2, LL)
         IF ( ( (M1 .EQ. N1) .AND. (M2 .EQ. N2)) .OR.
     &      ( (M2 .EQ. N1) .AND. (M1 .EQ. N2) ) ) THEN
            L = LL
            RETURN
         ENDIF
  100 CONTINUE
      L = 0
      ERR = .TRUE.
      WRITE ( * , 10000) K, N1, N2
10000 FORMAT (' IN FNDLNK, NO LINE CAN BE FOUND FOR K, N1, N2: ', 3I5)
      RETURN
      END
