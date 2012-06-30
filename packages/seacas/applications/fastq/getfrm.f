C $Id: getfrm.f,v 1.1 1990/11/30 11:08:05 gdsjaar Exp $
C $Log: getfrm.f,v $
C Revision 1.1  1990/11/30 11:08:05  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]GETFRM.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE GETFRM (MXND, LINES, NL, NXL, NODE, N0, N2, NFROM)
C***********************************************************************
C
C  SUBROUTINE GETFRM = GETS THE NODES THAT THE CURRENT NODE CAME FROM
C
C***********************************************************************
C
      DIMENSION NXL(2, 3*MXND), LINES(NL)
C
      NFROM = 0
C
      IF (NL .EQ. 3) THEN
         DO 100 IL = 1, NL
            ILL = LINES (IL)
            IF (NXL (1, ILL) .EQ. NODE) THEN
               NTEST = NXL (2, ILL)
            ELSEIF (NXL (2, ILL) .EQ. NODE) THEN
               NTEST = NXL (1, ILL)
            ELSE
               CALL MESAGE ('** PROBLEMS IN GETFRM **')
               GOTO 110
            ENDIF
            IF ((NTEST .NE. N0) .AND. (NTEST .NE. N2)) THEN
               NFROM = NTEST
               GOTO 110
            ENDIF
  100    CONTINUE
      ENDIF
  110 CONTINUE
C
      RETURN
C
      END
