C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE GETFRM (MXND, LINES, NL, NXL, NODE, N0, N2, NFROM)
C***********************************************************************

C  SUBROUTINE GETFRM = GETS THE NODES THAT THE CURRENT NODE CAME FROM

C***********************************************************************

      DIMENSION NXL(2, 3*MXND), LINES(NL)

      NFROM = 0

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

      RETURN

      END
