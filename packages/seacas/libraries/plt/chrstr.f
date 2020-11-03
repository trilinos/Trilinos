C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CHRSTR(LINE1,LINE2,L)
      CHARACTER*(*) LINE1,LINE2
      CHARACTER CH

      K = LEN(LINE1)
      J = 1
 2150 IF (.NOT. (J.LT.K)) GO TO 2170
      CH = LINE1(J:J)
      IF (CH.EQ.' ') THEN
         GO TO 2160

      END IF

      IF (CH.EQ.CHAR(9)) THEN
         GO TO 2160

      END IF

      GO TO 2170

 2160 J = J + 1
      GO TO 2150

 2170 CONTINUE
      IF (J.GT.K) THEN
         L = 0
         LINE2 = ' '
         RETURN

      END IF

      LINE2(1:K-J+1) = LINE1(J:K)
      CALL CHRTRM(LINE2(1:K-J+1),L)
      RETURN

      END
