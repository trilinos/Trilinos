C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION LXGTQT(STR,NS,CH)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT
      CHARACTER*(*) STR
      CHARACTER CH

      CH = ILINE(JLINE:JLINE)
      LXGTQT = (CH.EQ.CHAR(39)) .OR. (CH.EQ.CHAR(34))
      IF (.NOT.LXGTQT) THEN
         RETURN

      END IF

      NS = 0
      J = JLINE
      LP = 1
 2470 IF (.NOT. (LP.GT.0)) GO TO 2480
      J = J + 1
      CH = ILINE(J:J)
      IF (CH.EQ.CHAR(39)) THEN
         LP = LP - 1

      ELSE IF (CH.EQ.CHAR(34)) THEN
         LP = LP - 1

      ELSE IF (CH.EQ.CHAR(0)) THEN
         LXGTQT = .FALSE.
         RETURN

      END IF

      NS = NS + 1
      STR(NS:NS) = CH
      GO TO 2470

 2480 CONTINUE
      NS = NS - 1
      JLINE = J + 1
      RETURN

      END
