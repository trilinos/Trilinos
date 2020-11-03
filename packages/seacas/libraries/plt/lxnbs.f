C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION LXNBS(LINE,L)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT
      CHARACTER*(*) LINE
      CHARACTER CH
      LOGICAL LDUM
      LOGICAL LXGTWH

      LDUM = LXGTWH(CH)
      CH = ILINE(JLINE:JLINE)
      J0 = JLINE
 2520 IF (.NOT. (INDEX(' '//CHAR(9)//CHAR(0),CH).EQ.0)) GO TO 2530
      JLINE = JLINE + 1
      CH = ILINE(JLINE:JLINE)
      GO TO 2520

 2530 CONTINUE
      L = JLINE - J0
      IF (J0.EQ.JLINE) THEN
         LXNBS = .FALSE.
         RETURN

      END IF

      LINE = ILINE(J0:JLINE-1)
      LXNBS = .TRUE.
      RETURN

      END
