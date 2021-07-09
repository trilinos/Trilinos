C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION LXSYMB(SYM,NS,CH)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT
      CHARACTER*(*) SYM,CH
      CHARACTER*26 ALPHA,BETA
      DATA ALPHA/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      DATA BETA/'abcdefghijklmnopqrstuvwxyz'/

      SYM = ' '
      LS = LEN(SYM)
      CH = ILINE(JLINE:JLINE)
      IF (INDEX(ALPHA//BETA,CH).EQ.0) THEN
         LXSYMB = .FALSE.
         RETURN

      END IF

      LXSYMB = .TRUE.
      J0 = JLINE
 2540 IF (.NOT. (INDEX(ALPHA//BETA//'0123456789_$',CH).NE.0)) GO TO 2550
      JLINE = JLINE + 1
      CH = ILINE(JLINE:JLINE)
      GO TO 2540

 2550 CONTINUE
      SYM = ILINE(J0:JLINE-1)
      NS = MIN(LS,JLINE-J0)
      RETURN

      END
