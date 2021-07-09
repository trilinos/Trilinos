C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION LXGTWH(CH)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT
      CHARACTER CH

      J0 = JLINE
      CH = ILINE(JLINE:JLINE)
 2490 IF (.NOT. (CH.EQ.' '.OR.CH.EQ.CHAR(9))) GO TO 2510
      JLINE = JLINE + 1
      CH = ILINE(JLINE:JLINE)
      GO TO 2490

 2510 CONTINUE
      LXGTWH = JLINE .NE. J0
      RETURN

      END
