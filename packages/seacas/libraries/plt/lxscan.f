C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION LXSCAN(DELIM,REM,L,CH)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT
      CHARACTER*(*) DELIM
      CHARACTER*(*) REM
      INTEGER L
      CHARACTER CH

      REM = ' '
      L = 0
      LXSCAN = .FALSE.
 2410 CONTINUE
      CH = ILINE(JLINE:JLINE)
      IF (CH.EQ.CHAR(0)) THEN
         RETURN

      END IF

      JLINE = JLINE + 1
      LXSCAN = (INDEX(DELIM,CH).NE.0)
      IF (LXSCAN) THEN
         RETURN

      END IF

      L = L + 1
      IF (L.GT.LEN(REM)) THEN
         GO TO 2420

      END IF

      REM(L:L) = CH
 2420 GO TO 2410

      END
