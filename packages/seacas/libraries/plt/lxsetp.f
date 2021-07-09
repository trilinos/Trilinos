C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE LXSETP(LINE)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT
      CHARACTER*(*) LINE
      CHARACTER*255 LOCLIN
      CHARACTER*80 TMPLIN
      INTEGER CHRLEN

      IF (LXINIT.NE.12345) THEN
         CALL LXRST
      END IF

      LOCLIN = LINE
      L = CHRLEN(LOCLIN)
      K = JLINE - L - 1
      IF (K.LE.0) THEN
         TMPLIN = 'Buffer overflow in lxsetp: '//LINE(1:L)
         CALL LXERR(TMPLIN,3)
         RETURN

      END IF

      IF (L.GT.0) THEN
         ILINE(K:JLINE-1) = LINE(1:L)//CHAR(0)

      ELSE
         ILINE(K:JLINE-1) = CHAR(0)
      END IF

      JLINE = K
      RETURN

      END
