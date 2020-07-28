C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE LXCLN
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT

      K = INDEX(ILINE(JLINE:),CHAR(0))
      IF (K.EQ.0) THEN
         JLINE = 504

      ELSE
         JLINE = MIN(JLINE+K,504)
      END IF

      RETURN

      END
