C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE LXREM(LINE,L)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT
      CHARACTER*(*) LINE
      CHARACTER*80 TMPLIN

      K = INDEX(ILINE(JLINE:),CHAR(0))
      IF (K.EQ.0) THEN
         JLINE = 504
         L = 0
         LINE = ' '

      ELSE
         L = K - 1
         IF (L.GT.LEN(LINE)) THEN
            L = LEN(LINE)
            LINE = ILINE(JLINE:JLINE+L-1)
            TMPLIN = 'Remainder truncated:'//LINE
            CALL LXERR(TMPLIN,1)

         ELSE IF (L.GT.0) THEN
            LINE = ILINE(JLINE:JLINE+K-1)
         END IF

         JLINE = JLINE + K
      END IF

      RETURN

      END
