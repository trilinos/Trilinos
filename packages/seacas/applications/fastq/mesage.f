C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE MESAGE (PROMPT)
C***********************************************************************

C  SUBROUTINE MESAGE = PRINTS A MESSAGE ONTO THE SCREEN

C***********************************************************************

      CHARACTER * (*) PROMPT

      IF (PROMPT .EQ. ' ') THEN
         WRITE (*, 10000)
      ELSE
         WRITE (*, 10010)PROMPT
      ENDIF
      RETURN

10000 FORMAT ( / )
10010 FORMAT (' ', A)
      END
