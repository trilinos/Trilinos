C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE DPREAD (X, Y, BUTTON)
C***********************************************************************

C  SUBROUTINE DPREAD = READS INPUT FROM A DIGIPAD DIGITIZING TABLET

C***********************************************************************

      CHARACTER*1 BUTTON, DUMMY*5

C  SWITCH THE TERMINAL TO PASSTHROUGH MODE <ESC>[5i

      DUMMY (1:1) = '+'
      DUMMY (2:2) = CHAR (27)
      DUMMY (3:5) = '[5i'
      WRITE (*, ' (A)')DUMMY

C  INPUT THE BUTTON AND X, Y PAIR FROM THE PAD

      BUTTON = ' '
      READ (*, 10000, END = 100)BUTTON, IX, IY

C  CONVERT THE BUTTON

      IF (BUTTON .EQ. ':') THEN
         BUTTON = 'A'
      ELSEIF (BUTTON .EQ. ';') THEN
         BUTTON = 'B'
      ELSEIF (BUTTON .EQ. '<') THEN
         BUTTON = 'C'
      ELSEIF (BUTTON .EQ. ' = ') THEN
         BUTTON = 'D'
      ELSEIF (BUTTON .EQ. '>') THEN
         BUTTON = 'E'
      ELSEIF (BUTTON .EQ. '?') THEN
         BUTTON = 'F'
      ELSEIF (BUTTON .EQ. ' ') THEN
         BUTTON = 'E'
      END IF

C CONVERT  (X,  Y) LOCATION

      X = IX
      Y = IY

  100 CONTINUE

C  SWITCH THE TERMINAL OUT OF PASSTHROUGH MODE <ESC>[4i

      WRITE (*, ' (A)')' '//CHAR (27)//'[4i'
      RETURN

10000 FORMAT (A1, I5, 1X, I5)

      END
