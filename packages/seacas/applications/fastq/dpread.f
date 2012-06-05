C $Id: dpread.f,v 1.1 1990/11/30 11:06:19 gdsjaar Exp $
C $Log: dpread.f,v $
C Revision 1.1  1990/11/30 11:06:19  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]DPREAD.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE DPREAD (X, Y, BUTTON)
C***********************************************************************
C
C  SUBROUTINE DPREAD = READS INPUT FROM A DIGIPAD DIGITIZING TABLET
C
C***********************************************************************
C
      CHARACTER*1 BUTTON, DUMMY*5
C
C  SWITCH THE TERMINAL TO PASS-THRU MODE <ESC>[5i
C
      DUMMY (1:1) = '+'
      DUMMY (2:2) = CHAR (27)
      DUMMY (3:5) = '[5i'
      WRITE (*, ' (A)')DUMMY
C
C  INPUT THE BUTTON AND X, Y PAIR FROM THE PAD
C
      BUTTON = ' '
      READ (*, 10000, END = 100)BUTTON, IX, IY
C
C  CONVERT THE BUTTON
C
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
C
C CONVERT  (X,  Y) LOCATION
C
      X = IX
      Y = IY
C
  100 CONTINUE
C
C  SWITCH THE TERMINAL OUT OF PASS-THRU MODE <ESC>[4i
C
      WRITE (*, ' (A)')' '//CHAR (27)//'[4i'
      RETURN
C
10000 FORMAT (A1, I5, 1X, I5)
C
      END
