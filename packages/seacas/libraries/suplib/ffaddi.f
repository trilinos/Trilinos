C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE FFADDI (IVAL, LINE)
C=======================================================================

C   --*** FFADDI *** (FFLIB) Add integer to line
C   --   Written by Amy Gilkey - revised 11/16/87
C   --
C   --FFADDI adds an integer (as a character string) to a line.
C   --
C   --Parameters:
C   --   IVAL - IN - the integer to add
C   --   LINE - IN/OUT - the line being built

      INTEGER IVAL
      CHARACTER*(*) LINE

      CHARACTER*10 STR

      IF (LINE .EQ. ' ') THEN
         I = 0
      ELSE
         I = LENSTR (LINE) + 1
         IF (I .LE. LEN (LINE)) LINE(I:I) = ' '
      END IF
      IF (I .LT. LEN (LINE)) THEN
         CALL INTSTR (1, 0, IVAL, STR, L)
         LINE(I+1:) = STR
      END IF

      RETURN
      END
