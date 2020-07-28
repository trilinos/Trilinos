C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE FFADDR (RVAL, LINE)
C=======================================================================

C   --*** FFADDR *** (FFLIB) Add real to line
C   --   Written by Amy Gilkey - revised 11/16/87
C   --
C   --FFADDR adds a real (as a character string) to a line.
C   --
C   --Parameters:
C   --   RVAL - IN - the real to add
C   --   LINE - IN/OUT - the line being built

      REAL RVAL
      CHARACTER*(*) LINE

      CHARACTER*20 STR

      IF (LINE .EQ. ' ') THEN
         I = 0
      ELSE
         I = LENSTR (LINE) + 1
         IF (I .LE. LEN (LINE)) LINE(I:I) = ' '
      END IF
      IF (I .LT. LEN (LINE)) THEN
         CALL NUMSTR (1, 6, RVAL, STR, L)
         LINE(I+1:) = STR
      END IF

      RETURN
      END
