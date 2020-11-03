C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE FFADDO (ISON, LINE)
C=======================================================================

C   --*** FFADDO *** (FFLIB) Add ON/OFF to line
C   --   Written by Amy Gilkey - revised 11/16/87
C   --
C   --FFADDO adds either ON or OFF to a line.
C   --
C   --Parameters:
C   --   ISON - IN - the character string to add (ON if true, OFF if false)
C   --   LINE - IN/OUT - the line being built

      LOGICAL ISON
      CHARACTER*(*) LINE

      IF (LINE .EQ. ' ') THEN
         I = 0
      ELSE
         I = LENSTR (LINE) + 1
         IF (I .LE. LEN (LINE)) LINE(I:I) = ' '
      END IF
      IF (I .LT. LEN (LINE)) THEN
         IF (ISON) THEN
            LINE(I+1:) = 'ON'
         ELSE
            LINE(I+1:) = 'OFF'
         END IF
      END IF

      RETURN
      END
