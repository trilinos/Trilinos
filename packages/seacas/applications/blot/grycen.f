C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE GRYCEN (CHLSIZ, TOPLIN, BOTLIN, NUMLIN, NUMOVR)
C=======================================================================

C   --*** GRYCEN *** (GRPLIB) Find center line of text area
C   --   Written by Amy Gilkey - revised 02/20/87
C   --
C   --GRYCEN finds the center of a section of screen and returns the
C   --top and bottom line coordinates for the given number of lines.
C   --
C   --Parameters:
C   --   CHLSIZ - IN - the size of a line of text
C   --   TOPLIN, BOTLIN - IN/OUT - the device coordinates of the bottom of
C   --      the top and bottom lines of text
C   --   NUMLIN - IN/OUT - the number of lines requested, reduced by
C   --      NUMOVR if too long
C   --   NUMOVR - OUT - if the number of lines requested is greater than
C   --      the number of lines allowed, NUMOVR is the number of lines
C   --      deleted

      MAXLIN = INT((TOPLIN - BOTLIN) / CHLSIZ + 1 + 0.25)
      IF (NUMLIN .GT. MAXLIN) THEN
         NUMOVR = NUMLIN - MAXLIN
         NUMLIN = MAXLIN
      ELSE
         NUMOVR = 0
      END IF
      CEN = 0.5 * (TOPLIN + BOTLIN)
      TOPLIN = CEN + 0.5 * (NUMLIN-1) * CHLSIZ
      BOTLIN = TOPLIN - (NUMLIN-1) * CHLSIZ

      RETURN
      END
