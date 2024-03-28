C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE GRCOLU (MAP)
C=======================================================================

C   --*** GRCOLU *** (GRPLIB) Set color table to use (PLT)
C   --   Written by Amy Gilkey - revised 04/11/88
C   --
C   --GRCOLU sets the color table to use.
C   --
C   --Parameters:
C   --   MAP - IN - the color table map to use:
C   --      'STANDARD'  = standard color map
C   --      'ALTERNATE' = alternate color map (if defined)
C   --
C   --Common Variables:
C   --   Uses DEVOK, ICURDV, NUMCOL, MAPALT of /GRPCOM/
C   --   Sets MAPUSE of /GRPCOM/

      include 'grpcom.blk'

      CHARACTER*(*) MAP

      IF (MAP .EQ. 'ALTERNATE') THEN
         MAPUSE(ICURDV) = MAPALT(ICURDV)
      ELSE
         MAPUSE(ICURDV) = 0
      END IF

      RETURN
      END
