C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Log: grcolu.f,v $
C Revision 1.2  2009/03/25 12:36:44  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:02:16  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:51:33  gdsjaar
c Added RCS Id and Log to all files
c
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

      COMMON /GRPCOC/ DEVNAM(2), DEVCOD(2)
      CHARACTER*3 DEVNAM
      CHARACTER*8 DEVCOD
      COMMON /GRPCOM/ ICURDV, ISHARD, DEVOK(2), TALKOK(2),
     &   NSNAP(2), IFONT(2), SOFTCH(2), AUTOPL(2),
     &   MAXCOL(2), NUMCOL(0:1,2), MAPALT(2), MAPUSE(2)
      LOGICAL ISHARD, DEVOK, TALKOK, SOFTCH, AUTOPL

      CHARACTER*(*) MAP

      IF (MAP .EQ. 'ALTERNATE') THEN
         MAPUSE(ICURDV) = MAPALT(ICURDV)
      ELSE
         MAPUSE(ICURDV) = 0
      END IF

      RETURN
      END
