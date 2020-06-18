C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Log: grexit.f,v $
C Revision 1.2  2009/03/25 12:36:44  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:02:18  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:51:35  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE GREXIT
C=======================================================================

C   --*** GREXIT *** (GRPLIB) Cleanup graphics for exit (PLT)
C   --   Written by Amy Gilkey - revised 12/18/87
C   --
C   --GREXIT terminates graphic output.
C   --
C   --Common Variables:
C   --   Uses DEVOK of /GRPCOM/

C   --Routines Called:
C   --   GRSNAP - (GRPLIB) Initialize frame snap
C   --   PLTEND - (PLTLIB) Terminate graphics

      PARAMETER (KDVDI=10000)

      COMMON /GRPCOC/ DEVNAM(2), DEVCOD(2)
      CHARACTER*3 DEVNAM
      CHARACTER*8 DEVCOD
      COMMON /GRPCOM/ ICURDV, ISHARD, DEVOK(2), TALKOK(2),
     &   NSNAP(2), IFONT(2), SOFTCH(2), AUTOPL(2),
     &   MAXCOL(2), NUMCOL(0:1,2), MAPALT(2), MAPUSE(2)
      LOGICAL ISHARD, DEVOK, TALKOK, SOFTCH, AUTOPL

      IF (DEVOK(1)) CALL GRSNAP ('EXIT', 1)
      IF (DEVOK(2)) CALL GRSNAP ('EXIT', 2)

C   --Turn on both graphics devices to ensure proper monitoring
      IDEV = 0
      IF (DEVOK(1)) IDEV = IDEV + 1
      IF (DEVOK(2)) IDEV = IDEV + 2
      CALL VDESCP (KDVDI + IDEV, 0, 0)

      CALL PLTEND

      RETURN
      END
