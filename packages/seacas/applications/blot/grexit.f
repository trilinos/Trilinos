C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

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

      include 'grpcom.blk'

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
