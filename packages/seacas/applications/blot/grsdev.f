C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE GRSDEV (INDEV)
C=======================================================================

C   --*** GRSDEV *** (GRPLIB) Select device
C   --   Written by Amy Gilkey, revised 08/24/87
C   --
C   --GRSDEV selects a device and changes all the device parameters.
C   --
C   --Parameters:
C   --   INDEV - IN - the device to be selected
C   --
C   --Common Variables:
C   --   Uses ICURDV, DEVOK, IFONT, NUMCOL of /GRPCOM/
C   --   Sets ICURDV of /GRPCOM/

C   --Routines Called:
C   --   GRCOLT - (GRPLIB) Set color table
C   --   GRFONT - (GRPLIB) Set font

      PARAMETER (KDVDI=10000)

      include 'grpcom.blk'

C   --If the device number is not given, choose the first available device
      IF ((INDEV .NE. 1) .AND. (INDEV .NE. 2)) THEN
         IF (DEVOK(1)) THEN
            IDEV = 1
         ELSE
            IDEV = 2
         END IF
      ELSE
         IDEV = INDEV
      END IF

C   --Skip if invalid parameter
      IF (.NOT. DEVOK(IDEV)) GOTO 100

C   --Skip if device already selected
      IF (IDEV .EQ. ICURDV) GOTO 100

C   --Turn off old device and turn on new device
      CALL VDESCP (KDVDI + IDEV, 0, 0)

      ICURDV = IDEV

C   --Set color table
      CALL GRCOLT

C   --Set font
      CALL GRFONT (IFONT(ICURDV))

C   --Set number of frames to snap
      CALL GRSNAP ('INIT', ICURDV)

C   --Set line widths
      CALL GRLWID

C   --Reset the single hardcopy flag if terminal device selected
      IF (ICURDV .EQ. 1) ISHARD = .FALSE.

  100 CONTINUE
      RETURN
      END
