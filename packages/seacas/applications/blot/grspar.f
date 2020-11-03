C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE GRSPAR (PARTYP, INDEV, IPARM, ERRMSG)
C=======================================================================

C   --*** GRSPAR *** (GRPLIB) Set graphics parameter
C   --   Written by Amy Gilkey - revised 04/11/88
C   --
C   --GRSPAR sets a graphics parameter specified by PARTYP:
C   --   DEVICE - select current device
C   --   SNAP - set the number of frame snaps
C   --   FONT - select font to use
C   --   SOFTCHAR - select software or hardware characters
C   --   COLOR - select number of colors to use (normal color map)
C   --   SPECTRUM - select spectrum color map
C   --   AUTO - set automatic plotting
C   --
C   --Parameters:
C   --   PARTYP - IN - the parameter type (as above)
C   --   INDEV - IN - the device to be set; 0 for current device;
C   --      <0 for both devices
C   --      for DEVICE - 0 for terminal if defined, else hardcopy
C   --   IPARM - IN - the parameter value (dependent on PARTYP)
C   --      for SNAP - the number of frames to snap
C   --      for FONT - the font to use (as in /GRPCOM/)
C   --      for SOFTCHAR - true iff software characters
C   --      for COLOR - the number of colors to use
C   --      for SPECTRUM - the number of colors to use
C   --      for AUTO - true iff automatic plotting
C   --   ERRMSG - OUT - the returned error message; ' ' if no error
C   --
C   --Common Variables:
C   --   Uses ICURDV, DEVOK, TALKOK, MAXCOL of /GRPCOM/
C   --   Sets IFONT, SOFTCH, AUTOPL, NUMCOL, MAPALT of /GRPCOM/

C   --Routines Called:
C   --   GRSNAP - (GRPLIB) Initialize frame snap
C   --   GRFONT - (GRPLIB) Set font
C   --   GRCOLT - (GRPLIB) Set color table

      COMMON /GRPCOC/ DEVNAM(2), DEVCOD(2)
      CHARACTER*3 DEVNAM
      CHARACTER*8 DEVCOD
      COMMON /GRPCOM/ ICURDV, ISHARD, DEVOK(2), TALKOK(2),
     &   NSNAP(2), IFONT(2), SOFTCH(2), AUTOPL(2),
     &   MAXCOL(2), NUMCOL(0:1,2), MAPALT(2), MAPUSE(2)
      LOGICAL ISHARD, DEVOK, TALKOK, SOFTCH, AUTOPL

      CHARACTER*(*) PARTYP
      CHARACTER*(*) ERRMSG

      LOGICAL ISON

      ERRMSG = ' '

      IF (PARTYP .EQ. 'DEVICE') THEN
         IF (INDEV .EQ. 0) THEN
            ISDEV = 1
            IF (.NOT. DEVOK(ISDEV)) ISDEV = 2
         ELSE
            ISDEV = INDEV
         END IF

         IF ((ISDEV .NE. 1) .AND. (ISDEV .NE. 2)) THEN
            ERRMSG = 'Invalid device number'
            GOTO 170
         END IF
         IF (.NOT. DEVOK(ISDEV)) THEN
            ERRMSG = 'Device not defined'
            GOTO 170
         END IF

         CALL GRSDEV (ISDEV)

      ELSE
         IF (INDEV .LT. 0) THEN
            ISDEV = 1
            IF (.NOT. DEVOK(ISDEV)) ISDEV = 2
            IEDEV = 2
            IF (.NOT. DEVOK(IEDEV)) IEDEV = 1
         ELSE IF (INDEV .EQ. 0) THEN
            ISDEV = ICURDV
            IEDEV = ICURDV
         ELSE
            ISDEV = INDEV
            IEDEV = INDEV
         END IF

         IF ((ISDEV .NE. 1) .AND. (ISDEV .NE. 2)) THEN
            ERRMSG = 'Invalid device number'
            GOTO 170
         END IF
         IF (.NOT. DEVOK(ISDEV)) THEN
            ERRMSG = 'Device not defined'
            GOTO 170
         END IF

         IF (PARTYP .EQ. 'SNAP') THEN
            IF ((NSNAP(ISDEV) .LT. 0) .AND. (NSNAP(IEDEV) .LT. 0)) THEN
               ERRMSG = 'Device cannot be snapped'
               GOTO 170
            END IF
            N = MAX (0, IPARM)
            DO 100 IDEV = ISDEV, IEDEV
               IF (NSNAP(IDEV) .GE. 0) THEN
                  NSNAP(IDEV) = N
                  IF (ICURDV .EQ. IDEV) CALL GRSNAP ('INIT', IDEV)
               END IF
  100       CONTINUE

         ELSE IF (PARTYP .EQ. 'FONT') THEN
            IF (IPARM .LE. 0) THEN
               IFNT = 1
            ELSE
               IFNT = IPARM
            END IF
            IF ((IFNT .LT. 1) .OR. (IFNT .GT. 3)) THEN
               ERRMSG = 'Invalid font type'
               GOTO 170
            END IF
            DO 110 IDEV = ISDEV, IEDEV
               IFONT(IDEV) = IFNT
               IF (ICURDV .EQ. IDEV) CALL GRFONT (IFNT)
  110       CONTINUE

         ELSE IF (PARTYP .EQ. 'SOFTCHAR') THEN
            CALL CPYINT (1, IPARM, ISON)
            DO 120 IDEV = ISDEV, IEDEV
               SOFTCH(IDEV) = ISON
  120       CONTINUE

         ELSE IF (PARTYP .EQ. 'COLOR') THEN
            NCOL = IPARM
            IF ((NCOL .GT. MAXCOL(ISDEV)) .AND.
     &         (NCOL .GT. MAXCOL(IEDEV))) THEN
               ERRMSG = 'Number of colors is greater than the maximum'
            END IF
            IF (NCOL .GT. 6) THEN
               ERRMSG = 'Number of standard colors is limited to 6'
               NCOL = 6
            END IF
            DO 130 IDEV = ISDEV, IEDEV
               IF (NCOL .GT. 0)
     &            NUMCOL(0,IDEV) = MIN (NCOL, MAXCOL(IDEV))
               IF (ICURDV .EQ. IDEV) CALL GRCOLT
  130       CONTINUE

         ELSE IF (PARTYP .EQ. 'SPECTRUM') THEN
            NCOL = IPARM
            IF (NCOL .LE. 0) THEN
               DO 140 IDEV = ISDEV, IEDEV
                  IF (MAPALT(IDEV) .GT. 0) THEN
                     MAPALT(IDEV) = 0
                     NUMCOL(1,IDEV) = 0
                     IF (ICURDV .EQ. IDEV) CALL GRCOLT
                  END IF
  140          CONTINUE
            ELSE
               IF ((NCOL .GT. MAXCOL(ISDEV)-6) .AND.
     &            (NCOL .GT. MAXCOL(IEDEV)-6)) THEN
                  ERRMSG =
     &               'Number of colors is greater than the maximum'
               END IF
               DO 150 IDEV = ISDEV, IEDEV
                  IF (MAXCOL(IDEV) .GT. 6) THEN
                     MAPALT(IDEV) = 1
                     NUMCOL(1,IDEV) = MIN (NCOL, MAXCOL(IDEV)-6)
                     IF (ICURDV .EQ. IDEV) CALL GRCOLT
                  END IF
  150          CONTINUE
            END IF

         ELSE IF (PARTYP .EQ. 'AUTO') THEN
            CALL CPYINT (1, IPARM, ISON)
            IF ((.NOT. ISON) .AND.
     &         (.NOT. TALKOK(ISDEV)) .AND. (.NOT. TALKOK(IEDEV))) THEN
               ERRMSG = 'Device cannot be user-directed'
               GOTO 170
            END IF
            DO 160 IDEV = ISDEV, IEDEV
               IF (TALKOK(IDEV)) AUTOPL(IDEV) = ISON
  160       CONTINUE

         ELSE
            ERRMSG = 'Invalid parameter type'
            GOTO 170
         END IF
      END IF

  170 CONTINUE
      RETURN
      END
