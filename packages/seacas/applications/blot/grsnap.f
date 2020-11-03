C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE GRSNAP (CMD, INDEV)
C=======================================================================

C   --*** GRSNAP *** (GRPLIB) Perform movie snap operations
C   --   Written by Amy Gilkey, revised 01/25/88
C   --
C   --GRSNAP sets up a device for snapping multiple frames for movies.
C   --This code is very device-dependent.
C   --
C   --At the present time, the only device that may be snapped is a camera
C   --(attached to a terminal) or a Dicomed.
C   --
C   --An attached camera device must be pre-initialized and NSNAP must
C   --be set non-negative (zero ok).
C   --
C   --Parameters:
C   --   CMD - IN - the snap operation:
C   --      QUERY = initialize device code, nsnap < 0 if cannot snap frames
C   --      INIT  = initialize device to snap nsnap frames
C   --      ABORT = abort plot
C   --      START = start plot
C   --      STOP  = end plot and snap frames
C   --      EXIT  = deactivate device
C   --   INDEV - IN - the device; 0 for current device (only QUERY allowed
C   --      for non-current device)
C   --
C   --Common Variables:
C   --   Uses ICURDV, DEVCOD, NSNAP of /GRPCOM/

C   --Routines Called:
C   --   VDESCP - (VDI) Send escape sequence to device

      COMMON /GRPCOC/ DEVNAM(2), DEVCOD(2)
      CHARACTER*3 DEVNAM
      CHARACTER*8 DEVCOD
      COMMON /GRPCOM/ ICURDV, ISHARD, DEVOK(2), TALKOK(2),
     &   NSNAP(2), IFONT(2), SOFTCH(2), AUTOPL(2),
     &   MAXCOL(2), NUMCOL(0:1,2), MAPALT(2), MAPUSE(2)
      LOGICAL ISHARD, DEVOK, TALKOK, SOFTCH, AUTOPL

      CHARACTER*(*) CMD
      INTEGER INDEV

      INTEGER RBUF(2)

      CHARACTER*8 IDSNAP(2)
      SAVE IDSNAP
C      --IDSNAP - an artificial device code related to the device code

      DATA IDSNAP / ' ', ' ' /

      IF ((INDEV .NE. 1) .AND. (INDEV .NE. 2)) THEN
         IDEV = ICURDV
      ELSE
         IDEV = INDEV
      END IF

C   --Initialize device code
      IF (IDSNAP(IDEV) .EQ. ' ') THEN
         IDSNAP(IDEV) = 'NONE'
         IF (NSNAP(IDEV) .GE. 0) THEN
            IDSNAP(IDEV) = 'CAMERA'
         ELSE IF (DEVCOD(IDEV) .EQ. 'DICOMED') THEN
            IDSNAP(IDEV) = DEVCOD(IDEV)
         END IF

C      --Initialize the device
         IF (IDSNAP(IDEV) .EQ. 'CAMERA') THEN
C         --Initialized before the routine was called
            CONTINUE
         ELSE IF (IDSNAP(IDEV) .EQ. 'DICOMED') THEN
            CONTINUE
         END IF
      END IF

C   --Set nsnap = -999 if cannot snap device
      IF (IDSNAP(IDEV) .EQ. 'NONE') NSNAP(IDEV) = -999

      IF (CMD .NE. 'QUERY') THEN

C      --Skip if not current device
         IF (IDEV .NE. ICURDV) GOTO 110

C      --Skip if snaps are not requested on device
         IF (NSNAP(ICURDV) .LE. 0) GOTO 110

         CALL PLTFLU

         IF (CMD .EQ. 'INIT') THEN

C         --Initialize number of snaps

            IF (IDSNAP(IDEV) .EQ. 'CAMERA') THEN
               CONTINUE

            ELSE IF (IDSNAP(IDEV) .EQ. 'DICOMED') THEN
               CONTINUE
            END IF

         ELSE IF (CMD .EQ. 'ABORT') THEN

C         --Abort plot

            IF (IDSNAP(IDEV) .EQ. 'CAMERA') THEN
               CONTINUE

            ELSE IF (IDSNAP(IDEV) .EQ. 'DICOMED') THEN
C            --Close segment and delete all segments
               RBUF(1) = 0
               CALL VDESCP (201, 0, RBUF)
               CALL VDESCP (203, 1, RBUF)
            END IF

         ELSE IF (CMD .EQ. 'START') THEN

C         --Start plot

            IF (IDSNAP(IDEV) .EQ. 'CAMERA') THEN
               CONTINUE

            ELSE IF (IDSNAP(IDEV) .EQ. 'DICOMED') THEN
C            --Open segment
               RBUF(1) = 1
               RBUF(2) = 1
               CALL VDESCP (200, 2, RBUF)
            END IF

         ELSE IF (CMD .EQ. 'STOP') THEN

C         --End plot and snap frames

            IF (IDSNAP(IDEV) .EQ. 'CAMERA') THEN
C            --Snap n frames
               CONTINUE

            ELSE IF (IDSNAP(IDEV) .EQ. 'DICOMED') THEN
C            --Close segment
               RBUF(1) = 0
               CALL VDESCP (201, 0, RBUF)

C            --Snap n-1 frames (plot segment with newpage)
               DO 100 I = 1, NSNAP(ICURDV)-1
                  CALL VDNWPG
  100          CONTINUE

C            --Delete all segments
               RBUF(1) = 0
               CALL VDESCP (203, 1, RBUF)
            END IF

         ELSE IF (CMD .EQ. 'EXIT') THEN

C         --Deactivate device

            IF (IDSNAP(IDEV) .EQ. 'CAMERA') THEN
               CONTINUE

            ELSE IF (IDSNAP(IDEV) .EQ. 'DICOMED') THEN
               CONTINUE
            END IF

         ELSE
            GOTO 110
         END IF
      END IF

  110 CONTINUE
      RETURN
      END
