C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE GRINIT (DBORD0, CHLSIZ)
C=======================================================================

C   --*** GRINIT *** (GRPLIB) Initialize graphics (PLT)
C   --   Written by Amy Gilkey - revised 04/27/88
C   --
C   --GRINIT initializes the graphics and sets graphics options and
C   --parameters common to the entire run.
C   --
C   --Alphanumeric mode is set upon exit from this routine.
C   --
C   --Parameters:
C   --   DBORD0 - OUT - the display / label area boundary (device units)
C   --      (1=left, 2=right, 3=bottom, 4=top)
C   --   CHLSIZ - OUT - the character line size (device units)
C   --
C   --Common Variables:
C   --   Sets ICURDV, DEVOK, DEVNAM, DEVCOD, TALKOK, NSNAP, IFONT, SOFTCH,
C   --      AUTOPL, MAXCOL, NUMCOL, MAPALT of /GRPCOM/

C   --Routines Called:
C   --   EXPARM - (SUPES) Get system-dependent parameters
C   --   PLTINT - (PLTLIB) Initialize Graphics Status Area (GSA)
C   --   PLTIQD - (PLTLIB) Obtain device information
C   --   PLTSTT - (PLTLIB) Set text parameter
C   --      1 = (KHCHSZ) hardware character size
C   --      2 = (KSCHSZ) software character size
C   --      6 = (KCHLSZ) vertical line size
C   --   GRSDEV - (GRPLIB) Select device

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4)

      PARAMETER (KDVDI=10000)
      PARAMETER (KHCHSZ=1, KSCHSZ=2, KCHLSZ=6)

      include 'grpcom.blk'

      REAL DBORD0(KTOP)
      REAL CHLSIZ

      LOGICAL LDUM, PLTSTT
      LOGICAL BATCH
      LOGICAL HRDCPY
      REAL VDINFO(23)

      character*2048 scratch
      character*32  device

      PARAMETER (MXDCOD = 12)

      CHARACTER*3 NDCOD(MXDCOD)
      CHARACTER*8 DCOD(MXDCOD)
      SAVE NDCOD, DCOD
C      --NDCOD - the device name for special devices
C      --DCOD - the device code for special devices

      DATA NDCOD /   '16B',      '16C',
     &   '35A',      '35B',      '35C',      '3MB',      '3MC',
     &   '24L',      '48L',      'BSQ',      'CSQ',
     &   'TK4' /
      DATA DCOD  /   'DICOMED ', 'DICOMED ',
     &   'DICOMED ', 'DICOMED ', 'DICOMED ', 'DICOMED ', 'DICOMED ',
     &   'DICOMED ', 'DICOMED ', 'DICOMED ', 'DICOMED ',
     &   'WAIT    ' /

C   --Set up graph borders, etc

      CHLSIZ = .014
      DBORD0(KLFT) = 0.00 + 0.0001
      DBORD0(KRGT) = 1.00 - 0.0001
      DBORD0(KBOT) = 0.00 + 0.0001
      DBORD0(KTOP) = 0.75 - 0.0001
C   --Open unit 6 to tty (for PLT)

c      OPEN (UNIT=6, FILE='tty', ERR=100)
c  100 CONTINUE

C   --Get graphic output devices ready

C ... Get graphics device name from executable (follows . in executable name)
      CALL GET_ARGUMENT(0,scratch, lfil)
      last = indexr(scratch, '_')
      if (last .gt. 2) then
        device = scratch(last+1:lfil)
        if (device(:lenstr(device)) .eq. 'dual') then
          devnam(1) = 'x11'
          devnam(2) = 'met'
        else if (device(:lenstr(device)) .eq. 'xcps') then
          devnam(1) = 'x11'
          devnam(2) = 'cps'
        else
          devnam(1) = device(:lenstr(device))
          devnam(2) = ' '
        end if
      else
        last = indexr(scratch, '.')
        if (last .gt. 2) then
          device = scratch(last+1:lfil)
          if (device(:lenstr(device)) .eq. 'dual') then
            devnam(1) = 'x11'
            devnam(2) = 'met'
          else if (device(:lenstr(device)) .eq. 'xcps') then
            devnam(1) = 'x11'
            devnam(2) = 'cps'
          else
            devnam(1) = device(:lenstr(device))
            devnam(2) = ' '
          end if
        else
          call PRTERR('ERROR',
     *      'Could not determine graphics device type.')
          call PRTERR('CMDSPEC',scratch(:lfil))
C ... Assume single device, x11
          devnam(1) = 'x11'
          devnam(2) = ' '
        end if
      end if

      CALL VDIQES (KDVDI, ION)
      IF (ION .EQ. 1) THEN
         DO 110 IDEV = 1, 2
            CALL VDIQES (KDVDI + IDEV, ION)
            DEVOK(IDEV) = (ION .EQ. 1)
            IF (.NOT. DEVOK(IDEV)) THEN
              DEVNAM(IDEV) = ' '
            END IF
  110    CONTINUE
      ELSE
         DO 120 IDEV = 1, 2
            DEVOK(IDEV) = (DEVNAM(IDEV) .NE. ' ')
  120    CONTINUE
         IF (.NOT. (DEVOK(1) .OR. DEVOK(2))) THEN
            DEVOK(1) = .TRUE.
            DEVNAM(1) = '***'
         END IF
      END IF

C ... Temporary kludge until dual device available GDS
      IF (DEVNAM(2) .eq. '***' ) then
         DEVOK(2) = .FALSE.
         DEVNAM(2) = ' '
      end if

      IF (.NOT. (DEVOK(1) .OR. DEVOK(2))) THEN
         CALL PRTERR ('FATAL', 'No device is active')
         STOP
      END IF

C   --Initialize the attached camera device (before the graphics device)

      NSNAP(1) = -999
      NSNAP(2) = -999

C   --Initialize graphics devices

      IDEV = 0
      IF (DEVOK(1)) IDEV = IDEV + 1
      IF (DEVOK(2)) IDEV = IDEV + 2
      CALL VDESCP (KDVDI + IDEV, 0, 0)
      CALL PLTINT

C   --Set device dependent parameters, for both devices

      DO 140 IDEV = 1, 2
         IF (DEVOK(IDEV)) THEN
            CALL VDESCP (KDVDI + IDEV, 0, 0)

            CALL PLTIQD (VDINFO)

            DEVCOD(IDEV) = ' '
            DO 130 I = 1, MXDCOD
               IF (NDCOD(I) .EQ. DEVNAM(IDEV)) DEVCOD(IDEV) = DCOD(I)
  130       CONTINUE
            HRDCPY = (NINT(VDINFO(1)) .EQ. 0)
            TALKOK(IDEV) = (.NOT. HRDCPY) .AND. (.NOT. BATCH ())
C         --Old: high resolution font if (1.75*VDINFO(15)) .GT. 1600
            CALL GRSNAP ('QUERY', IDEV)
            IFONT(IDEV) = 1
            SOFTCH(IDEV) = HRDCPY .OR. (NINT(VDINFO(7)) .EQ. 0)
            MAXCOL(IDEV) = MAX (0, NINT(VDINFO(4))-2)
            NUMCOL(0,IDEV) = MIN (MAXCOL(IDEV), 6)
            NUMCOL(1,IDEV) = 0
            MAPALT(IDEV) = 0
            MAPUSE(IDEV) = 0
            AUTOPL(IDEV) = .NOT. TALKOK(IDEV)

         ELSE
            HRDCPY = .FALSE.
            TALKOK(IDEV) = .FALSE.
            NSNAP(IDEV) = -999
            IFONT(IDEV) = 1
            SOFTCH(IDEV) = .FALSE.
            MAXCOL(IDEV) = 0
            NUMCOL(0,IDEV) = 0
            NUMCOL(1,IDEV) = 0
            MAPALT(IDEV) = 0
            MAPUSE(IDEV) = 0
            AUTOPL(IDEV) = .NOT. TALKOK(IDEV)
         END IF
  140 CONTINUE

C   --Set current device parameters

      ICURDV = 0
      CALL GRSDEV (0)

C   --Set character sizes

      RAT = 4.0 / 3.0
      VCS = CHLSIZ / RAT
      LDUM = PLTSTT (KSCHSZ, VCS)
      LDUM = PLTSTT (KHCHSZ, VCS)
      LDUM = PLTSTT (KCHLSZ, RAT)

      RETURN
      END
