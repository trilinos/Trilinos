C Copyright(C) 2009-2017 National Technology & Engineering Solutions of
C Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C 
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.
C     * Neither the name of NTESS nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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

      COMMON /GRPCOC/ DEVNAM(2), DEVCOD(2)
      CHARACTER*3 DEVNAM
      CHARACTER*8 DEVCOD
      COMMON /GRPCOM/ ICURDV, ISHARD, DEVOK(2), TALKOK(2),
     &   NSNAP(2), IFONT(2), SOFTCH(2), AUTOPL(2),
     &   MAXCOL(2), NUMCOL(0:1,2), MAPALT(2), MAPUSE(2)
      LOGICAL ISHARD, DEVOK, TALKOK, SOFTCH, AUTOPL
C      --For all device-dependent parameters:
C      --   (1) terminal, (2) hardcopy (file)
C      --ICURDV - the selected device number (1 or 2)
C      --ISHARD - true iff a single hardcopy plot is being done
C      --DEVNAM - the device name
C      --DEVOK - true iff the device is defined
C      --DEVCOD - a code associated with a class of devices:
C      --   DICOMED  = Dicomed
C      --   CAMERA   = Raster Tech with video
C      --   WAIT     = TK4 or any device that requires a wait after graph drawn
C      --            = other
C      --TALKOK - true iff interactive graphics device
C      --NSNAP - the number of frames to snap for the device
C      --IFONT - font (1=stick, 2=sanserif, 3=Roman)
C      --SOFTCH - true iff software characters are to be used
C      --AUTOPL - true iff automatic plotting (no response requested at end)
C      --MAXCOL - the maximum number of colors available on the graphics device
C      --   (excluding black and white)
C      --NUMCOL - the number of colors to be used (excluding black and white)
C      --   (0,x) - on standard rainbow
C      --   (1,x) - on alternate color map
C      --MAPALT - the alternate color map type:
C      --   0 = standard rainbow (no alternate)
C      --   1 = spectrum
C      --MAPUSE - the color map type to use (as in MAPALT)

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
