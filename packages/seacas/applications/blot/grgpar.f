C Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software.
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
C     * Neither the name of Sandia Corporation nor the names of its
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

C $Log: grgpar.f,v $
C Revision 1.2  2009/03/25 12:36:44  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:02:26  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:51:39  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE GRGPAR (PARTYP, INDEV, IPARMS, IPSTR)
C=======================================================================

C   --*** GRGPAR *** (GRPLIB) Get graphics parameters
C   --   Written by Amy Gilkey - revised 04/11/88
C   --
C   --GRGPAR gets information specified by PARTYP about a graphics device.
C   --   DEVICE - select current device
C   --   SNAP - set the number of frame snaps
C   --   FONT - select font to use
C   --   SOFTCHAR - select software or hardware characters
C   --   COLOR - select number of standard colors to use
C   --   SPECTRUM - select number of spectrum colors to use
C   --   AUTO - set automatic plotting
C   --
C   --Parameters:
C   --   INDEV - IN - the device for which information is requested;
C   --      0 for current device
C   --   PARTYP - IN - the parameter type (as above)
C   --   IPARMS - OUT - the parameter values (dependent on PARTYP)
C   --      for DEVICE - true iff device exists
C   --      for SNAP - the number of frames to snap
C   --      for FONT - the font to use
C   --      for SOFTCHAR - true iff software characters
C   --      for COLOR - the number of standard colors to use and the
C   --         maximum number of colors
C   --      for SPECTRUM - the number of alternate colors to use and the
C   --         maximum number of colors and the map type
C   --      for AUTO - true iff automatic plotting
C   --   IPARMS - OUT - the parameter string (dependent on PARTYP)
C   --      for DEVICE - device name
C   --      for SNAP - not set
C   --      for FONT - the font type
C   --      for SOFTCHAR - "software" or "hardware"
C   --      for COLOR, SPECTRUM - the color map type "standard" or "spectrum"
C   --      for AUTO - "automatic" or "user-directed"
C   --
C   --Common Variables:
C   --   Uses ICURDV, DEVOK, DEVNAM, NSNAP, IFONT, SOFTCH, AUTOPL,
C   --      MAXCOL, NUMCOL, MAPALT of /GRPCOM/

      COMMON /GRPCOC/ DEVNAM(2), DEVCOD(2)
      CHARACTER*3 DEVNAM
      CHARACTER*8 DEVCOD
      COMMON /GRPCOM/ ICURDV, ISHARD, DEVOK(2), TALKOK(2),
     &   NSNAP(2), IFONT(2), SOFTCH(2), AUTOPL(2),
     &   MAXCOL(2), NUMCOL(0:1,2), MAPALT(2), MAPUSE(2)
      LOGICAL ISHARD, DEVOK, TALKOK, SOFTCH, AUTOPL

      CHARACTER*(*) PARTYP
      INTEGER INDEV
      INTEGER IPARMS(*)
      CHARACTER*(*) IPSTR

      IF ((INDEV .NE. 1) .AND. (INDEV .NE. 2)) THEN
         IDEV = ICURDV
      ELSE
         IDEV = INDEV
      END IF

      IF (PARTYP .EQ. 'DEVICE') THEN
         IF (DEVOK(IDEV)) THEN
            CALL CPYLOG (1, .TRUE., IPARMS(1))
            IPSTR = DEVNAM(IDEV)
         ELSE
            CALL CPYLOG (1, .FALSE., IPARMS(1))
            IPSTR = ' '
         END IF

      ELSE IF (PARTYP .EQ. 'SNAP') THEN
         IPARMS(1) = NSNAP(IDEV)

      ELSE IF (PARTYP .EQ. 'FONT') THEN
         IPARMS(1) = IFONT(IDEV)
         IF (IFONT(IDEV) .EQ. 1) THEN
            IPSTR = 'stick'
         ELSE IF (IFONT(IDEV) .EQ. 2) THEN
            IPSTR = 'sanserif'
         ELSE IF (IFONT(IDEV) .EQ. 3) THEN
            IPSTR = 'Roman'
         ELSE
            IPSTR = ' '
         END IF

      ELSE IF (PARTYP .EQ. 'SOFTCHAR') THEN
         CALL CPYLOG (1, SOFTCH(IDEV), IPARMS(1))
         IF (SOFTCH(IDEV)) THEN
            IPSTR = 'software'
         ELSE
            IPSTR = 'hardware'
         END IF

      ELSE IF (PARTYP .EQ. 'COLOR') THEN
         IPARMS(1) = NUMCOL(0,IDEV)
         IPARMS(2) = MIN (6, MAXCOL(IDEV))
         IPARMS(3) = 0
         IPSTR = 'standard'

      ELSE IF (PARTYP .EQ. 'SPECTRUM') THEN
         IPARMS(3) = MAPALT(IDEV)
         IF (MAPALT(IDEV) .EQ. 0) THEN
            IPARMS(1) = NUMCOL(0,IDEV)
            IPARMS(2) = MIN (6, MAXCOL(IDEV))
            IPSTR = 'standard'
         ELSE IF (MAPALT(IDEV) .EQ. 1) THEN
            IPARMS(1) = NUMCOL(1,IDEV)
            IPARMS(2) = MAXCOL(IDEV) - 6
            IPSTR = 'spectrum'
         ELSE
            IPARMS(1) = 0
            IPARMS(2) = MAXCOL(IDEV)
            IPSTR = 'invalid'
         END IF

      ELSE IF (PARTYP .EQ. 'AUTO') THEN
         CALL CPYLOG (1, AUTOPL(IDEV), IPARMS(1))
         IF (AUTOPL(IDEV)) THEN
            IPSTR = 'automatic'
         ELSE
            IPSTR = 'user-directed'
         END IF
      END IF

      RETURN
      END
