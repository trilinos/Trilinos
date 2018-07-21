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

C $Log: grsdev.f,v $
C Revision 1.2  2009/03/25 12:36:44  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:02:44  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:51:52  gdsjaar
c Added RCS Id and Log to all files
c
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

      COMMON /GRPCOC/ DEVNAM(2), DEVCOD(2)
      CHARACTER*3 DEVNAM
      CHARACTER*8 DEVCOD
      COMMON /GRPCOM/ ICURDV, ISHARD, DEVOK(2), TALKOK(2),
     &   NSNAP(2), IFONT(2), SOFTCH(2), AUTOPL(2),
     &   MAXCOL(2), NUMCOL(0:1,2), MAPALT(2), MAPUSE(2)
      LOGICAL ISHARD, DEVOK, TALKOK, SOFTCH, AUTOPL

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
