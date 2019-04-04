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
