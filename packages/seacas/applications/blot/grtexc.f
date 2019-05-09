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

C $Log: grtexc.f,v $
C Revision 1.2  2009/03/25 12:36:45  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:02:53  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:52:01  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE GRTEXC (DX, DY, STRING)
C=======================================================================

C   --*** GRTEXC *** (GRPLIB) Write centered text (PLT)
C   --   Written by Amy Gilkey - revised 03/22/88
C   --
C   --GRTEXC writes a software or hardware character string centered on
C   --a location.
C   --
C   --Parameters:
C   --   DX, DY - IN - the horizontal and vertical string location
C   --      (in device coordinates)
C   --   STRING - IN - the string to be written, may be truncated
C   --
C   --Common Variables:
C   --   Uses ICURDV, SOFTCH of /GRPCOM/

C   --Routines Called:
C   --   PLTXSL - (PLTLIB) Find the software string length
C   --   PLTXTC1 - (PLTLIB) Display a centered software string
C   --   PLTXTH - (PLTLIB) Display a hardware string
C   --   PLTXTL - (PLTLIB) Find the hardware string length
C   --   LENSTR - (STRLIB) Find string length

      COMMON /GRPCOC/ DEVNAM(2), DEVCOD(2)
      CHARACTER*3 DEVNAM
      CHARACTER*8 DEVCOD
      COMMON /GRPCOM/ ICURDV, ISHARD, DEVOK(2), TALKOK(2),
     &   NSNAP(2), IFONT(2), SOFTCH(2), AUTOPL(2),
     &   MAXCOL(2), NUMCOL(0:1,2), MAPALT(2), MAPUSE(2)
      LOGICAL ISHARD, DEVOK, TALKOK, SOFTCH, AUTOPL

      REAL DX, DY
      CHARACTER*(*) STRING

      LSTR = LENSTR(STRING)
      IF (STRING(LSTR:LSTR) .EQ. ' ') RETURN

      IF (SOFTCH(ICURDV)) THEN
         CALL PLTXSL (STRING(:LSTR), SLEN)
         XLEFT = DX
         IF (XLEFT*2 .LT. SLEN) XLEFT = 0.5*SLEN
         CALL PLTXTC1 (XLEFT, DY, STRING(:LSTR))
      ELSE
         CALL PLTXHL (STRING(:LSTR), SLEN)
         XLEFT = DX - .5*SLEN
         IF (XLEFT .LT. 0.0) XLEFT = 0.0000
         CALL PLTXTH (XLEFT, DY, STRING(:LSTR))
      END IF

      RETURN
      END
