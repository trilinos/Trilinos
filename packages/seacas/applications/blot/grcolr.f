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
      SUBROUTINE GRCOLR (INDX)
C=======================================================================

C   --*** GRCOLR *** (GRPLIB) Set color (PLT)
C   --   Written by Amy Gilkey - revised 04/11/88
C   --
C   --GRCOLR sets the color of lines depending on the passed index.
C   --The colors are chosen consecutively up to the maximum set,
C   --skipping black and white, and wrapping around if necessary.
C   --The colors of lines and symbols drawn by PLTGRH are also set.
C   --
C   --Parameters:
C   --   INDX - IN - the color index:
C   --      -1 = black
C   --      =0 = white
C   --      +n = color number n
C   --
C   --Common Variables:
C   --   Uses ICURDV, MAXCOL, NUMCOL, MAPUSE of /GRPCOM/

C   --Routines Called:
C   --   PLTICL - (PLTLIB) Get color index from name
C   --   PLTSTD - (PLTLIB) Set device parameter
C   --      1 = (KCOLOR) set color
C   --   PLTSTG - (PLTLIB) Set graph parameter
C   --      6 = (KCOLIN) set line color for PLTGRH lines
C   --      44 = (KCOSYM) set symbol color for PLTGRH lines

      PARAMETER (KCOLOR=1)
      PARAMETER (KCOLIN=6, KCOSYM=44)
      include 'params.blk'
      include 'cmap-lst.blk'

      COMMON /GRPCOC/ DEVNAM(2), DEVCOD(2)
      CHARACTER*3 DEVNAM
      CHARACTER*8 DEVCOD
      COMMON /GRPCOM/ ICURDV, ISHARD, DEVOK(2), TALKOK(2),
     &   NSNAP(2), IFONT(2), SOFTCH(2), AUTOPL(2),
     &   MAXCOL(2), NUMCOL(0:1,2), MAPALT(2), MAPUSE(2)
      LOGICAL ISHARD, DEVOK, TALKOK, SOFTCH, AUTOPL

      include 'grcol.blk'

      INTEGER INDX

      LOGICAL LDUM, PLTSTG, PLTSTD, PLTICL
      CHARACTER*8 BLCOLR, WHCOLR
      SAVE BLCOLR, WHCOLR

      DATA BLCOLR, WHCOLR / 'BLACK   ', 'WHITE   ' /
      LSTDEV = 0

      IF (ICURDV .NE. LSTDEV) THEN
         LSTDEV = ICURDV
         LSTMAP = 0

C      --Reset last color
         LSTCOL = -999

C      --Set device black and white
         IF (MAXCOL(ICURDV) .GT. 0) THEN
            LDUM = PLTICL (WHCOLR, RWHITE)
            IWHITE = INT(RWHITE)
            LDUM = PLTICL (BLCOLR, RBLACK)
            IBLACK = INT(RBLACK)
         ELSE
            IBLACK = 0
            IWHITE = 1
         END IF

      ELSE IF (MAPUSE(ICURDV) .NE. LSTMAP) THEN

C      --Reset last color
         LSTCOL = -999
         LSTMAP = MAPUSE(ICURDV)
      END IF

      NCOL = NUMCOL(MAPUSE(ICURDV),ICURDV)

      IF ((INDX .GT. 0) .AND. (NCOL .GT. 0)) THEN
         ICOLOR = MOD (INDX-1, NCOL)
         IF (MAPUSE(ICURDV) .EQ. 0) THEN
            IF (ICOLOR .GE. IBLACK) ICOLOR = ICOLOR + 1
            IF (ICOLOR .GE. IWHITE) ICOLOR = ICOLOR + 1
            IF (ICOLOR .GT. IBLACK .AND.
     &             ICOLOR .LT. IWHITE) ICOLOR = INT(COLMAP(ICOLOR))
         ELSE
            ICOLOR = ICOLOR + (6 + 2)
         END IF
      ELSE IF (INDX .EQ. -1) THEN
         ICOLOR = IBLACK
      ELSE
         ICOLOR = IWHITE
      END IF

      IF (ICOLOR .NE. LSTCOL) THEN
         LDUM = PLTSTD (KCOLOR, REAL (ICOLOR))
         LDUM = PLTSTG (KCOLIN, REAL (ICOLOR))
         LDUM = PLTSTG (KCOSYM, REAL (ICOLOR))
         LSTCOL = ICOLOR
      END IF

      RETURN
      END
