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
      SUBROUTINE SPLABN (IPVAR, TIME, NENUM, NAMES,
     &   PLTITL, TXLAB, TYLAB, MAPEL, MAPND)
C=======================================================================

C   --*** SPLABN *** (SPLOT) Get neutral file plot labels
C   --   Written by Amy Gilkey - revised 03/10/86
C   --
C   --SPLABN makes up the plot titles and labels for the neutral file.
C   --
C   --Parameters:
C   --   IPVAR - IN - the /SPVARS/ index of the starting plot variable
C   --   TIME - IN - the plot time
C   --   NENUM - IN - the node/element numbers
C   --   NAMES - IN - the variable names
C   --   PLTITL - OUT - the plot title describing the curves to be
C   --      plotted (e.g. "TIME vs SIGXX at ELEMENT 30")
C   --   TXLAB, TYLAB - OUT - the X and Y axis labels, either the
C   --      user-input labels or the plot variable descriptions
C   --
C   --Common Variables:
C   --   Uses NODVAR, NNENUM of /SELNE/
C   --   Uses ISVID of /SPVARS/
C   --   Uses XLAB, YLAB of /XYLAB/

      include 'params.blk'
      include 'dbnums.blk'
      include 'selne.blk'
      include 'spvars.blk'
      include 'xylab.blk'

      INTEGER NENUM(NNENUM)
      CHARACTER*(*) NAMES(*)
      CHARACTER*(*) PLTITL
      CHARACTER*(*) TXLAB, TYLAB
      INTEGER MAPEL(*), MAPND(*)

      CHARACTER*32 STRNUM
      CHARACTER*32 STRTIM
      CHARACTER*(MXNAME) NAM

C   --Get the plot legend

      NAM = NAMES(ISVID(IPVAR))

      if (nodvar) then
        WRITE (STRNUM, 10000, IOSTAT=IDUM)
     *    MAPND(NENUM(1)), MAPND(NENUM(NNENUM))
      else
        WRITE (STRNUM, 10000, IOSTAT=IDUM)
     *    MAPEL(NENUM(1)), MAPEL(NENUM(NNENUM))
      end if
10000  FORMAT (I9, '..', I9)
      CALL PCKSTR (1, STRNUM)

      CALL NUMSTR (1, 4, TIME, STRTIM, LSTR)

      IF (NODVAR) THEN
         PLTITL = 'DISTANCE vs ' // NAM(:LENSTR(NAM))
     &      // ' NODES ' // STRNUM(:LENSTR(STRNUM))
     &      // ' at TIME ' // STRTIM(:LSTR)
      ELSE
         PLTITL = 'DISTANCE vs ' // NAM(:LENSTR(NAM))
     &      // ' ELEMENTS ' // STRNUM(:LENSTR(STRNUM))
     &      // ' at TIME ' // STRTIM(:LSTR)
      END IF

C   --Get the axis labels

      IF (XLAB .NE. ' ') THEN
         TXLAB = XLAB
      ELSE
         TXLAB = 'DISTANCE'
      END IF

      IF (YLAB .NE. ' ') THEN
         TYLAB = YLAB
      ELSE
         TYLAB = NAM
      END IF

      RETURN
      END
