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

C=======================================================================
      SUBROUTINE TPLABN (IPVAR, TIMLIM, NAMES, PLTITL, TXLAB, TYLAB,
     *  MAPEL, MAPND)
C=======================================================================

C   --*** TPLABN *** (TPLOT) Get neutral file plot labels
C   --   Written by Amy Gilkey - revised 03/23/87
C   --
C   --TPLABN makes up the plot titles and labels for the neutral file.
C   --
C   --Parameters:
C   --   IPVAR - IN - the /TPVARS/ index of the starting plot variable
C   --   TIMLIM - IN - the starting and ending times for a
C   --      variable-versus-variable curve
C   --   NAMES - IN - the variable names
C   --   PLTITL - OUT - the plot title describing the curves to be
C   --      plotted (e.g. "TIME vs SIGXX at ELEMENT 30" or
C   --      "LOAD vs SIGXX at ELEMENT 30 for times 0.000 to 15.000")
C   --   TXLAB, TYLAB - OUT - the X and Y axis labels, either the
C   --      user-input labels or the plot variable descriptions
C   --
C   --Common Variables:
C   --   Uses TIMPLT, ITVID, ITVNE of /TPVARS/
C   --   Uses XLAB, YLAB of /XYLAB/

      include 'params.blk'
      include 'tpvars.blk'
      include 'xylab.blk'

      REAL TIMLIM(2)
      CHARACTER*(*) NAMES(*)
      CHARACTER*(*) PLTITL
      CHARACTER*(*) TXLAB, TYLAB
      INTEGER MAPEL(*), MAPND(*)

      CHARACTER*(1024) PV1, PV2
      CHARACTER*20 RSTR(2)

C   --Get the plot legend

      N = IPVAR

      IF (TIMPLT) THEN
         PV1 = 'TIME'
         CALL TPLABV (-1, ITVID(N), NAMES(ITVID(N)), ITVNE(N), PV2,
     *     MAPEL, MAPND)
         PLTITL = PV1(:LENSTR(PV1)) // ' vs ' // PV2(:LENSTR(PV2))
         write (*,*) pltitl(:lenstr(pltitl))
      ELSE
         CALL TPLABV (-1, ITVID(N), NAMES(ITVID(N)), ITVNE(N), PV1,
     *    MAPEL, MAPND)
         N = N + 1
         CALL TPLABV (-1, ITVID(N), NAMES(ITVID(N)), ITVNE(N), PV2,
     *     MAPEL, MAPND)
         CALL NUMSTR (2, 4, TIMLIM, RSTR, LSTR)
         PLTITL = PV1(:LENSTR(PV1)) // ' vs ' // PV2(:LENSTR(PV2))
     &      // ' for times ' // RSTR(1)(:LENSTR(RSTR(1)))
     &      // ' to ' // RSTR(2)(:LSTR)
      END IF

C   --Get the axis labels

      IF (XLAB .NE. ' ') THEN
         TXLAB = XLAB
      ELSE
         TXLAB = PV1
      END IF

      IF (YLAB .NE. ' ') THEN
         TYLAB = YLAB
      ELSE
         TYLAB = PV2
      END IF

      RETURN
      END
