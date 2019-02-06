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
      SUBROUTINE TPMAIN (A, NEUTRL, NAMES,
     &   NPTIMS, IPTIMS, TIMES, WHOTIM, BLKCOL,
     &   IDELB, MAPEL, MAPND)
C=======================================================================

C   --*** TPMAIN *** (TPLOT) TPLOT main plot routine
C   --   Written by Amy Gilkey - revised 02/29/88
C   --
C   --TPMAIN reads the plot variables, and displays the time-history or
C   --the variable-versus-variable plots.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   NEUTRL - IN  - the type of neutral file to write.
C   --   NAMES - IN - the variable names
C   --   NPTIMS - IN - the number of selected time steps
C   --   IPTIMS - IN - the selected time steps
C   --   TIMES - IN - the database times
C   --   WHOTIM - IN - true iff whole (versus history) time step
C   --   BLKCOL - IN/OUT - the user selected colors of the element blocks.
C   --                    BLKCOL(0) = 1 if the user defined material
C   --                                colors should be used in mesh plots.
C   --                              = -1 if program selected colors should
C   --                                be used.
C   --                    BLKCOL(i) = the user selected color of element
C   --                               block i:
C   --                                  -2 - no color selected by user.
C   --                                  -1 - black
C   --                                   0 - white
C   --                                   1 - red
C   --                                   2 - green
C   --                                   3 - yellow
C   --                                   4 - blue
C   --                                   5 - cyan
C   --                                   6 - magenta
C   --
C   --Common Variables:
C   --   Uses NTPVAR of /TPVARS/

      include 'params.blk'

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4)
C      --These parameters define the indices of 2D limit arrays

      include 'dbnums.blk'
      include 'tpvars.blk'

      DIMENSION A(*)
      CHARACTER*(*) NAMES(*)
      INTEGER IPTIMS(NPTIMS)
      REAL TIMES(*)
      LOGICAL WHOTIM(*)
      INTEGER BLKCOL(0:NELBLK)
      INTEGER IDELB(*)
      INTEGER MAPEL(*), MAPND(*)

      REAL TIMLIM(2)

C   --Use the selected color table
      CALL GRCOLU ('ALTERNATE')

C   --Reserve storage for plot variables

C   --*NOTE*  PLTVAL(x,NTPVAR+1) holds the time data
C   --        PLTVAL(x,NTPVAR+2) holds the compressed time data
      LPLVAR = NPTIMS * NTPVAR
      L = 0
      IF (TIMPLT) THEN
         L = NPTIMS + NPTIMS
      END IF
      CALL MDRSRV ('PLTVAL', KPLVAL, LPLVAR + L)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 100

C   --Read plot variables

      CALL TPREAD (A, NPTIMS, IPTIMS, TIMES, WHOTIM,
     &   A(KPLVAL+LPLVAR), A(KPLVAL))

      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 100

C   --Compress plot variables

      CALL MDRSRV ('NPTCRV', KNPTS, NTPVAR)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 100

      CALL SQZTPV (NPTIMS, IPTIMS, WHOTIM, A(KNPTS), A(KPLVAL))

C   --Plot

      TIMLIM(1) = TIMES(IPTIMS(1))
      TIMLIM(2) = TIMES(IPTIMS(NPTIMS))

      CALL TPPLOT (NEUTRL, NPTIMS, A(KNPTS), TIMLIM, A(KPLVAL),
     &   NAMES, BLKCOL, MAPEL, MAPND)

      CALL MDDEL ('PLTVAL')
      CALL MDDEL ('NPTCRV')

  100 CONTINUE
      RETURN
      END
