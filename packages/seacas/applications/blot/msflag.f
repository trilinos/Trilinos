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
      SUBROUTINE MSFLAG (ANYDEF, ANYUND,
     &   DOIXF, DON2B, DOELED, DOELEU, DODEAD, DONPS, DOESS, DOSCAL,
     &   MINMSH, MAXMSH, MAXHID)
C=======================================================================

C   --*** MSFLAG *** (MESH) Set mesh plot flags
C   --   Written by Amy Gilkey - revised 05/31/88
C   --
C   --MSFLAG sets the flags needed to plot a mesh.
C   --
C   --Parameters:
C   --   ANYDEF - OUT - true iff any deformed mesh is to be plotted
C   --   ANYUND - OUT - true iff any undeformed mesh is to be plotted
C   --   DOIXF - OUT - true iff the IXFAC array is needed
C   --   DON2B - OUT - true iff the IN2ELB array is needed
C   --   DOELED - OUT - true iff the deformed element quarilateral centers
C   --      are needed
C   --   DOELEU - OUT - true iff the undeformed element quarilateral centers
C   --      are needed
C   --   DODEAD - OUT - true iff dead nodes are needed
C   --   DONPS - OUT - true iff node set information is needed
C   --   DOESS - OUT - true iff side set information is needed
C   --   DOSCAL - OUT - true iff the zoom window limits need to be calculated
C   --   MINMSH, MAXMSH - OUT - the minimum and maximum mesh line types
C   --      to be displayed
C   --   MAXHID - OUT - the maximum hidden line option
C   --
C   --Common Variables:
C   --   Uses IS3DIM of /D3NUMS/
C   --   Uses DFAC of /DEFORM/
C   --   Uses MSHDEF, MSHNUM, MSHLIN, IHIDOP, NALVAR, DEADNP of /MSHOPT/
C   --   Uses MSCTYP of /MSHLIM/

      PARAMETER (MSHNON=0, MSHBOR=1, MSHDIV=2, MSHSEL=3, MSHALL=4)

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)

      include 'dbnums.blk'
      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM
      COMMON /DEFORM/ DEFPRO, DEFOK, DEFFAC, DDFAC, DFAC,
     &   IXDEF, IYDEF, IZDEF
      LOGICAL DEFPRO, DEFOK

      include 'mshopt.blk'

      COMMON /MSHLIM/ UNMESH(KFAR), ALMESH(KFAR),
     &   ZMMESH(KTOP), RDMESH(KTOP), TICMSH, SQMESH
      LOGICAL SQMESH
      COMMON /MSHLIC/ MSCTYP
      CHARACTER*8 MSCTYP

      LOGICAL ANYDEF, ANYUND
      LOGICAL DOIXF, DON2B, DOELED, DOELEU, DODEAD, DONPS, DOESS, DOSCAL

      INTEGER NUMMOD, NDEFVW, IXVW

C   --Set up calculation flags

      ANYDEF = (NUMMOD (MSHDEF, ' ', 'DEFORM', ' ') .GE. 1)
      ANYUND = (NUMMOD (MSHDEF, ' ', 'UNDEFORM', ' ') .GE. 1)
      IF (DFAC .EQ. 0.0) THEN
         ANYUND = .TRUE.
         ANYDEF = .FALSE.
      END IF

      DOIXF = .FALSE.

      DON2B = (NUMMOD (MSHDEF, MSHNUM, 'DEFORM', 'NODE') .GE. 1)
     &   .OR. (NUMMOD (MSHDEF, MSHNUM, 'DEFORM', 'ALL') .GE. 1)
     &   .OR. (NUMMOD (MSHDEF, MSHNUM, 'UNDEFORM', 'NODE') .GE. 1)
     &   .OR. (NUMMOD (MSHDEF, MSHNUM, 'UNDEFORM', 'ALL') .GE. 1)

      DOELED = (NUMMOD (MSHDEF, MSHNUM, 'DEFORM', 'ELEMENT') .GE. 1)
     &   .OR.  (NUMMOD (MSHDEF, MSHNUM, 'DEFORM', 'ALL') .GE. 1)
      DOELEU = (NUMMOD (MSHDEF, MSHNUM, 'UNDEFORM', 'ELEMENT') .GE. 1)
     &   .OR.  (NUMMOD (MSHDEF, MSHNUM, 'UNDEFORM', 'ALL') .GE. 1)
      IF (.NOT. ANYDEF) THEN
         IF (DOELED) DOELEU = .TRUE.
         DOELED = .FALSE.
      END IF

      DODEAD = (NALVAR .GT. 0) .AND. DEADNP

      DONPS = .FALSE.
      DOESS = .FALSE.
      DO 100 IVW = 1, NDEFVW (.FALSE.)
         IVIEW = IXVW (.FALSE., IVW)
         DONPS = DONPS .OR. (NNPSET(IVIEW) .GT. 0)
         DOESS = DOESS .OR. (NESSET(IVIEW) .GT. 0)
  100 CONTINUE

      DOSCAL = (MSCTYP .EQ. 'EACH')

      IF (IS3DIM) THEN
         MAXHID = IHIDOP
         MINMSH = 999
         MAXMSH = 0
         DO 110 IVW = 1, NDEFVW (.FALSE.)
            IVIEW = IXVW (.FALSE., IVW)
            MINMSH = MIN (MINMSH, MSHLIN(IVIEW))
            MAXMSH = MAX (MAXMSH, MSHLIN(IVIEW))
  110    CONTINUE
      END IF

      RETURN
      END
