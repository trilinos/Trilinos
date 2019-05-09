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

C $Log: maksur.f,v $
C Revision 1.3  2009/03/25 12:36:45  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.2  2009/01/22 21:34:21  gdsjaar
C There were several inline dbnums common blocks. Replaced with the
C include so they all have the same size with the added variable types.
C
C Added minor support for nodeset and sideset variables.
C
C It can print the count and the names, but that is all
C at this time.
C
C Revision 1.1  1994/04/07 20:04:56  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:53:23  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE MAKSUR (LENF, NLNKF, LINKF, DODEAD, IDN2B, NPSURF)
C=======================================================================

C   --*** MAKSUR *** (MESH) Determine which nodes are on surface
C   --   Written by Amy Gilkey - revised 03/10/88
C   --
C   --MAKSUR calculates the nodes-on-surface indicator (NPSURF).
C   --
C   --Parameters:
C   --   LENF - IN - the cumulative face counts by element block
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF - IN - the connectivity for all faces
C   --   DODEAD - IN - add in dead nodes from IDN2B iff true
C   --   IDN2B - IN - the element block for each dead node; dead if >= 0
C   --   NPSURF - OUT - the node numbers of the surface nodes
C   --
C   --Common Variables:
C   --   Uses NDIM, NELBLK of /DBNUMS/
C   --   Sets NNPSUR, NUMNPF of /D3NUMS/

      include 'dbnums.blk'
      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM

      INTEGER LENF(0:NELBLK+2)
      INTEGER NLNKF(NELBLK)
      INTEGER LINKF(*)
      LOGICAL DODEAD
      INTEGER IDN2B(NUMNPF)
      INTEGER NPSURF(*)

      CALL INIINT (NUMNPF, 0, NPSURF)

C   --Mark all nodes that are in surface faces

      DO 120 IELB = 1, NELBLK
         IXL0 = IDBLNK (IELB, 0, LENF, NLNKF) - 1
         DO 110 IFAC = LENF(IELB-1)+1, LENF(IELB)
            DO 100 K = 1, NLNKF(IELB)
               NPSURF(LINKF(IXL0+K)) = 1
  100       CONTINUE
            IXL0 = IXL0 + NLNKF(IELB)
  110    CONTINUE
  120 CONTINUE

C   --Add in dead nodes

      IF (DODEAD) THEN
         DO 130 INP = 1, NUMNPF
            IF (IDN2B(INP) .GE. 0) NPSURF(INP) = 1
  130    CONTINUE
      END IF

C   --Make up the list of surface nodes

      NNPSUR = 0
      DO 140 INP = 1, NUMNPF
         IF (NPSURF(INP) .GT. 0) THEN
            NNPSUR = NNPSUR + 1
            NPSURF(NNPSUR) = INP
         END IF
  140 CONTINUE

      RETURN
      END
