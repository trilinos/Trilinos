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

C $Log: surf2d.f,v $
C Revision 1.3  2009/03/25 12:36:48  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.2  2009/01/22 21:34:22  gdsjaar
C There were several inline dbnums common blocks. Replaced with the
C include so they all have the same size with the added variable types.
C
C Added minor support for nodeset and sideset variables.
C
C It can print the count and the names, but that is all
C at this time.
C
C Revision 1.1  1994/04/07 20:15:58  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:58:50  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE SURF2D (LENE, NLNKE, LINKE,
     &   LENF, NLNKF, LINKF, IF2EL, NFACES, LENLNK)
C=======================================================================

C   --*** SURF2D *** (MESH) Group 2D faces
C   --   Written by Amy Gilkey - revised 11/04/87
C   --
C   --SURF2D packs the elements by element block, eliminating empty elements.
C   --
C   --Parameters:
C   --   LENE - IN - the cumulative element counts by element block
C   --   NLNKE - IN - the number of nodes per element
C   --   LINKE - IN - the original connectivity; connectivity all zero
C   --      if element is undefined
C   --   LENF - OUT - the cumulative face counts by element block
C   --   NLNKF - OUT - the number of nodes per face
C   --   LINKF - OUT - the connectivity for all faces
C   --   IF2EL - OUT - the element number of each face
C   --   NFACES - OUT - the number of returned faces
C   --   LENLNK - OUT - the length of the returned LINKF array
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/

      common /debugc/ cdebug
      common /debugn/ idebug
      character*8 cdebug

      include 'dbnums.blk'
      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM

      INTEGER LENE(0:NELBLK), LINKE(*)
      INTEGER NLNKE(NELBLK)
      INTEGER LENF(0:NELBLK+2)
      INTEGER NLNKF(NELBLK)
      INTEGER LINKF(*)
      INTEGER IF2EL(*)

      LENF(0) = 0
      NFACES = 0
      IXL0 = 0
      DO 110 IELB = 1, NELBLK
         NLNKF(IELB) = NLNKE(IELB)
         IXE0 = IDBLNK (IELB, 0, LENE, NLNKE) - 1
         DO 100 IEL = LENE(IELB-1)+1, LENE(IELB)
            IF (LINKE(IXE0+1) .NE. 0) THEN
               NFACES = NFACES + 1
               CALL CPYINT (NLNKF(IELB), LINKE(IXE0+1), LINKF(IXL0+1))
               IF2EL(NFACES) = IEL
               IXL0 = IXL0 + NLNKF(IELB)
            END IF
            IXE0 = IXE0 + NLNKE(IELB)
  100    CONTINUE
         LENF(IELB) = NFACES
  110 CONTINUE

      LENF(NELBLK+1) = NFACES
      LENF(NELBLK+2) = NFACES

      LENLNK = IXL0
      RETURN
      END
