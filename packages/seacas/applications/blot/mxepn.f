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

C $Log: mxepn.f,v $
C Revision 1.4  2009/03/25 12:36:46  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.3  2009/01/22 21:34:21  gdsjaar
C There were several inline dbnums common blocks. Replaced with the
C include so they all have the same size with the added variable types.
C
C Added minor support for nodeset and sideset variables.
C
C It can print the count and the names, but that is all
C at this time.
C
C Revision 1.2  1999/03/09 22:02:50  gdsjaar
C Fixed problem with incorrect 'mesh contiguity' errors. Max number of
C faces connected to a node was being based on block-by-block counting
C of elements connected to a node which didn't work well for meshes
C containing lots of blocks. Changed to do count based on all element blocks.
C
C Revision 1.1  1994/04/07 20:06:04  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:54:12  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE  MXEPN (NUME, NLNKE, LINKE, NFPN)
C=======================================================================

C   --*** MXEPN *** (MESH) Count number of elements per node
C   --   Written by Amy Gilkey - revised 10/26/87
C   --
C   --MXEPN counts the maximum number of elements that a node is in.
C   --
C   --Parameters:
C   --   NUME - IN - the number of elements
C   --   NLNKE - IN - the number of nodes per element
C   --   LINKE - IN - the original connectivity (some nodes may be 0)
C   --   NFPN - SCRATCH - size = 1+NUMNP
C   --
C   --Common Variables:
C   --   Uses NUMNP of /DBNUMS/

      common /debugc/ cdebug
      common /debugn/ idebug
      character*8 cdebug

      include 'dbnums.blk'
      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM

      INTEGER LINKE(NLNKE,NUME)
      INTEGER NFPN(0:NUMNP)

      DO 110 IEL = 1, NUME
         DO 100 ILINK = 1, NLNKE
            INP = LINKE(ILINK,IEL)
            NFPN(INP) = NFPN(INP) + 1
  100    CONTINUE
  110 CONTINUE
      RETURN
      END

      integer function mxepnmx(nfpn)
      include 'dbnums.blk'

      integer nfpn(0:numnp)
      MXEPNMX = 0
      DO 120 INP = 1, NUMNP
         MXEPNMX = MAX (MXEPNMX, NFPN(INP))
  120 CONTINUE
      RETURN
      END
