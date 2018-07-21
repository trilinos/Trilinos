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

C $Log: fndsel.f,v $
C Revision 1.3  2009/03/25 12:36:44  gdsjaar
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
C Revision 1.1  1994/04/07 20:01:14  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:50:49  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE FNDSEL (NODVAR, NNENUM, NENUM,
     &   LENF, NLNKF, LINKF, IF2EL, IX2NP,
     &   NPSURF, NSEL, NPSEL)
C=======================================================================

C   --*** FNDSEL *** (MESH) Find selected nodes/elements
C   --   Written by Amy Gilkey - revised 03/31/88
C   --
C   --FNDSEL constructs a logical array of the selected nodes and elements
C   --that are on the surface of a 3D mesh.
C   --
C   --Parameters:
C   --   NODVAR - IN - true if nodal versus element plot
C   --   NNENUM - IN - the number of selected nodes/elements
C   --   NENUM - IN - the selected nodes/elements
C   --   LENF - IN - the cumulative face counts by element block (if element)
C   --   NLNKF - IN - the number of nodes per face (if element)
C   --   LINKF - IN - the connectivity for all faces (if element)
C   --   IF2EL - IN - the element number of each face (if element)
C   --   IX2NP - IN - the node number for each mesh index (if nodal)
C   --   NPSURF - IN - the node numbers of the surface nodes (if nodal)
C   --   NSEL - OUT - the number of selected nodes
C   --   NPSEL - OUT - the node numbers of the selected nodes
C   --
C   --Common Variables:
C   --   Uses NDIM of /DBNUMS/
C   --   Uses NNPSUR, NUMNPF of /D3NUMS/

      include 'dbnums.blk'
      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM

      LOGICAL NODVAR
      INTEGER NENUM(NNENUM)
      INTEGER LENF(0:NELBLK)
      INTEGER NLNKF(NELBLK)
      INTEGER LINKF(*)
      INTEGER IF2EL(*)
      INTEGER IX2NP(NUMNPF)
      INTEGER NPSURF(*)
      INTEGER NPSEL(*)

      INTEGER IFACES(10)

      IF (NODVAR) THEN
         IF (.NOT. IS3DIM) THEN
            CALL INIINT (NUMNPF, 0, NPSEL)
            DO 100 NNP = 1, NNENUM
               INP = NENUM(NNP)
               IF (INP .GT. 0) NPSEL(INP) = 1
  100       CONTINUE

         ELSE
            CALL INIINT (NUMNPF, -1, NPSEL)
            DO 110 IX = 1, NNPSUR
               INP = NPSURF(IX)
               NPSEL(INP) = 0
  110       CONTINUE
            DO 120 NNP = 1, NNENUM
               INP = NENUM(NNP)
               IF (INP .GT. 0) THEN
                  IF (NPSEL(INP) .EQ. 0) NPSEL(INP) = 1
               END IF
  120       CONTINUE
         END IF

      ELSE
         CALL INIINT (NUMNPF, 0, NPSEL)

         DO 150 NEL = 1, NNENUM
            IEL = NENUM(NEL)
            CALL FNDE2F (IEL, LENF, IF2EL, NFARY, IFACES, IELB)
            DO 140 N = 1, NFARY
               IFAC = IFACES(N)
               IELB = 0
               IXL0 = IDBLNK (IELB, IFAC, LENF, NLNKF) - 1
               DO 130 K = 1, NLNKF(IELB)
                  INP = LINKF(IXL0+K)
                  NPSEL(INP) = 1
  130          CONTINUE
  140       CONTINUE
  150    CONTINUE
      END IF

C   --Make up the list of selected nodes

      NSEL = 0
      DO 160 INP = 1, NUMNPF
         IF (NPSEL(INP) .GT. 0) THEN
            NSEL = NSEL + 1
            NPSEL(NSEL) = INP
         END IF
  160 CONTINUE

      RETURN
      END
