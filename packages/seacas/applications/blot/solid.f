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

C $Log: solid.f,v $
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
C Revision 1.1  1994/04/07 20:13:40  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.3  1992/06/05  22:23:58  gdsjaar
c Fixed problem with ugrcol -- now uses iblk instead of idelb(iblk)
c
c Revision 1.2  1990/12/14  08:58:01  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE SOLID (LENF, NLNKF, LINKF, NXFAC, IXFAC,
     &   XN, YN, ZN, IELBST, BLKCOL, IDELB, *)
C=======================================================================

C   --*** SOLID *** (DETOUR) Paint solid mesh (by index)
C   --   Modified by John Glick - 11/29/88
C   --   Written by Amy Gilkey - revised 10/27/87
C   --
C   --SOLID paints the mesh in the color of each element's element block.
C   --
C   --Parameters:
C   --   LENF - IN - the cumulative face counts by element block
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF - IN - the connectivity for all faces
C   --   NXFAC - IN - the number of ordered faces (if DOIXF)
C   --   IXFAC - IN - the indices of the ordered faces (if DOIXF)
C   --   XN, YN, ZN - IN - the nodal coordinates
C   --   IELBST - IN - the element block status (>0 if selected)
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
C   --   * - return statement if the cancel function is active
C   --
C   --Common Variables:
C   --   Uses NUMEL, NELBLK of /DBNUMS/
C   --   Uses IS3DIM of /D3NUMS/

      include 'dbnums.blk'
      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM

      INTEGER LENF(0:NELBLK)
      INTEGER NLNKF(NELBLK)
      INTEGER LINKF(*)
      INTEGER IXFAC(*)
      REAL XN(*), YN(*), ZN(*)
      INTEGER IELBST(NELBLK)
      INTEGER BLKCOL(0:NELBLK)
      INTEGER IDELB(*)

      LOGICAL GRABRT

      DO 100 IX = 1, NXFAC
         IFAC = IXFAC(IX)
         IELB = 0
         IXL = IDBLNK (IELB, IFAC, LENF, NLNKF)

         IF (GRABRT ()) RETURN 1
         IF (IELB .GT. 0) THEN
c            ITEMP = IDELB(IELB)
            ITEMP = IELB
         ELSE
            ITEMP = IELB
         ENDIF
         CALL UGRCOL (ITEMP, BLKCOL)
         IF ( (.NOT. IS3DIM) .AND. (NLNKF(IELB) .EQ. 9)) THEN
            NNPF = 8
         ELSE
            NNPF = NLNKF(IELB)
         ENDIF
         CALL SOLIDF (NNPF, LINKF(IXL), XN, YN, ZN)
  100 CONTINUE

      CALL PLTFLU

      RETURN
      END
