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

C $Log: deadup.F,v $
C Revision 1.4  2009/03/25 12:36:43  gdsjaar
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
C Revision 1.2  1997/06/26 17:51:32  caforsy
C Substituted if NeedsDoubleEscape for ifdef sun.  This define is needed on cougar, sun and ibm platforms
C
C Revision 1.1  1994/04/07 19:59:17  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
CRevision 1.3  1993/07/08  21:57:11  gdsjaar
CReplaced SUN with sun in all ifdef statements.
C
CRevision 1.2  1990/12/14  08:49:21  gdsjaar
CAdded RCS Id and Log to all files
C
C=======================================================================
      SUBROUTINE DEADUP (HIDENP, XN, YN, ZN, IDN2B, *)
C=======================================================================

C   --*** DEADUP *** (MESH) Plot "dead" nodes with symbols
C   --   Written by Amy Gilkey - revised 03/09/88
C   --
C   --DEADUP draws a symbol at each dead node.  A dead node is a node
C   --which is not attached to any alive elements.
C   --Only nodes of selected element blocks are drawn.
C   --
C   --Parameters:
C   --   HIDENP(i) - IN - true iff node i is hidden (3D only)
C   --   XN, YN, ZN - IN - the nodal coordinates
C   --   IDN2B - IN - the element block for each dead node; dead if >= 0;
C   --      =0 if in more than one selected element block
C   --   * - return statement if the cancel function is active
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/
C   --   Uses IS3DIM, NUMNPF of /D3NUMS/

      include 'dbnums.blk'
      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM

      LOGICAL HIDENP(NUMNPF)
      REAL XN(NUMNPF), YN(NUMNPF), ZN(NUMNPF)
      INTEGER IDN2B(NUMNPF)

      LOGICAL GRABRT

      DO 110 IELB = 0, NELBLK

C      --Set the node element block color
         CALL GRCOLR (IELB)

         DO 100 INP = 1, NUMNPF

            IF (IS3DIM) THEN
               IF (HIDENP(INP)) GOTO 100
            END IF

C         --If node is dead and of the appropriate element block, draw symbol

            IF (IDN2B(INP) .EQ. IELB) THEN
               IF (GRABRT ()) RETURN 1
#if NeedsDoubleEscape
               CALL MPD2SY (1, XN(INP), YN(INP), '\\CDO')
#else
               CALL MPD2SY (1, XN(INP), YN(INP), '\CDO')
#endif
            END IF

  100    CONTINUE

         CALL PLTFLU
  110 CONTINUE

      RETURN
      END
