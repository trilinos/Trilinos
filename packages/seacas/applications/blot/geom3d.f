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

C $Log: geom3d.f,v $
C Revision 1.4  2009/03/25 12:36:44  gdsjaar
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
C Revision 1.2  2004/10/18 16:30:00  gdsjaar
C Add capability to handle tet elements.
C
C There are still a few minor problems, but the display is substantially correct.
C
C Revision 1.1  1994/04/07 20:01:28  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:51:00  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE GEOM3D (LENF, NLNKF, LINKF, IF2EL,
     &   LENL, LINSET, LINDEF, LREF)
C=======================================================================

C   --*** GEOM3D *** (MESH) Group 3D lines connecting nodes
C   --   Written by Amy Gilkey - revised 03/29/88
C   --              Sam Key, 04/85
C   --
C   --GEOM3D finds all lines making up the surface faces.
C   --The lines are classified and sorted into the following parts:
C   --   -1) Edge
C   --    0) Element block boundary
C   --    n) Interior within element block 'n'
C   --
C   --Parameters:
C   --   LENF - IN - the cumulative face counts by element block
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF - IN - the connectivity for all faces
C   --   IF2EL - IN - the element number of each face
C   --   LENL - OUT - the cumulative line counts by element block
C   --   LINSET - OUT - the sorted line set, for surface faces only
C   --   LINDEF - SCRATCH - the line set scratch array;
C   --      size = 6 * maximum number of lines
C   --   LREF - SCRATCH - the last line set index for node, length = NUMNPF
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/
C   --   Uses NUMNPF, LLNSET of /D3NUMS/

      common /debugc/ cdebug
      common /debugn/ idebug
      character*8 cdebug

      include 'dbnums.blk'
      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM

      INTEGER LENF(0:NELBLK)
      INTEGER NLNKF(NELBLK)
      INTEGER LINKF(*)
      INTEGER IF2EL(*)
      INTEGER LENL(-2:NELBLK), LINSET(LLNSET,*)
      INTEGER LINDEF(0:5,*)
      INTEGER LREF(*)


      CALL INIINT (NUMNPF, 0, LREF)

C   --Define line set for surface faces only

      NHEAP = 0
      DO 170 IELB = 1, NELBLK
         IXL0 = IDBLNK (IELB, 0, LENF, NLNKF) - 1

C      --Handle normal multi-sided elements

         IF (NLNKF(IELB) .GT. 2) THEN
            DO 130 IFAC = LENF(IELB-1)+1, LENF(IELB)
               IEL = IF2EL(IFAC)

C            --Extract element lines

               N2 = LINKF(IXL0+NLNKF(IELB))
               LREFN2 = LREF(N2)
               LSTREF = LREFN2

               DO 120 ILINK = 1, NLNKF(IELB)
                  N1 = N2
                  N2 = LINKF(IXL0+ILINK)
                  if (n1 .eq. n2) goto 120
                  NMIN = MIN (N1, N2)
                  NMAX = MAX (N1, N2)

C               --Search for line in existing lines

                  LREFN1 = LREFN2
                  IF (ILINK .LT. NLNKF(IELB)) THEN
                     LREFN2 = LREF(N2)
                  ELSE
                     LREFN2 = LSTREF
                  END IF

                  IL = MIN (LREFN1, LREFN2)
  100             CONTINUE
                  IF (IL .GT. 0) THEN
                     IF (LINDEF(1,IL) .EQ. NMIN) THEN
                        IF (LINDEF(2,IL) .EQ. NMAX) GOTO 110
                        IL = LINDEF(4,IL)
                     ELSE IF (LINDEF(1,IL) .EQ. NMAX) THEN
                        IL = LINDEF(4,IL)
                     ELSE
                        IL = LINDEF(5,IL)
                     END IF
                     GOTO 100
                  END IF
  110             CONTINUE

                  IF (IL .GT. 0) THEN
                     IOLELB = IABS (LINDEF(0,IL))

                     IF (IELB .EQ. IOLELB) THEN
C                     --If both faces of the same element (i.e., an edge),
C                     --change from element block n set to edge set,
C                     --else change to element block n set
                        IF (IEL .EQ. LINDEF(3,IL)) THEN
                           LINDEF(0,IL) = -1
                        ELSE
                           LINDEF(0,IL) = IELB
                        END IF

                     ELSE
C                     --If not same element block, change line from
C                     --element block n set to element block boundary set
                        LINDEF(0,IL) = 0

                     END IF
                  ELSE

C                  --No match, add to edge set and point to face

                     NHEAP = NHEAP + 1
                     LINDEF(0,NHEAP) = -IELB
                     LINDEF(1,NHEAP) = NMIN
                     LINDEF(2,NHEAP) = NMAX
                     LINDEF(3,NHEAP) = IEL
                     LINDEF(4,NHEAP) = LREF(NMIN)
                     LINDEF(5,NHEAP) = LREF(NMAX)

                     LREF(N1) = NHEAP
                     LREF(N2) = NHEAP
                  END IF
  120          CONTINUE

               IXL0 = IXL0 + NLNKF(IELB)
  130       CONTINUE

C      --Handle 2-node truss elements

         ELSE IF (NLNKF(IELB) .EQ. 2) THEN
            DO 150 IFAC = LENF(IELB-1)+1, LENF(IELB)
               IEL = IF2EL(IFAC)

C            --Extract element lines for 2-node elements

               N2 = LINKF(IXL0+NLNKF(IELB))
               DO 140 ILINK = 1, NLNKF(IELB)
                  N1 = N2
                  N2 = LINKF(IXL0+ILINK)
                  NMIN = MIN (N1, N2)
                  NMAX = MAX (N1, N2)

C               --Add to edge set as element block n set

                  NHEAP = NHEAP + 1
                  LINDEF(0,NHEAP) = IELB
                  LINDEF(1,NHEAP) = NMIN
                  LINDEF(2,NHEAP) = NMAX
                  LINDEF(3,NHEAP) = IEL
                  LINDEF(4,NHEAP) = 0
                  LINDEF(5,NHEAP) = 0
  140          CONTINUE

               IXL0 = IXL0 + NLNKF(IELB)
  150       CONTINUE

         ELSE
            DO 160 IFAC = LENF(IELB-1)+1, LENF(IELB)
               IXL0 = IXL0 + NLNKF(IELB)
  160       CONTINUE
         END IF
  170 CONTINUE

C   --Count lines in each line set

      DO 180 IELB = -2, NELBLK
         LENL(IELB) = 0
  180 CONTINUE

      DO 190 IL = 1, NHEAP
         IF (LINDEF(0,IL) .LT. -1) LINDEF(0,IL) = -1
         IELB = LINDEF(0,IL)
         LENL(IELB) = LENL(IELB) + 1
  190 CONTINUE

C   --Get starting indices for each line set

      LENL(-2) = 0

      NHEAP = 0
      DO 200 IELB = -1, NELBLK
         L = LENL(IELB)
         LENL(IELB) = NHEAP
         NHEAP = NHEAP + L
  200 CONTINUE

C   --Move lines into appropriate line set

      DO 210 IL = 1, NHEAP
         IELB = LINDEF(0,IL)
         IX = LENL(IELB) + 1
         LENL(IELB) = IX
C...         CALL CPYINT (2, LINDEF(1,IL), LINSET(1,IX))
         linset(1,ix) = lindef(1,il)
         linset(2,ix) = lindef(2,il)
  210 CONTINUE


      RETURN
      END
