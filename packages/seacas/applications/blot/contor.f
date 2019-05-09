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

C $Log: contor.f,v $
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
C Revision 1.2  1997/11/11 14:55:54  gdsjaar
C Added 'external blkdat' to main program to ensure that the block data
C gets linked into the executable. Wasn't happening on dec alpha
C systems.
C
C Removed unreachable lines in several routines
C
C Fixed variable name spelling in contor.f
C
C Unsplit strings that were split across lines
C
C Removed old error variables left over from exodusIIv1
C
C Upped version number
C
C Revision 1.1  1994/04/07 19:57:18  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:49:04  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE CONTOR (VARNP, LENF, NLNKF, LINKF, HIDEF,
     &   XN, YN, ZN, XF, YF, ZF, LENL, LINSET,
     &   IN2ELB, ISVOK, FMIN, FMAX, *)
C=======================================================================

C   --*** CONTOR *** (DETOUR) Plot line contours
C   --   Written by Amy Gilkey - revised 10/29/87
C   --   D. P. Flanagan, 03/30/83
C   --
C   --CONTOR draws contour lines with alphabetic tags in a color
C   --sequence.  The contour algorithm assumes that the elements do not
C   --contain internal nodes.  The element interpolation field is
C   --approximated by logically drawing lines from each node to the element
C   --center and, thusly, dividing the element into triangles.  Contour
C   --lines are then drawn by connecting the intersection points of the
C   --sub-element edges and the contour plane.  The contour letter is
C   --displayed at points where the contour plane intersects a mesh boundary
C   --or element block boundary line.  Also, every "LABINC" interior mesh line
C   --which intersects the contour plane is labelled with the contour letter.
C   --
C   --The element block status indicates which element blocks are active
C   --for contours.  No contours are drawn in inactive elements.  Mesh lines
C   --are inactive if either of its nodes is not connected to an active
C   --element.  Inactive mesh lines are not labeled with contour letters.
C   --
C   --Parameters:
C   --   VARNP - IN - the contour function values
C   --   LENF - IN - the cumulative face counts by element block
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF - IN - the connectivity for all faces
C   --   HIDEF(i) - IN - true iff face i is hidden (3D only)
C   --   XN, YN, ZN - IN - the nodal coordinates
C   --   XF, YF, ZF - IN - the face center coordinates
C   --   LENL - IN - the cumulative line counts by element block
C   --   LINSET - IN - the sorted line set
C   --   IN2ELB - IN - the element block for each node;
C   --      <0 if not in any selected element block
C   --      =0 if in more than one selected element block
C   --   ISVOK - IN - ISVOK(i) is true iff the contour variable is defined
C   --      for element block i
C   --   FMIN, FMAX - IN - the minimum and maximum function values
C   --   * - return statement if the cancel function is active
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/
C   --   Uses IS3DIM, LLNSET of /D3NUMS/
C   --   Uses NCNTR, LABINC, MAXMIN, MAXMAX of /CNTR/

      include 'dbnums.blk'
      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM
      COMMON /CNTR/   CINTOK, LINCON, NCNTR, CMIN, CMAX, DELC,
     &   CINTV(256), NOCMIN, NOCMAX, LABINC, MAXMIN, MAXMAX
      LOGICAL CINTOK, LINCON, NOCMIN, NOCMAX

      REAL VARNP(*)
      INTEGER LENF(0:NELBLK)
      INTEGER NLNKF(NELBLK)
      INTEGER LINKF(*)
      LOGICAL HIDEF(*)
      REAL XN(*), YN(*), ZN(*)
      REAL XF(*), YF(*), ZF(*)
      INTEGER LENL(-2:NELBLK), LINSET(LLNSET,*)
      INTEGER IN2ELB(*)
      LOGICAL ISVOK(NELBLK)

      LOGICAL GRABRT

      DO 150 NC = 1, NCNTR
         CNTR = CNTRI (NC)

C      --Skip this contour if the values are outside the contour value
         IF ((FMIN .LE. CNTR) .AND. (FMAX .GE. CNTR)) THEN

            CALL GRCOLR (NC)

            DO 130 IELB = 1, NELBLK
               IF (ISVOK(IELB) .AND. (NLNKF(IELB) .GT. 2)) THEN

                  DO 100 IFAC = LENF(IELB-1)+1, LENF(IELB)
                     IF (IS3DIM) THEN
                        IF (HIDEF(IFAC)) GOTO 100
                     END IF

                     IF (GRABRT ()) RETURN 1
                     IXL = IDBLNK (IELB, IFAC, LENF, NLNKF)
                     CALL CONTF (CNTR, VARNP, NLNKF(IELB), LINKF(IXL),
     &                  XN, YN, ZN, XF(IFAC), YF(IFAC), ZF(IFAC))
  100             CONTINUE

C               --Label interior lines

                  IF (LABINC .GT. 0) THEN
                     NHIT = 0
                     DO 110 IL = LENL(IELB-1)+1, LENL(IELB)
                        CALL CONLAB (NC, CNTR, NHIT, LABINC,
     &                     LINSET(1,IL), VARNP, XN, YN, ZN, IN2ELB,
     &                     *160)
  110                CONTINUE
                  END IF

C               --Label interior lines which are on the edge

                  IF (IS3DIM .AND. (LABINC .GE. 0)) THEN
                     NHIT = 0
                     DO 120 IL = LENL(IELB-1)+1, LENL(IELB)
C                     --Label edge line only
                        IF (LINSET(3,IL) .LT. 0) THEN
                           CALL CONLAB (NC, CNTR, NHIT, 1, LINSET(1,IL),
     &                        VARNP, XN, YN, ZN, IN2ELB, *160)
                        END IF
  120                CONTINUE
                  END IF

               END IF
  130       CONTINUE

C         --Label mesh boundary and element block boundary lines

            IF (LABINC .GE. 0) THEN
               NHIT = 0
               DO 140 IL = 1, LENL(0)
                  CALL CONLAB (NC, CNTR, NHIT, 1, LINSET(1,IL),
     &               VARNP, XN, YN, ZN, IN2ELB, *160)
  140          CONTINUE
            END IF

            CALL PLTFLU
         END IF
  150 CONTINUE

      RETURN

  160 CONTINUE
      RETURN 1
      END
