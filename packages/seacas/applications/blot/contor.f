C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

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
      include 'd3nums.blk'
      include 'cntr.blk'

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
