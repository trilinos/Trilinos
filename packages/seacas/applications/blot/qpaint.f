C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE QPAINT (VARNP, LENF, NLNKF, LINKF, HIDEF,
     &   XN, YN, ZN, XF, YF, ZF, ISVOK, FMIN, FMAX, *)
C=======================================================================

C   --*** QPAINT *** (DETOUR) Paint contours (quick)
C   --   Modified by John H. Glick - 10/26/88
C   --   Written by Amy Gilkey - revised 03/14/88
C   --
C   --QPAINT paints contour sections in a color sequence.  The contour
C   --algorithm assumes that the elements do not contain internal nodes.
C   --The element interpolation field is approximated by logically drawing
C   --lines from each node to the element center and, thusly, dividing the
C   --element into triangles.  Contour sections are then drawn by connecting
C   --the intersection points of the sub-element edges and the contour
C   --plane.
C   --
C   --The element block status indicates which element blocks are active
C   --for contours.  No contours are drawn in inactive elements.
C   --
C   --Parameters:
C   --   VARNP - IN - the contour function values
C   --   LENF - IN - the cumulative face counts by element block
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF - IN - the connectivity for all faces
C   --   HIDEF(i) - IN - true iff face i is hidden (3D only)
C   --   XN, YN, ZN - IN - the nodal coordinates
C   --   XF, YF, ZF - IN - the face center coordinates
C   --   ISVOK - IN - ISVOK(i) is true iff the contour variable is defined
C   --      for element block i (always true if nodal variable)
C   --   FMIN, FMAX - IN - the minimum and maximum function value
C   --   * - return statement if the cancel function is active
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/
C   --   Uses IS3DIM of /D3NUMS/
C   --   Uses NCNTR, NOCMIN, NOCMAX of /CNTR/

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
      LOGICAL ISVOK(NELBLK)

C   --Set the contour minimum and maximum (in case NOCMIN, etc.)
      CNTRMN = FMIN - ABS (DELC)
      CNTRMX = FMAX + ABS (DELC)

      DO 120 NC = 1, NCNTR
         IF (DELC .GE. 0.0) THEN
            CNTR0 = CNTRI (NC)
            IF (NOCMIN .AND. (NC .EQ. 1)) CNTR0 = CNTRMN
            CNTR1 = CNTRI (NC+1)
            IF (NOCMAX .AND. (NC .EQ. NCNTR)) CNTR1 = CNTRMX
         ELSE
            CNTR1 = CNTRI (NC)
            IF (NOCMIN .AND. (NC .EQ. 1)) CNTR1 = CNTRMX
            CNTR0 = CNTRI (NC+1)
            IF (NOCMAX .AND. (NC .EQ. NCNTR)) CNTR0 = CNTRMN
         END IF

C      --Skip this contour if the values are outside the contour range

         IF ((FMIN .LT. CNTR1) .AND. (FMAX .GE. CNTR0)) THEN

            CALL GRCOLR (NC)

            DO 110 IELB = 1, NELBLK

               IF ((.NOT. IS3DIM)  .AND.  (NLNKF(IELB) .EQ. 9)) THEN
                  NNPF = 8
               ELSE
                  NNPF = NLNKF(IELB)
               ENDIF

               IF (ISVOK(IELB) .AND. (NLNKF(IELB) .GT. 2)) THEN

                  DO 100 IFAC = LENF(IELB-1)+1, LENF(IELB)
                     IF (IS3DIM) THEN
                        IF (HIDEF(IFAC)) GOTO 100
                     END IF

                     IXL = IDBLNK (IELB, IFAC, LENF, NLNKF)

C                  --Compute the minimum and maximum values for the face
                     if (nlnkf(ielb) .eq. 4) then
                       femax = max(varnp(linkf(ixl)),
     *                   varnp(linkf(ixl+1)), varnp(linkf(ixl+2)),
     *                   varnp(linkf(ixl+3)))
                       femin = min(varnp(linkf(ixl)),
     *                   varnp(linkf(ixl+1)), varnp(linkf(ixl+2)),
     *                   varnp(linkf(ixl+3)))
                     else
                       CALL CFMAX (VARNP, NLNKF(IELB), LINKF(IXL),
     &                   FEMIN, FEMAX)
                     end if

                     IF ((FEMIN .GE. CNTR0) .AND. (FEMAX .LT. CNTR1))
     &                  THEN

C                     --Face is entirely inside the contour area
                        CALL SOLIDF (NNPF, LINKF(IXL),
     &                     XN, YN, ZN)

                     ELSE IF ((FEMIN .LT. CNTR1)
     &                  .AND. (FEMAX .GE. CNTR0)) THEN

C                     --Face is partially inside the contour area
                        CALL PAINTF (CNTR0, CNTR1,
     &                     VARNP, NLNKF(IELB), LINKF(IXL),
     &                     XN, YN, ZN, XF(IFAC), YF(IFAC), ZF(IFAC))
                     END IF
  100             CONTINUE

               END IF
  110       CONTINUE

            CALL PLTFLU
         END IF
  120 CONTINUE

      RETURN
      END
