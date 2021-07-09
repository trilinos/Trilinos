C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CONTF (CNTR, VARNP, NLNKF, LINKF1,
     &   XN, YN, ZN, XF, YF, ZF)
C=======================================================================

C   --*** CONTF *** (DETOUR) Plot line contours for a face
C   --   Written by Amy Gilkey - revised 10/22/87
C   --   D. P. Flanagan, 03/30/83
C   --
C   --CONTF draws a contour line for a face.  The contour algorithm
C   --assumes that the elements do not contain internal nodes.
C   --The element interpolation field is approximated by logically drawing
C   --lines from each node to the element center and, thusly, dividing the
C   --element into triangles.  Contour lines are then drawn by connecting
C   --the intersection points of the sub-element edges and the contour plane.
C   --
C   --Parameters:
C   --   CNTR - the contour value
C   --   VARNP - IN - the contour function values
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF1 - IN - the connectivity for the face
C   --   XN, YN, ZN - IN - the nodal coordinates
C   --   XF, YF, ZF - IN - the face center coordinates
C   --   * - return statement if the cancel function is active
C   --
C   --Common Variables:
C   --   Uses IS3DIM of /D3NUMS/

      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM

      REAL VARNP(*)
      INTEGER LINKF1(NLNKF)
      REAL XN(*), YN(*), ZN(*)

      LOGICAL INTERP_BL
      LOGICAL CROSSL, CROSSN

C   --Compute the minimum and maximum values for the face

      FEMAX = VARNP(LINKF1(1))
      FEMIN = VARNP(LINKF1(1))
      DO 100 K = 2, NLNKF
         FEMAX = MAX (FEMAX, VARNP(LINKF1(K)))
         FEMIN = MIN (FEMIN, VARNP(LINKF1(K)))
  100 CONTINUE

      IF ((FEMIN .LE. CNTR) .AND. (FEMAX .GE. CNTR)) THEN

C      --Contour line runs through the face

         IF (NLNKF .EQ. 3) THEN

C         --Special case - triangle element

            N = LINKF1(2)
            FM = VARNP(N)
            XM = XN(N)
            YM = YN(N)
            NTRI = 1
            NNPF = NLNKF

         ELSE IF ((.NOT. IS3DIM) .AND. (NLNKF .EQ. 8)) THEN

C         --Compute center function value and coordinates

            FM = 0.0
            DO 110 K = 1, NLNKF, 2
               FM = FM - .25 * VARNP(LINKF1(K))
  110       CONTINUE
            DO 120 K = 2, NLNKF, 2
               FM = FM + .5 * VARNP(LINKF1(K))
  120       CONTINUE
            FM = FM
            XM = XF
            YM = YF
            NTRI = NLNKF
            NNPF = NLNKF

         ELSE IF ((.NOT. IS3DIM) .AND. (NLNKF .EQ. 9)) THEN

C         --Get center function value and coordinates

            FM = VARNP(LINKF1(9))
            XM = XN(LINKF1(9))
            YM = YN(LINKF1(9))
            NTRI = NLNKF - 1
            NNPF = NLNKF - 1

         ELSE

C         --Compute center function value and coordinates

            FM = 0.0
            DO 130 K = 1, NLNKF
               FM = FM + VARNP(LINKF1(K))
  130       CONTINUE
            FM = FM / NLNKF
            XM = XF
            YM = YF
            NTRI = NLNKF
            NNPF = NLNKF
         END IF

         N = LINKF1(NNPF)
         CROSSN = (INTERP_BL (CNTR, FM, VARNP(N), PSI))
         IF (CROSSN) THEN
            XCN = XM * (1.-PSI) + XN(N) * PSI
            YCN = YM * (1.-PSI) + YN(N) * PSI
         END IF

         DO 140 K = 1, NTRI
            L = N
            N = LINKF1(K)
            CROSSL = CROSSN
            IF (CROSSL) THEN
               XCL = XCN
               YCL = YCN
            END IF
            CROSSN = (INTERP_BL (CNTR, FM, VARNP(N), PSI))
            IF (CROSSN) THEN
               XCN = XM * (1.-PSI) + XN(N) * PSI
               YCN = YM * (1.-PSI) + YN(N) * PSI
            END IF

            IF (CROSSL .OR. CROSSN) THEN
               IF (INTERP_BL (CNTR, VARNP(N), VARNP(L), PSI)) THEN
                  XC0 = XN(N) * (1.-PSI) + XN(L) * PSI
                  YC0 = YN(N) * (1.-PSI) + YN(L) * PSI
                  IF (CROSSL .AND. CROSSN) THEN
C                  --There are three cases where the contour crosses the
C                  --side and both middle lines: contour = middle point
C                  --(go either way), contour = side last corner (go to
C                  --next), contour = side next corner (go to last)
                     IF ((ABS (XC0 - XCL) + ABS (YC0 - YCL)) .LT.
     &                  (ABS (XC0 - XCN) + ABS (YC0 - YCN)))
     &                  CROSSL = .FALSE.
                  END IF
                  IF (CROSSL) THEN
                     CALL MPD2VC (1, XCL, YCL, XC0, YC0)
                  ELSE
                     CALL MPD2VC (1, XC0, YC0, XCN, YCN)
                  END IF
               ELSE IF (CROSSL .AND. CROSSN) THEN
                  CALL MPD2VC (1, XCL, YCL, XCN, YCN)
               ELSE
C               --If the contour crosses the next middle line but not
C               --the side or the last middle line (impossible with
C               --absolute accuracy), ignore crossing.
                  CROSSN = .FALSE.
               END IF
            END IF
  140    CONTINUE
      END IF

      RETURN
      END
