C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE QEPAIN (VARFAC, LENF, NLNKF, LINKF, HIDEF,
     &   XN, YN, ZN, ISVOK, FMIN, FMAX, *)
C=======================================================================

C   --*** QEPAIN *** (DETOUR) Paint element contours (quick)
C   --   Written by Amy Gilkey - revised 03/02/88
C   --
C   --EPAINT paints contour sections in a color sequence.  Each element
C   --is assigned a single contour value (based on the element variable
C   --value.
C   --
C   --The element block status indicates which element blocks are active
C   --for contours.  No contours are drawn in inactive elements.
C   --
C   --Parameters:
C   --   VARFAC - IN - the contour function values
C   --   LENF - IN - the cumulative face counts by element block
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF - IN - the connectivity for all faces
C   --   HIDEF(i) - IN - true iff face i is hidden (3D only)
C   --   XN, YN, ZN - IN - the nodal coordinates
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

      REAL VARFAC(*)
      INTEGER LENF(0:NELBLK)
      INTEGER NLNKF(NELBLK)
      INTEGER LINKF(*)
      LOGICAL HIDEF(*)
      REAL XN(*), YN(*), ZN(*)
      LOGICAL ISVOK(NELBLK)

      LOGICAL GRABRT

      DO 120 NC = 1, NCNTR
         IF (DELC .GE. 0.0) THEN
            CNTR0 = CNTRI (NC)
            CNTR1 = CNTRI (NC+1)
            IF (NOCMIN .AND. (NC .EQ. 1)) CNTR0 = FMIN - DELC
            IF (NOCMAX .AND. (NC .EQ. NCNTR)) CNTR1 = FMAX + DELC
         ELSE
            CNTR1 = CNTRI (NC)
            CNTR0 = CNTRI (NC+1)
            IF (NOCMIN .AND. (NC .EQ. 1)) CNTR1 = FMAX - DELC
            IF (NOCMAX .AND. (NC .EQ. NCNTR)) CNTR0 = FMIN + DELC
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

C                  --If face is inside the contour area, paint it

                     IF ((VARFAC(IFAC) .GE. CNTR0)
     &                  .AND. (VARFAC(IFAC) .LT. CNTR1)) THEN
                        IF (GRABRT ()) RETURN 1
                        IXL = IDBLNK (IELB, IFAC, LENF, NLNKF)
                        CALL SOLIDF (NNPF, LINKF(IXL), XN, YN, ZN)
                     END IF
  100             CONTINUE

               END IF
  110       CONTINUE

            CALL PLTFLU
         END IF
  120 CONTINUE

      RETURN
      END
