C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE EPAINT (VARFAC, LENF, NLNKF, LINKF, NXFAC, IXFAC,
     &   XN, YN, ZN, ISVOK, FMIN, FMAX, *)
C=======================================================================

C   --*** EPAINT *** (DETOUR) Paint element contours
C   --   Written by Amy Gilkey - revised 03/07/88
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
C   --   NXFAC - IN - the number of ordered faces
C   --   IXFAC - IN - the indices of the ordered faces
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
      INTEGER IXFAC(*)
      REAL XN(*), YN(*), ZN(*)
      LOGICAL ISVOK(NELBLK)

      LOGICAL GRABRT

C   --Fill the CINTV array with the contour intervals

      IF (.NOT. CINTOK) THEN
         DO 100 NC = 1, NCNTR+1
            CINTV(NC) = CNTRI (NC)
  100    CONTINUE
      END IF

      DO 110 IX = 1, NXFAC
         IFAC = IXFAC(IX)
         IELB = 0
         IXL = IDBLNK (IELB, IFAC, LENF, NLNKF)
         IF (NLNKF(IELB) .GT. 2) THEN

           IF ((.NOT. IS3DIM)  .AND.  (NLNKF(IELB) .EQ. 9)) THEN
             NNPF = 8
           ELSE
             NNPF = NLNKF(IELB)
           ENDIF
           IF (ISVOK(IELB)) THEN

C            --Find contour interval
               IF (DELC .GE. 0.0) THEN
                  NC = LOCREA (VARFAC(IFAC), NCNTR+1, CINTV)
                  IF (VARFAC(IFAC) .LT. CINTV(NC)) NC = NC - 1
               ELSE
                  NC = LOCREA (VARFAC(IFAC), NCNTR+1, CINTV)
                  IF (VARFAC(IFAC) .GE. CINTV(NC)) NC = NC - 1
               END IF

               IF (NOCMIN .AND. (NC .LT. 1)) NC = 1
               IF (NOCMAX .AND. (NC .GT. NCNTR)) NC = NCNTR
               IF ((NC .LT. 1) .OR. (NC .GT. NCNTR)) NC = -1

            ELSE
C            --Not selected, paint face black
               NC = -1
            END IF

C         --Face is entirely inside one contour area

            IF (GRABRT ()) RETURN 1
            CALL GRCOLR (NC)
            CALL SOLIDF (NNPF, LINKF(IXL), XN, YN, ZN)

         END IF
  110 CONTINUE

      CALL PLTFLU

      RETURN
      END
