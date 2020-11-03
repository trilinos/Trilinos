C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SCAELB (A, USESEL, IELBST,
     &   VALMN, NUMMN, XYZMN, ISTMN, VALMX, NUMMX, XYZMX, ISTMX,
     &   VALMIN, NUMMIN, XYZMIN, ISTMIN, VALMAX, NUMMAX, XYZMAX, ISTMAX)
C=======================================================================

C   --*** SCAELB *** (BLOT) Scale element variable for selected blocks
C   --   Written by Amy Gilkey - revised 04/01/88
C   --
C   --SCAELB returns the minimum and maximum values for the selected
C   --element blocks.  The minimum and maximums for each block have already
C   --been calculated and stored.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   USESEL - IN - use the element blocks selected array iff true,
C   --      else all selected
C   --   IELBST - IN - the element block status (>0 if selected)
C   --   VALMN, VALMX - IN - the minimum and maximum value for each
C   --      element block
C   --   NUMMN, NUMMX - IN - the element number of the minimum and
C   --      maximum value for each element block
C   --   XYZMN, XYZMX - IN - the coordinates of NUMMN, NUMMX
C   --   ISTMN, ISTMX - IN - the step number of the minimum and maximum
C   --      value for each element block
C   --   VALMIN, VALMAX - OUT - the minimum and maximum value
C   --      (with selected element block and birth/death)
C   --   NUMMIN, NUMMAX - OUT - the element number of the minimum and
C   --      maximum value
C   --   XYZMIN, XYZMAX - OUT - the coordinates of NUMMIN, NUMMAX
C   --   ISTMIN, ISTMAX - OUT - the step number of the minimum and maximum value
C   --
C   --Common Variables:
C   --   Uses NUMEL, NELBLK, NVAREL, NSTEPS of /DBNUMS/

      include 'dbnums.blk'

      DIMENSION A(*)
      LOGICAL USESEL
      INTEGER IELBST(*)
      REAL VALMN(0:NELBLK), VALMX(0:NELBLK)
      INTEGER NUMMN(0:NELBLK), NUMMX(0:NELBLK)
      REAL XYZMN(3,0:NELBLK), XYZMX(3,0:NELBLK)
      INTEGER ISTMN(0:NELBLK), ISTMX(0:NELBLK)
      REAL XYZMIN(3), XYZMAX(3)

      LOGICAL INIT

C   --Calculate min/max for element by element block

      IF (USESEL) THEN
         INIT = .TRUE.
         DO 100 IELB = 1, NELBLK
            IF ((.NOT. USESEL) .OR. (IELBST(IELB) .GT. 0)) THEN
               IF (ISTMN(IELB) .GT. 0) THEN
                  IF (INIT) THEN
                     IMIN = IELB
                     IMAX = IELB
                     INIT = .FALSE.
                  ELSE
                     IF (VALMN(IMIN) .GT. VALMN(IELB)) THEN
                        IMIN = IELB
                     END IF
                     IF (VALMX(IMAX) .LT. VALMX(IELB)) THEN
                        IMAX = IELB
                     END IF
                  END IF
               END IF
            END IF
  100    CONTINUE
      ELSE
         INIT = .FALSE.
         IMIN = 0
         IMAX = 0
      END IF

      IF (INIT) THEN
         VALMIN = 0.0
         NUMMIN = 0
         DO 110 I = 1, 3
            XYZMIN(I) = 0.0
  110    CONTINUE
         ISTMIN = 0
         VALMAX = 0.0
         NUMMAX = 0
         DO 120 I = 1, 3
            XYZMAX(I) = 0.0
  120    CONTINUE
         ISTMAX = 0
      ELSE
         VALMIN = VALMN(IMIN)
         NUMMIN = NUMMN(IMIN)
         DO 130 I = 1, 3
            XYZMIN(I) = XYZMN(I,IMIN)
  130    CONTINUE
         ISTMIN = ISTMN(IMIN)
         VALMAX = VALMX(IMAX)
         NUMMAX = NUMMX(IMAX)
         DO 140 I = 1, 3
            XYZMAX(I) = XYZMX(I,IMAX)
  140    CONTINUE
         ISTMAX = ISTMX(IMAX)
      END IF

      RETURN
      END
