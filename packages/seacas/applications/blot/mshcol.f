C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MSHCOL (FNCT, IELB, MLNTYP, WIDLIN, BLKCOL,
     &   IDELB)
C=======================================================================

C   --*** MSHCOL *** (MESH) Set line color, type, and width
C   --   Modified by John Glick - 1/17/89
C   --   Written by Amy Gilkey - revised 03/29/88
C   --
C   --MSHCOL sets the line color, type, and width.  The lines have been
C   --divided into parts as follows:
C   --   -1) Mesh boundary
C   --    0) Element block boundary
C   --    n) Interior within element block 'n'
C   --
C   --Be sure to call this routine to reset the line type and line width
C   --to the default before exiting the mesh-drawing routine.
C   --
C   --Parameters:
C   --   FNCT - IN - the type of lines being drawn:
C   --      'RESET'    - reset line type, etc.
C   --      'INTERIOR' - the interior mesh lines (small)
C   --      'BOUNDARY' - mesh boundaries defined by hidden lines
C   --      'LINE'     - normal mesh lines
C   --      'DEBUG'    - lines for debugging (blue)
C   --   IELB - IN - the line "part" (as above)
C   --   MLNTYP - IN - the line type (and color) of lines (as in /MSHOPT/)
C   --   WIDLIN - IN - true iff mesh lines should be wide versus narrow
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

      PARAMETER (KWIDLN=2)

      include 'dbnums.blk'

      COMMON /LINTHC/ MSHBND, BLKBND, ELEBND, THKNSS
      REAL MSHBND, BLKBND, ELEBND
C      --      Line thickness specification for lines appearing
C      --      on mesh plots.  Specification is a real value in the
C      --      range 0. - 1000., with 0. being the thinest line and
C      --      1000. being the thickest.
C      -- MSHBND - Thickness of lines forming the mesh boundary.
C      -- BLKBND - Thickness of lines forming element block boundaries.
C      -- ELEBND - Thickness of lines forming element boundaries.
      REAL THKNSS(3)
C     --       Line thickness specifications for THICK, MEDIUM,
C     --       and THIN keywords.

      CHARACTER*(*) FNCT
      INTEGER MLNTYP(-1:1)
      LOGICAL WIDLIN
      INTEGER BLKCOL(0:NELBLK)
      INTEGER IDELB(*)

      LOGICAL PLTGTV, PLTSTV, LDUM

      LOGICAL FIRST
      SAVE FIRST

      DATA FIRST / .TRUE. /

      IF (FIRST) THEN
C      --Get the default line width
         LDUM = PLTGTV (KWIDLN, WIDLN)
         FIRST = .FALSE.
      END IF

      IF ((FNCT .EQ. 'RESET') .OR. (IELB .LT. -1)) THEN
C      --Reset to default color, line type, and line width
         ICOLOR = 0
         LTYP = 1
         WID = 160.
C         WID = 1.00
      ELSE IF (IELB .EQ. -1) THEN
         ICOLOR = 0
         IF (MLNTYP(-1) .LT. 0) ICOLOR = -1
         LTYP = IABS (MLNTYP(-1))
C         WID = 1.75
         WID = MSHBND
      ELSE IF (IELB .EQ. 0) THEN
         ICOLOR = 0
         IF (MLNTYP(0) .LT. 0) ICOLOR = -1
         LTYP = IABS (MLNTYP(0))
C         WID = 1.25
         WID = BLKBND
      ELSE
         ICOLOR = IELB
         IF (MLNTYP(1) .LT. 0) ICOLOR = -1
         LTYP = IABS (MLNTYP(1))
         IF (WIDLIN) THEN
C            WID = 1.00
            WID = ELEBND
         ELSE
C            WID = 0.50
            WID = ELEBND / 2
         END IF
      END IF

      IF (FNCT .EQ. 'INTERIOR') THEN
C      --Interior lines as thin lines
         WID = 0.5 * WID
      ELSE IF (FNCT .EQ. 'BOUNDARY') THEN
C      --Element lines on hidden line boundaries - line type as
C      --boundary lines, black => white as boundary
         IF (IELB .GT. 0) THEN
            IF ((ICOLOR .EQ. -1) .AND. (MLNTYP(0) .GT. 0)) ICOLOR = 0
            LTYP = IABS (MLNTYP(0))
         END IF
      ELSE IF (FNCT .EQ. 'LINE') THEN
C      --Normal mesh lines are ok
         CONTINUE
      ELSE IF (FNCT .EQ. 'DEBUG') THEN
C      --Lines for debugging are blue
         CALL GRCOLR (4)
      END IF

C   --Make line size of dotted lines smaller
      IF (LTYP .NE. 1) THEN
         IF ((IELB .LE. 0) .OR. WIDLIN) THEN
            WID = 0.75 * WID
         END IF
      END IF

C      IF (ICOLOR .GT. 0) ICOLOR = IDELB(IELB)
C ... Test by Sjaardema, should pass ielb?
      IF (ICOLOR .GT. 0) ICOLOR = IELB
      CALL UGRCOL (ICOLOR, BLKCOL)
      LDUM = PLTSTV (1, 1.0 * LTYP)
      LDUM = PLTSTV (KWIDLN, WID)

      RETURN
      END
