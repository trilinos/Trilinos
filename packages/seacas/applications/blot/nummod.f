C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      INTEGER FUNCTION NUMMOD (MODE, MODSUB, MATMOD, MATSUB)
C=======================================================================

C   --*** NUMMOD *** (MESH) Count the number of matching modes
C   --   Written by Amy Gilkey - revised 10/05/87
C   --
C   --NUMMOD returns the number of defined modes that match both mode
C   --and sub-mode (if given).
C   --
C   --Parameters:
C   --   MODE - the modes
C   --   MODSUB - the sub-modes
C   --   MATMOD - IN - the mode to match
C   --   MATSUB - IN - the sub-mode to match (if not ' ')

      CHARACTER*(*) MODE(4), MODSUB(4), MATMOD, MATSUB

      NUMMOD = 0

      IF (MATSUB .EQ. ' ') THEN
         DO 100 IVIEW = 1, 4
            IF (MODE(IVIEW) .EQ. MATMOD) NUMMOD = NUMMOD + 1
  100    CONTINUE
      ELSE
         DO 110 IVIEW = 1, 4
            IF ((MODE(IVIEW) .EQ. MATMOD)
     &         .AND. (MODSUB(IVIEW) .EQ. MATSUB))
     &         NUMMOD = NUMMOD + 1
  110    CONTINUE
      END IF

      RETURN
      END
