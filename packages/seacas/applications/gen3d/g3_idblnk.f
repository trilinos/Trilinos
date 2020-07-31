C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      INTEGER FUNCTION IDBLNK (IELBLK, IEL, IXELB, NUMLNK)
C=======================================================================

C   --*** IDBLNK *** (EXOLIB) Return link index
C   --   Written by Amy Gilkey - revised 12/03/87
C   --
C   --IDBLNK returns the link index of the element.
C   --
C   --Parameters:
C   --   IELBLK - IN/OUT - the element block number for IEL;
C   --      set if IELBLK <= 0
C   --   IEL - IN - the element number; assumed first in block if <= 0
C   --   IXELB - IN - the cumulative element counts by element block
C   --   NUMLNK - IN - the number of nodes per element

      INTEGER IXELB(0:*)
      INTEGER NUMLNK(*)

      IDBLNK = 1
      IF ((IELBLK .LE. 0) .AND. (IEL .LE. 0)) RETURN

      IF (IELBLK .LE. 0) THEN
         IELB = 0
   10    CONTINUE
         IF (IEL .GT. IXELB(IELB)) THEN
            IELB = IELB + 1
            GOTO 10
         END IF
         IELBLK = IELB
      END IF

      IDBLNK = 1
      DO 20 IELB = 1, IELBLK-1
         NEL = IXELB(IELB) - IXELB(IELB-1)
         IDBLNK = IDBLNK + NEL * NUMLNK(IELB)
   20 CONTINUE
      IF (IEL .GT. 0) THEN
         NEL = IEL - IXELB(IELBLK-1) - 1
         IDBLNK = IDBLNK + NEL * NUMLNK(IELBLK)
      END IF

      RETURN
      END
