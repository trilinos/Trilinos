C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SPSEGM (NNE, ISEGEL, ISTART, IEND)
C=======================================================================

C   --*** SPSEGM *** (SPLOT) Find a segment of defined values
C   --   Written by Amy Gilkey - revised 11/05/87
C   --
C   --SPSEGM finds the starting and ending indices of a segment of
C   --defined values.  The segment is defined as having consecutive indices.
C   --
C   --Parameters:
C   --   NNE - IN - the number of defined elements
C   --   ISEGEL - IN - the indices of the defined elements
C   --   ISTART - OUT - the starting index of the next segment
C   --   IEND - IN/OUT - the ending index of the last segment on input;
C   --      the ending index of the next segment on output

      INTEGER ISEGEL(*)

      ISTART = IEND + 1

      DO 100 IEND = ISTART+1, NNE
         IF (ISEGEL(IEND-1)+1 .NE. ISEGEL(IEND)) GOTO 110
  100 CONTINUE
  110 CONTINUE
      IEND = IEND - 1

      RETURN
      END
