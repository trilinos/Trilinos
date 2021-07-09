C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE CHKRNG (IARRAY, NITEMS, MAXVAL, NZERO, NERR)
C=======================================================================

C   --*** CHKRNG *** (EXPLORE) Check count of items
C   --
C   --CHKRNG checks that each item in an array of items is within the
C   --specified range.
C   --
C   --Parameters:
C   --   IARRAY - IN  - the array
C   --   NITEMS - IN  - the number of items in IARRAY
C   --   MAXVAL - IN  - the maximum range value, minimum = 1
C   --   NZERO  - OUT - the number of 'zero' out-of-range values
C   --   NERR   - OUT - the number of out-of-range values

      INTEGER IARRAY(*)

      NERR  = 0
      NZERO = 0
      DO 100 I = 1, NITEMS
         N = IARRAY(I)
         IF (N .LT. 1)      NZERO = NZERO + 1
         IF (N .GT. MAXVAL) NERR  = NERR + 1
  100 CONTINUE

      RETURN
      END
