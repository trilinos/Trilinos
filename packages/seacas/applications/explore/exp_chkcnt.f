C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE CHKCNT (ICHECK, LENCHK, NZERO, NMULT)
C=======================================================================

C   --*** CHKCNT *** (EXPLORE) Check count of items
C   --
C   --CHKCNT checks an array of counted items for counts of zero or counts
C   --greater than one.
C   --
C   --Parameters:
C   --   ICHECK - IN - the counts
C   --   LENCHK - IN - the length of ICHECK
C   --   NZERO - OUT - the number of zero counts
C   --   NMULT - OUT - the number of counts greater than one

      INTEGER ICHECK(*)

      NZERO = 0
      NMULT = 0
      DO 100 I = 1, LENCHK
         IF (ICHECK(I) .LE. 0) NZERO = NZERO + 1
         IF (ICHECK(I) .GT. 1) NMULT = NMULT + 1
  100 CONTINUE

      RETURN
      END
