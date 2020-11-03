C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SCAINI (NELB, NVAR, ISTMN)
C=======================================================================

C   --*** SCAINI *** (BLOT) Initialize variable min/max
C   --   Written by Amy Gilkey - revised 04/01/88
C   --
C   --SCAINI initializes the variable minimums and maximums to uncalculated.
C   --ISTMN is the "state" of the minimums and maximums: <0 if not calculated,
C   --0 if no elements in element block, +n if calculated
C   --
C   --Parameters:
C   --   NELB - IN - the number of element blocks (NELBLK for element,
C   --      else 0)
C   --   NVAR - IN - the number of variables of this type
C   --   ISTMN - OUT - initialized to -999

      INTEGER ISTMN(0:NELB,*)

      DO 110 IXVAR = 1, NVAR
         DO 100 IELB = 0, NELB
            ISTMN(IELB,IXVAR) = -999
  100    CONTINUE
  110 CONTINUE

      RETURN
      END
