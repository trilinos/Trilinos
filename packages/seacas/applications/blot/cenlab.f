C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CENLAB (LAB, NSNUM, SNUM, LMAX, RETLAB)
C=======================================================================

C   --*** CENLAB *** (BLOT) Center label over number strings
C   --   Written by Amy Gilkey - revised 02/06/85
C   --
C   --CENLAB centers the label over the number strings, enlarging either
C   --either as needed.
C   --
C   --Parameters:
C   --   LAB - IN - the left-justified label
C   --   NSNUM - IN - the number of number strings
C   --   SNUM - IN/OUT - the number strings, returned centered
C   --   LMAX - IN/OUT - the input number string lengths, returned centered
C   --      string lengths
C   --   RETLAB - OUT - the centered label

      CHARACTER*(*) RETLAB, LAB, SNUM(NSNUM)

      CHARACTER*20 TMPSTR
      CHARACTER*20 BLANKS

      DATA BLANKS / '                    ' /

      LL = LENSTR (LAB)

      I = INT (IABS(LMAX-LL) / 2)
      IF (I .EQ. 0) THEN
         RETLAB = LAB
      ELSE IF (LL .LT. LMAX) THEN
         RETLAB = BLANKS(:I) // LAB
      ELSE
         DO 100 J = 1, NSNUM
            TMPSTR = SNUM(J)
            SNUM(J) = BLANKS(:I) // TMPSTR
  100    CONTINUE
      END IF

      LMAX = MAX (LMAX, LL)

      RETURN
      END
