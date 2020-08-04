C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE LOWSTR (LCSTR, UCSTR)
C=======================================================================

C   --*** LOWSTR *** (STRLIB) Convert string to lower-case
C   --   Written by Amy Gilkey - revised 08/06/87
C   --
C   --LOWSTR converts the passed string to lower-case letters.
C   --
C   --Parameters:
C   --   LCSTR - OUT - the returned lower-case string
C   --   UCSTR - IN - the input string

      CHARACTER*(*) LCSTR, UCSTR

      CHARACTER*26 UPPER, LOWER
      SAVE UPPER, LOWER

      DATA UPPER / 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' /
      DATA LOWER / 'abcdefghijklmnopqrstuvwxyz' /

      LCSTR = UCSTR
      DO 10 I = 1, LENSTR (LCSTR)
         K = INDEX (UPPER, LCSTR(I:I))
         IF (K .GE. 1) LCSTR(I:I) = LOWER(K:K)
   10 CONTINUE

      RETURN
      END
