C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE ENGSTR (NNUM, NSIG, RNUM, RSTR, LSTR)
C=======================================================================

C   --*** ENGSTR *** (STRLIB) Convert real numbers to strings
C   --   Written by Amy Gilkey - revised 02/14/86
C   --
C   --ENGSTR converts a set of real numbers into engineering notation
C   --strings.  The real numbers are not compared in any way to get a
C   --consistent set of exponents, etc.
C   --
C   --Parameters:
C   --   NNUM - IN - the number of numbers in the array set
C   --   NSIG - IN - the maximum number of significant digits, max of 8
C   --   RNUM - IN - the array of real numbers to be converted
C   --   RSTR - OUT - the set of real number strings
C   --   LSTR - OUT - the maximum length of the number strings

C   --Routines Called:
C   --   IENGRX - (included) Get engineering notation exponent

      INTEGER NNUM
      INTEGER NSIG
      REAL RNUM(*)
      CHARACTER*(*) RSTR(*)
      INTEGER LSTR

      CHARACTER*10 SCRFMT
      CHARACTER*20 SCRSTR
      CHARACTER*15 FFMT

C   --Do not try to use a common exponent, but use engineering notation

      WRITE (SCRFMT, 10, IOSTAT=IDUM) NSIG+7, NSIG
   10 FORMAT ('(0PE', I2.2, '.', I2.2, ')')

      LSTR = 0
      DO 40 I = 1, NNUM
         WRITE (SCRSTR(1:NSIG+7), SCRFMT, IOSTAT=IDUM) RNUM(I)
         READ (SCRSTR(NSIG+5:NSIG+7), '(I3)', IOSTAT=IDUM) IE
         ISIGN = 0
         IF (RNUM(I) .LT. 0.0) ISIGN = 1

         NEWEXP = IENGRX (IE, IE)

         EXPDIV = 10.0 ** NEWEXP

         NWHOLE = MAX (0, IE - NEWEXP)
         NFRAC = MAX (0, MIN (NEWEXP - IE + NSIG,
     &      NSIG - (IE - NEWEXP)))
         NTOTAL = ISIGN + NWHOLE + 1 + NFRAC

         WRITE (FFMT, 20, IOSTAT=IDUM) NTOTAL, NFRAC
   20    FORMAT ('(F', I2.2, '.', I2.2, ')')
         WRITE (FFMT(8:15), 30, IOSTAT=IDUM) NEWEXP
   30    FORMAT (',''E', SP, I3.2, ''')')
         LSTR = MAX (LSTR, NTOTAL+4)

         WRITE (RSTR(I), FFMT, IOSTAT=IDUM) RNUM(I)/EXPDIV
   40 CONTINUE

      RETURN
      END
