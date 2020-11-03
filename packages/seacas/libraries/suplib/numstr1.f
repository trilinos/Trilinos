C     Copyright(C) 1999-2020 National Technology & Engineering Solutions
C     of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C     NTESS, the U.S. Government retains certain rights in this software.
C
C     See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE NUMSTR1 (NSIG, RNUM, RSTR, LSTR)
C=======================================================================

C     --*** NUMSTR *** (STRLIB) Convert real numbers to strings
C     --
C     --NUMSTR converts a set of real numbers into a consistent set of
C     --strings.  It will convert to engineering notation with all
C     --exponents the same, if possible.
C     --
C     --Parameters:
C     --   NSIG - IN - the maximum number of significant digits, max of 8
C     --   RNUM - IN - the array of real numbers to be converted
C     --   RSTR - OUT - the set of real number strings
C     --   LSTR - OUT - the maximum length of the number strings

C     --Routines Called:
C     --   IENGRX - (included) Get engineering notation exponent

      INTEGER NSIG
      REAL RNUM
      CHARACTER*(*) RSTR
      INTEGER LSTR

      CHARACTER*20 BLANKS
      CHARACTER*10 SCRFMT
      CHARACTER*20 SCRSTR
      CHARACTER*20 TMPSTR
      CHARACTER*15 FFMT

C     --Convert all to E notation and find the minimum and maximum exponent
C     --   MINE and MAXE are the minimum and maximum exponents
C     --   ISIGN is the number of digits for the sign
C     --      (0 if all positive, 1 if any number negative)

      BLANKS = ' '

      WRITE (SCRFMT, 10000, IOSTAT=IDUM) NSIG+7, NSIG
10000 FORMAT ('(0PE', I2.2, '.', I2.2, ')')

      ISIGN = 0
      MINE  = 9999
      MINE2 = 9999
      MAXE  = -9999
      MAXES = MAXE
      IF (RNUM .NE. 0.0) THEN
         WRITE (SCRSTR(1:NSIG+7), SCRFMT, IOSTAT=IDUM) RNUM
         READ (SCRSTR(NSIG+5:NSIG+7), '(I3)', IOSTAT=IDUM) IE
         IF (MINE .GT. IE) MINE2 = MINE
         MINE = MIN (MINE, IE)
         MAXE = MAX (MAXE, IE)
         IF (RNUM .LT. 0.0) THEN
            ISIGN = 1
            MAXES = MAX (MAXES, IE)
         END IF
      END IF

C     --Correct for one very small number (should be zero)

      IF ((MINE2 .LT. 1000) .AND. ((MINE2 - MINE) .GE. 6)) MINE = MINE2

C     --Handle all zero case

      IF (MINE .GT. MAXE) THEN
         MINE = 0
         MAXE = 0
         MAXES = 0
      END IF

C     --Determine the new exponent NEWEXP (use engineering notation)

      NEWEXP = IENGRX (MAXE, MINE)
      IF (ISIGN .EQ. 1) THEN
         IF (MAX (1, MAXE - NEWEXP) .GT. MAX (1, MAXES - NEWEXP))
     &        ISIGN = 0
      END IF

C     --Check if the numbers can all be sensibly converted to a common exponent

      IF (((MAXE - NEWEXP) .LE. 4)
     &     .AND. ((NEWEXP - MINE) .LE. 2)
     &     .AND. (-MINE .LT. (NSIG - MAXE))) THEN

C     --Determine the new F format
C     --   EXPDIV is the number to divide by to get the number
C     --      without an exponent
C     --   NWHOLE is the number of digits before the decimal
C     --   NFRAC is the number of digits after the decimal
C     --   NTOTAL is the total number of digits
C     --The new exponent is tagged on the end of the F-format number

         EXPDIV = 10.0 ** NEWEXP

         NWHOLE = MAX (1, MAXE - NEWEXP)
         NFRAC = MAX (0, MIN (NEWEXP - MINE + NSIG,
     &        NSIG - (MAXE - NEWEXP)))
         NTOTAL = ISIGN + NWHOLE + 1 + NFRAC
         IF (EXPDIV .NE. 0.0) THEN
            WRITE (FFMT, 10010, IOSTAT=IDUM) NTOTAL, NFRAC
10010       FORMAT ('(F', I2.2, '.', I2.2, ')')
         ELSE
            WRITE (FFMT, 10020, IOSTAT=IDUM) NTOTAL
10020       FORMAT ('(A', I2.2, 3X, ')')
         END IF

         IF (NEWEXP .EQ. 0) THEN
            LSTR = NTOTAL
         ELSE IF ((NEWEXP .LE. -10) .OR. (NEWEXP .GE. 10)) THEN
            WRITE (FFMT(8:15), 10030, IOSTAT=IDUM) NEWEXP
10030       FORMAT (',''E', SP, I3.2, ''')')
            LSTR = NTOTAL + 4
         ELSE
            WRITE (FFMT(8:15), 10040, IOSTAT=IDUM) NEWEXP
10040       FORMAT (',''E', SP, I2.1, ''')')
            LSTR = NTOTAL + 3
         END IF

C     --Convert all numbers to the new exponent by using the F format

         IF (EXPDIV .NE. 0.0) THEN
            WRITE (RSTR, FFMT, IOSTAT=IDUM) RNUM/EXPDIV
            if (rstr(:1) .eq. '*') then
C     ... Roundoff occurred. Adjust format and try again...
               IF (EXPDIV .NE. 0.0) THEN
                  WRITE (FFMT(:7), 10010, IOSTAT=IDUM) NTOTAL,
     $                 NFRAC-1
                  WRITE (RSTR, FFMT, IOSTAT=IDUM) RNUM/EXPDIV
               end if
            end if
         ELSE
            WRITE (RSTR, FFMT, IOSTAT=IDUM) '********************'
         END IF

      ELSE

C     --Do not try to use a common exponent, but use engineering notation;
C     --Algorithm as above

         LSTR = 0
         MINEXP = IENGRX (MINE, MINE)
         MAXEXP = IENGRX (MAXE, MAXE)

         WRITE (SCRSTR(1:NSIG+7), SCRFMT, IOSTAT=IDUM) RNUM
         READ (SCRSTR(NSIG+5:NSIG+7), '(I3)', IOSTAT=IDUM) IE
         ISIGN = 0
         IF (RNUM .LT. 0.0) ISIGN = 1

         NEWEXP = IENGRX (IE, IE)

         EXPDIV = 10.0 ** NEWEXP

         NWHOLE = MAX (1, IE - NEWEXP)
         NFRAC = MAX (0, MIN (NEWEXP - IE + NSIG,
     &        NSIG - (IE - NEWEXP)))
         IF ((RNUM .EQ. 0.0) .AND. (MINE .GE. 0))
     &        NFRAC = NFRAC - 1
         NTOTAL = ISIGN + NWHOLE + 1 + NFRAC
         IF (EXPDIV .NE. 0.0) THEN
            WRITE (FFMT, 10010, IOSTAT=IDUM) NTOTAL, NFRAC
         ELSE
            WRITE (FFMT, 10020, IOSTAT=IDUM) NTOTAL
         END IF

         IF ((MINEXP .LE. -10) .OR. (MAXEXP .GE. 10)) THEN
            WRITE (FFMT(8:15), 10030, IOSTAT=IDUM) NEWEXP
            LSTR = MAX (LSTR, NTOTAL + 4)
         ELSE
            WRITE (FFMT(8:15), 10040, IOSTAT=IDUM) NEWEXP
            LSTR = MAX (LSTR, NTOTAL + 3)
         END IF

         IF (EXPDIV .NE. 0.0) THEN
            WRITE (RSTR, FFMT, IOSTAT=IDUM) RNUM/EXPDIV
            if (rstr(:1) .eq. '*') then
C     ... Roundoff occurred. Adjust format and try again...
               IF (EXPDIV .NE. 0.0) THEN
                  WRITE (FFMT(:7), 10010, IOSTAT=IDUM) NTOTAL,
     $                 NFRAC-1
                  WRITE (RSTR, FFMT, IOSTAT=IDUM) RNUM/EXPDIV
               end if
            end if
         ELSE
            WRITE (RSTR, FFMT, IOSTAT=IDUM) '********************'
         END IF

C     --Adjust the strings so that they are right-justified at
C     --a common length

         IB = INDEX (RSTR(:LSTR), ' ')
         IF (IB .GT. 0) THEN
            NB = LSTR - IB + 1
            TMPSTR = RSTR(:IB-1)
            RSTR = BLANKS(:NB) // TMPSTR
         END IF

      END IF

      RETURN
      END
