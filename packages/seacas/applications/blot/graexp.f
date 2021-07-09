C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE GRAEXP (AXTYP, NNUM, RNUM, ATIC)
C=======================================================================

C   --*** GRAEXP *** (GRPLIB) Choose axis exponent (PLT)
C   --   Written by Amy Gilkey - revised 02/20/87
C   --
C   --GRAEXP chooses an exponent from the axes limits and sets this
C   --value and the number of decimal digits for the given axis.
C   --It will convert to engineering notation if possible.
C   --
C   --Parameters:
C   --   AXTYP - IN - the axis (' ' for both, 'X' for X, 'Y' for Y)
C   --   NNUM - IN - the number of numbers in the array set;
C   --      0 for no exponent
C   --   RNUM - IN - the array of axes limits
C   --   ATIC - IN - the axis tick mark interval

C   --Routines Called:
C   --   PLTSTG - (PLTLIB) Set graph parameter
C   --      13, 14 = (KXEXP, KYEXP) X, Y axis exponent
C   --      19, 20 = (KXNDIG, KYNDIG) X, Y axis number of decimal digits

      PARAMETER (KXEXP=13, KYEXP=14, KXNDIG=19, KYNDIG=20)

      CHARACTER AXTYP
      INTEGER NNUM
      REAL RNUM(*)
      REAL ATIC

      LOGICAL LDUM, PLTSTG
      CHARACTER*12 SCRSTR

      IF (NNUM .GT. 0) THEN

C      --Convert all to E notation and find the minimum and maximum exponent
C      --   MINE and MAXE are the minimum and maximum exponents

         MINE = 999
         MAXE = -999
         DO 100 I = 1, NNUM
            IF (RNUM(I) .NE. 0.0) THEN
               WRITE (SCRSTR(1:8), '(0PE8.1)') RNUM(I)
               READ (SCRSTR(6:8), '(I3)') IE
               MINE = MIN (MINE, IE)
               MAXE = MAX (MAXE, IE)
            END IF
  100    CONTINUE

         IF (MAXE .GT. 0) THEN
            IEXP = INT ((MAXE - 1) / 3) * 3
         ELSE
            IEXP = INT ((MAXE - 2) / 3) * 3
         END IF

         WRITE (SCRSTR(1:12), '(0PE12.5)') ATIC
         READ (SCRSTR(10:12), '(I3)') ITICE
         DO 110 I = 8, 3, -1
            IF (SCRSTR(I:I) .NE. '0') GOTO 120
  110    CONTINUE
  120    CONTINUE
         NSIG = MAX (2, I-3 + MAXE - ITICE)
         NDIG = MAX (0, IEXP - MAXE + NSIG)

      ELSE

C      --Set exponent and number of digits to zero

         IEXP = 0
         NDIG = 0
      END IF

      IF (AXTYP .NE. 'Y') THEN
         LDUM = PLTSTG (KXEXP, REAL(IEXP))
         LDUM = PLTSTG (KXNDIG, REAL(NDIG))
      END IF
      IF (AXTYP .NE. 'X') THEN
         LDUM = PLTSTG (KYEXP, REAL(IEXP))
         LDUM = PLTSTG (KYNDIG, REAL(NDIG))
      END IF

      RETURN
      END
