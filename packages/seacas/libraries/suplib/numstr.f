C Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software.
C         
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C 
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.
C     * Neither the name of Sandia Corporation nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

C=======================================================================
      SUBROUTINE NUMSTR (NNUM, NSIG, RNUM, RSTR, LSTR)
C=======================================================================
C$Id: numstr.f,v 1.2 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: numstr.f,v $
CRevision 1.2  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:15:53  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:15:52  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:40  gdsjaar
c Initial revision
c 

C   --*** NUMSTR *** (STRLIB) Convert real numbers to strings
C   --   Written by Amy Gilkey - revised 12/11/87
C   --
C   --NUMSTR converts a set of real numbers into a consistent set of
C   --strings.  It will convert to engineering notation with all
C   --exponents the same, if possible.
C   --
C   --Parameters:
C   --   NNUM - IN - the number of real numbers in the set
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

      CHARACTER*20 BLANKS
      CHARACTER*10 SCRFMT
      CHARACTER*20 SCRSTR
      CHARACTER*20 TMPSTR
      CHARACTER*15 FFMT

C   --Convert all to E notation and find the minimum and maximum exponent
C   --   MINE and MAXE are the minimum and maximum exponents
C   --   ISIGN is the number of digits for the sign
C   --      (0 if all positive, 1 if any number negative)

      BLANKS = ' '

      WRITE (SCRFMT, 10000, IOSTAT=IDUM) NSIG+7, NSIG
10000  FORMAT ('(0PE', I2.2, '.', I2.2, ')')

      ISIGN = 0
      MINE  = 9999
      MINE2 = 9999
      MAXE  = -9999
      MAXES = MAXE
      DO 100 I = 1, NNUM
         IF (RNUM(I) .NE. 0.0) THEN
            WRITE (SCRSTR(1:NSIG+7), SCRFMT, IOSTAT=IDUM) RNUM(I)
            READ (SCRSTR(NSIG+5:NSIG+7), '(I3)', IOSTAT=IDUM) IE
            IF (MINE .GT. IE) MINE2 = MINE
            MINE = MIN (MINE, IE)
            MAXE = MAX (MAXE, IE)
            IF (RNUM(I) .LT. 0.0) THEN
               ISIGN = 1
               MAXES = MAX (MAXES, IE)
            END IF
         END IF
  100 CONTINUE

C   --Correct for one very small number (should be zero)

      IF ((MINE2 .LT. 1000) .AND. ((MINE2 - MINE) .GE. 6)) MINE = MINE2

C   --Handle all zero case

      IF (MINE .GT. MAXE) THEN
         MINE = 0
         MAXE = 0
         MAXES = 0
      END IF

C   --Determine the new exponent NEWEXP (use engineering notation)

      NEWEXP = IENGRX (MAXE, MINE)
      IF (ISIGN .EQ. 1) THEN
         IF (MAX (1, MAXE - NEWEXP) .GT. MAX (1, MAXES - NEWEXP))
     &      ISIGN = 0
      END IF

C   --Check if the numbers can all be sensibly converted to a common exponent

      IF (((MAXE - NEWEXP) .LE. 4)
     &   .AND. ((NEWEXP - MINE) .LE. 2)
     &   .AND. (-MINE .LT. (NSIG - MAXE))) THEN

C      --Determine the new F format
C      --   EXPDIV is the number to divide by to get the number
C      --      without an exponent
C      --   NWHOLE is the number of digits before the decimal
C      --   NFRAC is the number of digits after the decimal
C      --   NTOTAL is the total number of digits
C      --The new exponent is tagged on the end of the F-format number

         EXPDIV = 10.0 ** NEWEXP

         NWHOLE = MAX (1, MAXE - NEWEXP)
         NFRAC = MAX (0, MIN (NEWEXP - MINE + NSIG,
     &      NSIG - (MAXE - NEWEXP)))
         NTOTAL = ISIGN + NWHOLE + 1 + NFRAC
         IF (EXPDIV .NE. 0.0) THEN
            WRITE (FFMT, 10010, IOSTAT=IDUM) NTOTAL, NFRAC
10010        FORMAT ('(F', I2.2, '.', I2.2, ')')
         ELSE
            WRITE (FFMT, 10020, IOSTAT=IDUM) NTOTAL
10020        FORMAT ('(A', I2.2, 3X, ')')
         END IF

         IF (NEWEXP .EQ. 0) THEN
            LSTR = NTOTAL
         ELSE IF ((NEWEXP .LE. -10) .OR. (NEWEXP .GE. 10)) THEN
            WRITE (FFMT(8:15), 10030, IOSTAT=IDUM) NEWEXP
10030        FORMAT (',''E', SP, I3.2, ''')')
            LSTR = NTOTAL + 4
         ELSE
            WRITE (FFMT(8:15), 10040, IOSTAT=IDUM) NEWEXP
10040        FORMAT (',''E', SP, I2.1, ''')')
            LSTR = NTOTAL + 3
         END IF

C      --Convert all numbers to the new exponent by using the F format

         IF (EXPDIV .NE. 0.0) THEN
            DO 110 I = 1, NNUM
               WRITE (RSTR(I), FFMT, IOSTAT=IDUM) RNUM(I)/EXPDIV
  110       CONTINUE
         ELSE
            DO 120 I = 1, NNUM
               WRITE (RSTR(I), FFMT, IOSTAT=IDUM) '********************'
  120       CONTINUE
         END IF

      ELSE

C      --Do not try to use a common exponent, but use engineering notation;
C      --Algorithm as above

         LSTR = 0
         MINEXP = IENGRX (MINE, MINE)
         MAXEXP = IENGRX (MAXE, MAXE)

         DO 130 I = 1, NNUM
            WRITE (SCRSTR(1:NSIG+7), SCRFMT, IOSTAT=IDUM) RNUM(I)
            READ (SCRSTR(NSIG+5:NSIG+7), '(I3)', IOSTAT=IDUM) IE
            ISIGN = 0
            IF (RNUM(I) .LT. 0.0) ISIGN = 1

            NEWEXP = IENGRX (IE, IE)

            EXPDIV = 10.0 ** NEWEXP

            NWHOLE = MAX (1, IE - NEWEXP)
            NFRAC = MAX (0, MIN (NEWEXP - IE + NSIG,
     &         NSIG - (IE - NEWEXP)))
            IF ((RNUM(I) .EQ. 0.0) .AND. (MINE .GE. 0))
     &         NFRAC = NFRAC - 1
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
               WRITE (RSTR(I), FFMT, IOSTAT=IDUM) RNUM(I)/EXPDIV
            ELSE
               WRITE (RSTR(I), FFMT, IOSTAT=IDUM) '********************'
            END IF
  130    CONTINUE

C      --Adjust the strings so that they are right-justified at
C      --a common length

         DO 140 I = 1, NNUM
            IB = INDEX (RSTR(I)(:LSTR), ' ')
            IF (IB .GT. 0) THEN
               NB = LSTR - IB + 1
               TMPSTR = RSTR(I)(:IB-1)
               RSTR(I) = BLANKS(:NB) // TMPSTR
            END IF
  140    CONTINUE

      END IF

      RETURN
      END
