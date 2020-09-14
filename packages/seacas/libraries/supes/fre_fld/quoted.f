C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
      SUBROUTINE QUOTED ( LINE, ILEFT, IRIGHT )

      CHARACTER*(*) LINE

      CALL STRIPB ( LINE, ILEFT, IRIGHT )

C     The first character is required to be a quote, so remove it.

      LINE(1:1) = ' '
      ILEFT = 2
      IBEGIN = 2

C     Begin loop looking for more quotes.  There should be at least 1 more.

 100  CONTINUE
         IQUOTE = INDEX ( LINE(IBEGIN:IRIGHT), '''' )

C        Has the quote ended within this record?

         IF ( IQUOTE .EQ. 0 ) THEN
            IF ( ILEFT .GT. IRIGHT ) IRIGHT = 0
            RETURN
         END IF

         IQUOTE = IQUOTE + IBEGIN - 1
         IF ( IQUOTE .EQ. IRIGHT ) THEN

C           The quote is at the end of the record.

            LINE(IRIGHT:IRIGHT) = ' '
            IRIGHT = IQUOTE - 1
            IF ( ILEFT .GT. IRIGHT ) IRIGHT = 0
            RETURN
         END IF

C        The quote is internal -- check for double quote.

         IF ( LINE(IQUOTE+1:IQUOTE+1) .NE. '''' ) THEN

C           The quote is single, thus ending the quoted string.

            LINE(IQUOTE:IQUOTE) = ' '
            IRIGHT = IQUOTE - 1
            IF ( ILEFT .GT. IRIGHT ) IRIGHT = 0
            RETURN
         END IF

C        The quote is a double quote.  Remove the repeated quote and loop.

         DO 10 I = IQUOTE, ILEFT, -1
            LINE(I+1:I+1) = LINE(I:I)
 10      CONTINUE
         LINE(ILEFT:ILEFT) = ' '
         ILEFT = ILEFT + 1
         IBEGIN = IQUOTE + 2
      GO TO 100
      END
