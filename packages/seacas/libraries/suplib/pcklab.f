C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

      SUBROUTINE PCKLAB (CURVE,ACURVE,NCURVE)

C     THIS SUBROUTINE PACKS A 8 CHARACTER WORD AND A 'I8' MAXIMUM
C     INTEGER INTO A 16 CHARACTER WORD BY REMOVING
C     INCLUDED BLANKS -- USED TO CREATE VARIABLE 'CURVE' FOR PLTONE

      CHARACTER*8 ACURVE
      CHARACTER*16 CURVE
      CHARACTER*1 BLANK
      DATA BLANK/' '/

      IF (NCURVE .NE. 0) THEN
         WRITE (CURVE, 10) ACURVE, NCURVE
   10    FORMAT (A8,I8)
       ELSE
         CURVE(1:8)= ACURVE(1:8)
         CURVE(9:16)= '        '
       END IF
      L = 0
      DO 20 J=1,16
         IF (CURVE(J:J).NE.BLANK) THEN
            L = L + 1
            CURVE(L:L)= CURVE(J:J)
          END IF
   20    CONTINUE
      L1 = L + 1
      DO 30 J=L1,16
         CURVE(J:J)= BLANK
   30    CONTINUE
      RETURN

      END
