C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTESC(TEXT,I,ESC)
      CHARACTER*(*) TEXT,ESC

      IP = I + 1
      LT = LEN(TEXT)
      IF (TEXT(IP:IP).EQ.'^') THEN
         ESC = '^'
         I = I + 1
         RETURN

      ELSE IF (TEXT(IP:IP).EQ.'-') THEN
         ESC = '-'
         I = I + 1
         RETURN

      ELSE IF (TEXT(IP:IP).EQ.'_') THEN
         ESC = '_'
         I = I + 1
         RETURN

      END IF

      DO 2090 J = IP,LT
         IF (TEXT(J:J).EQ.' ' .OR. TEXT(J:J).EQ. CHAR(92)) THEN
            GO TO 2100

         END IF

 2090 CONTINUE
 2100 CONTINUE
      ESC = TEXT(IP:J-1)
      I = J
      IF (I.GE.LT) THEN
         RETURN

      END IF

      IF (TEXT(J:J).EQ. CHAR(92)) THEN
         I = J - 1
      END IF

      RETURN

      END
