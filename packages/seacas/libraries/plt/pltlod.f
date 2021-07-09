C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTLOD(LINE1,J,NUM)
      CHARACTER*(*) LINE1
      CHARACTER*10 LINE

      LINE1 = ' '
      LINE = ' '
      IF (NUM.GE.0) THEN
         WRITE (LINE,10,ERR=20) INT(DBLE(J)*10.**NUM)

   10    FORMAT (I10)

   20    DO 2740 I = 1,10
            IF (LINE(I:I).NE.' ') THEN
               GO TO 2750

            END IF

 2740    CONTINUE
 2750    CONTINUE
         LINE1 = LINE(I:)
         RETURN

      END IF

      LINE1(1:1) = '.'
      DO 2760 I = 1,ABS(NUM) - 1
         LINE1(I+1:I+1) = '0'
 2760 CONTINUE
      LINE1(I+1:I+1) = CHAR(J+48)
      RETURN

      END
