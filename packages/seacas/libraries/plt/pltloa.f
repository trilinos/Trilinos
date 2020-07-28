C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTLOA(LINE1,NUM,TYPE)
      CHARACTER*10 LINE
      CHARACTER*(*) LINE1
      INTEGER TYPE

      LINE1 = ' '
      LINE = ' '
      IF (TYPE.EQ.1) THEN
         IF (NUM.GE.0) THEN
            WRITE (LINE,10,ERR=20) INT(10.**NUM)

   10       FORMAT (I10)

   20       DO 2680 J = 1,10
               IF (LINE(J:J).NE.' ') THEN
                  GO TO 2690

               END IF

 2680       CONTINUE
 2690       CONTINUE
            LINE1 = LINE(J:)
            RETURN

         END IF

         LINE1(1:1) = '.'
         DO 2700 I = 1,ABS(NUM) - 1
            LINE1(I+1:I+1) = '0'
 2700    CONTINUE
         LINE1(I+1:I+1) = '1'

      ELSE
         WRITE (LINE,10,ERR=40) NUM
   40    DO 2720 J = 1,10
            IF (LINE(J:J).NE.' ') THEN
               GO TO 2730

            END IF

 2720    CONTINUE
 2730    CONTINUE
         LINE1 = LINE(J:)
      END IF

      RETURN

      END
