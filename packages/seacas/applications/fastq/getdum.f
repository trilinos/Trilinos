C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE GETDUM (I, DUMMY, LEN)
C***********************************************************************

C  SUBROUTINE GETDUM = GETS AN INTEGER INTO A DUMMY CHARACTER STRING

C***********************************************************************

      CHARACTER*72 DUMMY

      DUMMY = ' '
      IF (I .LT. -9999) THEN
         WRITE(DUMMY(1:6),10050)I
         LEN = 6
      ELSEIF (I .LT. -999) THEN
         WRITE (DUMMY(1:5), 10040) I
         LEN = 5
      ELSEIF (I .LT. -99) THEN
         WRITE (DUMMY(1:4), 10030) I
         LEN = 4
      ELSEIF (I .LT. -9) THEN
         WRITE (DUMMY(1:3), 10020) I
         LEN = 3
      ELSEIF (I .LT. 0) THEN
         WRITE (DUMMY(1:2), 10010) I
         LEN = 6
      ELSEIF(I .LT. 10) THEN
         WRITE (DUMMY(1:1), 10000) I
         LEN = 1
      ELSEIF(I.LT.100)THEN
         WRITE (DUMMY(1:2), 10010) I
         LEN = 2
      ELSEIF (I .LT. 1000) THEN
         WRITE (DUMMY(1:3), 10020) I
         LEN = 3
      ELSEIF (I .LT. 10000) THEN
         WRITE (DUMMY(1:4), 10030) I
         LEN = 4
      ELSE
         WRITE (DUMMY(1:5), 10040) I
         LEN = 5
      ENDIF
      RETURN

10000 FORMAT (I1)
10010 FORMAT (I2)
10020 FORMAT (I3)
10030 FORMAT (I4)
10040 FORMAT (I5)
10050 FORMAT (I6)

      END
