C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE GTXDUM (X, DUMMY, LEN)
C***********************************************************************

C  SUBROUTINE GTXDUM = GETS A REAL INTO A DUMMY CHARACTER STRING

C***********************************************************************

      CHARACTER*72 DUMMY

      DUMMY = ' '
      IF(X .LE. -10000.) THEN
         WRITE (DUMMY(1:10), 10000) X
         LEN = 10
      ELSEIF(X .LE. -1000.) THEN
         WRITE (DUMMY(1:6), 10010) X
         LEN = 6
      ELSEIF(X .LE. -100.) THEN
         WRITE (DUMMY(1:6), 10020) X
         LEN = 6
      ELSEIF(X .LE. -10.) THEN
         WRITE (DUMMY(1:6), 10030) X
         LEN = 6
      ELSEIF(X .LE. -1.) THEN
         WRITE (DUMMY(1:6), 10040) X
         LEN = 6
      ELSEIF(X .LT. 0.) THEN
         WRITE (DUMMY(1:6), 10050) X
         LEN = 6
      ELSEIF(X .GE. 10000.) THEN
         WRITE (DUMMY(1:10), 10060) X
         LEN = 9
      ELSEIF(X .GE. 1000.) THEN
         WRITE (DUMMY(1:5), 10070) X
         LEN = 5
      ELSEIF(X .GE. 100.) THEN
         WRITE (DUMMY(1:5), 10080) X
         LEN = 5
      ELSEIF(X .GE. 10.) THEN
         WRITE (DUMMY(1:5), 10090) X
         LEN = 5
      ELSEIF(X .GE. 1.) THEN
         WRITE (DUMMY(1:5), 10100) X
         LEN = 5
      ELSE
         WRITE (DUMMY(1:5), 10110) X
         LEN = 5
      ENDIF
      RETURN

10000 FORMAT (1PE10.3)
10010 FORMAT (F6.0)
10020 FORMAT (F6.1)
10030 FORMAT (F6.2)
10040 FORMAT (F6.3)
10050 FORMAT (F6.4)
10060 FORMAT (1PE10.3)
10070 FORMAT (F5.0)
10080 FORMAT (F5.1)
10090 FORMAT (F5.2)
10100 FORMAT (F5.3)
10110 FORMAT (F5.4)

      END
