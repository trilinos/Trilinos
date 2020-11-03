C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE NDSTAT (NODE, LXN, ANGLE, JSTAT)
C***********************************************************************

C  SUBROUTINE NDSTAT = DETERMINES THE MOST APPROPRIATE STATUS OF A
C                      GIVEN NODE

C***********************************************************************

      DIMENSION LXN(4)

C  AN UNDECIDED NODE HAS BEEN FOUND - TEST ANGLE AND CONNECTIVITY

C  THE NODE IS ON THE BOUNDARY
      IF (LXN (2) .LT. 0) THEN

C IF THE NODE HAS LESS THAN FOUR LINES ATTACHED
C  CUTOFFS ARE:    0 TO 135 DEGREES = ROW END
C                135 TO 225 DEGREES = ROW SIDE
C                225 TO 290 DEGREES = ROW CORNER
C                OVER 290 DEGREES   = ROW REVERSAL

         IF (LXN (4) .LE. 0) THEN
            IF (ANGLE .LT. 2.3561945) THEN
               JSTAT = 1
            ELSE IF (ANGLE .LT. 3.9269908) THEN
               JSTAT = 3
            ELSE IF (ANGLE .LT. 5.0614548) THEN
               JSTAT = 5
            ELSE
               JSTAT = 7
            ENDIF

C IF THE NODE HAS FOUR LINES ATTACHED
C  CUTOFFS ARE:    0 TO 110 DEGREES = ROW END
C                110 TO 225 DEGREES = ROW SIDE
C                OVER 225 DEGREES   = ROW CORNER (NEARLY IMPOSSIBLE)

         ELSE
            IF (ANGLE .LT. 1.9198622) THEN
               JSTAT = 1
            ELSE IF (ANGLE .LT. 3.9269908) THEN
               JSTAT = 3
            ELSE
               JSTAT = 5
            ENDIF
         ENDIF

C  THE NODE IS NOT ON THE BOUNDARY - CUTOFFS ARE ADJUSTED BASED
C  ON THE CONNECTIVITY AND THE ANGLE

      ELSE

C  ONLY TWO LINES ARE ATTACHED - LEAN TOWARDS A ROW CORNER NODE
C  OR A ROW END NODE

         IF (LXN(3) .EQ. 0) THEN

C  CUTOFFS ARE:    0 TO 135 DEGREES = ROW END
C                135 TO 210 DEGREES = ROW SIDE
C                210 TO 320 DEGREES = ROW CORNER
C                OVER 320 DEGREES   = ROW REVERSAL

            IF (ANGLE .LT. 2.3561945) THEN
               JSTAT = 1
            ELSE IF (ANGLE .LT. 3.6651914) THEN
               JSTAT = 3
            ELSE IF (ANGLE .LT. 5.5850536) THEN
               JSTAT = 5
            ELSE
               JSTAT = 7
            ENDIF

C  THREE LINES ARE ATTACHED - LEAN TOWARDS A ROW SIDE

         ELSEIF (LXN(4) .EQ. 0) THEN

C  CUTOFFS ARE:    0 TO 110 DEGREES = ROW END
C                110 TO 240 DEGREES = ROW SIDE
C                240 TO 320 DEGREES = ROW CORNER
C                OVER 320 DEGREES   = ROW REVERSAL (REALLY IMPOSSIBLE)

            IF (ANGLE .LT. 1.9198622) THEN
               JSTAT = 1
            ELSE IF (ANGLE .LT. 4.1887902) THEN
               JSTAT = 3
            ELSE IF (ANGLE .LT. 5.5850536) THEN
               JSTAT = 5
            ELSE
               JSTAT = 7
            ENDIF

C  FOUR LINES ARE ATTACHED - LEAN TOWARDS A ROW END NODE

         ELSE

C  CUTOFFS ARE:    0 TO 145 DEGREES = ROW END
C                145 TO 225 DEGREES = ROW SIDE
C                OVER 225 DEGREES   = ROW CORNER (REALLY IMPOSSIBLE)

            IF (ANGLE .LT. 2.5307274) THEN
               JSTAT = 1
            ELSE IF (ANGLE .LT. 3.9269908) THEN
               JSTAT = 3
            ELSE
               JSTAT = 5
            ENDIF
         ENDIF
      ENDIF

      RETURN

      END
