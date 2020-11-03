C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE KEEP34 (ITEST, LTEST, NBEGIN, NEND, ICHNG)
C***********************************************************************

C  SUBROTINE KEEP34 = GETS AN ACCEPTABLE SIDE FOR FILLING TO KEEP A
C                     TRIANGLE VALID OR CHANGING TO A RECTANGLE

C***********************************************************************

      DIMENSION ITEST (3), LTEST (3)

C  MAKE SURE THAT THE NBEGIN STARTS AT ONE OF THE CORNERS

      IF (NBEGIN .EQ. ITEST(1)) THEN
         NEND = ITEST(2)
      ELSEIF (NBEGIN .EQ. ITEST(2)) THEN
         NEND = ITEST(3)
      ELSEIF (NBEGIN .EQ. ITEST(3)) THEN
         NEND = ITEST(1)
      ELSE
         NBEGIN = ITEST(1)
         NEND = ITEST(3)
      ENDIF

C  FIND THE CORRECT ROW (THIS ALREADY ASSUMES THAT THE
C  SUM OF THE SMALLER TWO IS EQUAL TO THE LARGEST ONE)

      MMAX = MAX0 (LTEST(1), LTEST(2), LTEST(3))
      IF (LTEST(1) .EQ. MMAX) THEN
         IF (NBEGIN .EQ. ITEST(1)) THEN
            NEND = ICHNG
         ENDIF
      ELSEIF (LTEST(2) .EQ. MMAX) THEN
         IF (NBEGIN .EQ. ITEST(2)) THEN
            NEND = ICHNG
         ENDIF
      ELSE
         IF (NBEGIN .EQ. ITEST(3)) THEN
            NEND = ICHNG
         ENDIF
      ENDIF

      RETURN

      END
