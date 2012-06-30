C $Id: keep3.f,v 1.1 1990/11/30 11:10:41 gdsjaar Exp $
C $Log: keep3.f,v $
C Revision 1.1  1990/11/30 11:10:41  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]KEEP3.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE KEEP3 (ITEST, LTEST, NBEGIN, NEND)
C***********************************************************************
C
C  SUBROTINE KEEP3 = GETS AN ACCEPTABLE SIDE FOR FILLING TO KEEP A
C                    TRIANGLE VALID
C
C***********************************************************************
C
      DIMENSION ITEST (3), LTEST (3)
C
C  MAKE SURE THAT THE NBEGIN STARTS AT ONE OF THE CORNERS
C
      IF (NBEGIN .EQ. ITEST(1)) THEN
         NEND = ITEST(2)
      ELSEIF (NBEGIN .EQ. ITEST(2)) THEN
         NEND = ITEST(3)
      ELSEIF (NBEGIN .EQ. ITEST(3)) THEN
         NEND = ITEST(1)
      ELSE
         NBEGIN = ITEST(1)
         NEND = ITEST(2)
      ENDIF
C
C  FIND THE CORRECT ROW (THIS ALREADY ASSUMES THAT THE
C  SUM OF THE SMALLER TWO IS EQUAL TO THE LARGEST ONE)
C
      MMAX = MAX0 (LTEST(1), LTEST(2), LTEST(3))
      IF (LTEST(1) .EQ. MMAX) THEN
         IF (NBEGIN .EQ. ITEST(1)) THEN
            NBEGIN = ITEST(2)
            NEND = ITEST(3)
         ENDIF
      ELSEIF (LTEST(2) .EQ. MMAX) THEN
         IF (NBEGIN .EQ. ITEST(2)) THEN
            NBEGIN = ITEST(3)
            NEND = ITEST(1)
         ENDIF
      ELSE
         IF (NBEGIN .EQ. ITEST(3)) THEN
            NBEGIN = ITEST(1)
            NEND = ITEST(2)
         ENDIF
      ENDIF
C
      RETURN
C
      END
