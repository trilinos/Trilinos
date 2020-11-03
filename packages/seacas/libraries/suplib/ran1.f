C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      REAL FUNCTION RAN1(idum)

C  This function returns a pseudo-random number for each invocation.
C  It is a FORTRAN 77 adaptation of the "Integer Version 2" minimal
C  standard number generator whose Pascal code appears in the article:

C     Park, Steven K. and Miller, Keith W., "Random Number Generators:
C     Good Ones are Hard to Find", Communications of the ACM,
C     October, 1988.

      PARAMETER (MPLIER=16807,MODLUS=2147483647,MOBYMP=127773,
     +           MOMDMP=2836)

      data jseed /123456789/
      data ifrst /0/

      INTEGER HVLUE, LVLUE, TESTV, NEXTN
      SAVE    NEXTN

      IF (IFRST .EQ. 0) THEN
        if (idum .ne. 0) then
          nextn = idum
        else
          NEXTN = JSEED
        end if

        IFRST = 1
      ENDIF

      HVLUE = NEXTN / MOBYMP
      LVLUE = MOD(NEXTN, MOBYMP)
      TESTV = MPLIER*LVLUE - MOMDMP*HVLUE
      IF (TESTV .GT. 0) THEN
        NEXTN = TESTV
      ELSE
        NEXTN = TESTV + MODLUS
      ENDIF
      RAN1 = REAL(NEXTN)/REAL(MODLUS)

      RETURN
      END
