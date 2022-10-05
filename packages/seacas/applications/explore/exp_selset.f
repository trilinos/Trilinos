C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE SELSET (NUMSEL, IXSEL, NUMSET, LISNPS, LNPSNL,
     &  IDNPS, NNNPS, IXNNPS, LTNNPS, VALNAM)
C=======================================================================

      INTEGER IXSEL(*)
      INTEGER LISNPS(0:*)
      INTEGER IDNPS(*)
      INTEGER NNNPS(*)
      INTEGER IXNNPS(*)
      INTEGER LTNNPS(*)
      CHARACTER*(*) VALNAM
      CHARACTER*40 STRA
      CHARACTER*132 MSG

      NUMSEL = 0
      DO 100 IX = 1, LISNPS(0)
        INPS = LISNPS(IX)
        IS = IXNNPS(INPS)
        IE = IS + NNNPS(INPS) - 1
        do 90 i=is,ie
          NUMSEL = NUMSEL + 1
          IXSEL(NUMSEL) = ltnnps(i)
 90     continue
 100  CONTINUE
      write (stra, 10000) numsel
      call pckstr(1, stra)
      MSG = STRA(:lenstr(stra)) // ' ' // VALNAM // ' selected'
      call prterr('CMDSPEC', MSG)
10000 format(I12)
      RETURN
      END

