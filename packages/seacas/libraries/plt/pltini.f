C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTINI(MIN,MAX,START,REND,INTER,EXP,NMIN)
      REAL MIN,MAX,INTER
      INTEGER EXP

      DELTA = MAX - MIN
      IF (DELTA.LT.0.) THEN
         CALL PLTFLU
         CALL SIORPT('PLTINI',
     *               'Maximum value must be greater than minimum value.'
     *               ,2)
         RETURN

      END IF

      CALL PLTINO(MIN,MAX,START,REND,INTER,EXP,NMIN)
      TENEXP = 10.**EXP
      IEXP = NINT(LOG10(ABS(INTER))) - 2
      SMALL = 10.**IEXP
      IF (ABS(START-MIN/TENEXP).GT.SMALL) THEN
         START = START + INTER
      END IF

      IF (ABS(REND-MAX/TENEXP).GT.SMALL) THEN
         REND = REND - INTER
      END IF

      RETURN

      END
