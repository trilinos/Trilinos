C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTINO(MIN,MAX,START,REND,INTER,EXP,NMIN)
      REAL MIN,MAX,INTER,NINTER
      INTEGER EXP
      REAL NI
      DATA SMALL/1.E-4/,DESINT/5./

      DELTA = MAX - MIN
      TMAX = MAX
      TMIN = MIN
      IF (DELTA.NE.0) THEN
         EXP = NINT(LOG10(ABS(DELTA))) - 1

      ELSE
         IF (MIN.NE.0) THEN
            EXP = NINT(LOG10(ABS(MIN))) - 1
            EPS = .01*ABS(MIN)

         ELSE
            EXP = 0
            EPS = .1
         END IF

         TMAX = MAX + EPS
         TMIN = TMIN - EPS
      END IF

      TENEXP = 10.**EXP
      RMIN = 1.E20
      J = 1
      DO 2470 I = 1,5
         NINTER = DELTA/ (DBLE(I)*TENEXP)
         TEMP = ABS(DESINT-NINTER)
         IF (TEMP.LT.RMIN) THEN
            J = I
            RMIN = TEMP
         END IF

 2470 CONTINUE
      INTER = DBLE(J)
      IF (DELTA.EQ.0.) THEN
         INTER = 1.
      END IF

      IF (INTER.EQ.1.) THEN
         NMIN = 5

      ELSE IF (INTER.EQ.2.) THEN
         NMIN = 4

      ELSE IF (INTER.EQ.3.) THEN
         NMIN = 6

      ELSE IF (INTER.EQ.4.) THEN
         NMIN = 4

      ELSE IF (INTER.EQ.5.) THEN
         NMIN = 5
      END IF

      TENI = INTER*TENEXP
      IF (TMIN.GE.0.) THEN
         ADD = 0.

      ELSE
         ADD = -1.
      END IF

      RNI = TMIN/TENI
      NI = INT(RNI)
      IF (ABS(RNI-NI).LT.SMALL) THEN
         ADD = 0.
      END IF

      START = (NI+ADD)*TENI
      IF (TMAX.GE.0.) THEN
         ADD = 1.

      ELSE
         ADD = 0.
      END IF

      RNI = TMAX/TENI
      NI = INT(RNI)
      IF (ABS(RNI-NI).LT.SMALL) THEN
         ADD = 0.
      END IF

      REND = (NI+ADD)*TENI
      START = START/TENEXP
      REND = REND/TENEXP
      IF (REND.NE.0.) THEN
 2490    IF (.NOT. (ABS(REND).GT.10.)) GO TO 2500
         REND = REND/10.
         START = START/10.
         INTER = INTER/10.
         EXP = EXP + 1
         GO TO 2490

 2500    CONTINUE
      END IF

      IF (START.NE.0.) THEN
 2510    IF (.NOT. (ABS(START).LT.1.)) GO TO 2520
         REND = REND*10.
         START = START*10.
         INTER = INTER*10.
         EXP = EXP - 1
         GO TO 2510

 2520    CONTINUE
      END IF

      IF (START.EQ.0 .OR. TMIN.EQ.0) THEN
         RETURN

      END IF

      IF (START.EQ.REND) THEN
         REND = START + INTER
      END IF

      IF (ABS(START-REND).EQ.INTER) THEN
         NMIN = 10
      END IF

      RETURN

      END
