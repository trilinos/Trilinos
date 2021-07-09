C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE PRINC3(SK1, SK2, SK3, SK4, SK5, SK6, EV, INFO)

      REAL EV(3)
      REAL I1, I2, I3
      INTEGER INFO

      PARAMETER (THIRD = 0.33333 33333 33333 33333 33333)
      SQRT3  = SQRT(3.0)

      INFO = 0

C Find principal trial stresses and directions -

      I1 = ( SK1 + SK2 + SK3 )
      I2 = ( (SK1-SK2)**2 + (SK1-SK3)**2 + (SK2-SK3)**2 ) / 6.0
     *     + SK4**2 + SK5**2 + SK6**2

      S1 = ( (SK1 - SK2) + (SK1 - SK3) ) * THIRD
      S2 = (-(SK1 - SK2) + (SK2 - SK3) ) * THIRD
      S3 = (-(SK1 - SK3) - (SK2 - SK3) ) * THIRD

      I3 = S1 * S2 * S3 + 2.*SK4 * SK5 * SK6
     *     - S1 * SK5**2 - S2 *  SK6**2
     *     - S3 * SK4**2

C     Calculate constants for Malvern Method  (p92)

C .... This sign trick does not work on Ultrix (and others?) with IEEE
C      since IEEE has a signed 0...........
CC      FI2 = I2 + (SIGN (0.5, I2) + SIGN (0.5, -I2))

      IF (I2 .EQ. 0.0) THEN
         FI2 = 1.0
      ELSE
         FI2 = I2
      END IF

      COS3AL = SQRT3 * 1.5 * I3 / FI2 / SQRT(FI2)
      COS3AL = SIGN( MIN( 1.0, ABS(COS3AL) ),COS3AL )

C     ... TRIG FUNCTION USED

      CALPHA = COS( ACOS(COS3AL) / 3.0)
      SALPHA = SQRT(1.0 - CALPHA**2)

      T  = SQRT3 * SQRT(I2)
      P1 = (I1 + T *  2. * CALPHA             ) * THIRD
      P2 = (I1 - T * (CALPHA + SALPHA * SQRT3)) * THIRD
      P3 = (I1 - T * (CALPHA - SALPHA * SQRT3)) * THIRD

C ... Sort Into Correct Position (ev1 < ev2 < ev3)

      EV(1) = MIN(P1, P2, P3)
      EV(3) = MAX(P1, P2, P3)

C ... This is bad from a round-off perspective, but we don't use the
C     value in algebra.  Be warned if use where you need an accurate
C     ev(2).....

      EV(2) = P1 + P2 + P3 - EV(1) - EV(3)

      RETURN
      END
