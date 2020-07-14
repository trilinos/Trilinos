C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C    
C    See packages/seacas/LICENSE for details

C $Id: summry.f,v 1.2 2000/07/06 16:49:57 gdsjaar Exp $
C $Log: summry.f,v $
C Revision 1.2  2000/07/06 16:49:57  gdsjaar
C Changed real*4 to real
C
C Revision 1.1.1.1  1991/02/21 15:45:55  gdsjaar
C NUMBERS: Greg Sjaardema, initial Unix release
C
c Revision 1.1  1991/02/21  15:45:54  gdsjaar
c Initial revision
c
      SUBROUTINE SUMMRY (TYPE, NUM, SELECT, VALUE, SUMR, ISUMR, IOFF)
      CHARACTER*(*) TYPE
      LOGICAL SELECT(*)
      REAL VALUE(*), SUMR(*)
      INTEGER ISUMR(*)
C
C ... TYPE = 'A' - calculate stats based on absolute values
C              NOTE: VALUE will be modified if this option used
C          = ' ' - calculate stats based on true value
C
C ... SUMR(1) = MINIMUM   ISUMR(1) = ELEMENT NUMBER
C ... SUMR(2) = MAXIMUM   ISUMR(2) = ELEMENT NUMBER
C ... SUMR(3) = AVERAGE
C ... SUMR(4) = STD. DEV.
C
      SUMR(1) =  1.0E30
      SUMR(2) = -1.0E30
      NUMSEL = 0

      RMEAN  = 0.0
      STDDEV = 0.0

      IF (TYPE(:1) .EQ. 'A') THEN
         DO 5 I=1, NUM
            VALUE(I) = ABS(VALUE(I))
    5    CONTINUE
      END IF

      DO 10 I = 1, NUM
         IF (SELECT(I)) THEN
            NUMSEL = NUMSEL + 1
            TMEAN  = RMEAN + (VALUE(I) - RMEAN) / NUMSEL
            STDDEV = STDDEV + (VALUE(I) - RMEAN) * (VALUE(I) - TMEAN)
            RMEAN  = TMEAN

            IF (VALUE(I) .LT. SUMR(1)) THEN
               SUMR(1) = VALUE(I)
               ISUMR(1) = I
            END IF

            IF (VALUE(I) .GT. SUMR(2)) THEN
               SUMR(2) = VALUE(I)
               ISUMR(2) = I
            END IF

         END IF
   10 CONTINUE

      SUMR(3) = RMEAN
      SUMR(4) = SQRT(STDDEV / MAX(1.0, DBLE(NUMSEL-1)))

      ISUMR(1) = ISUMR(1) + IOFF
      ISUMR(2) = ISUMR(2) + IOFF

      RETURN
      END
