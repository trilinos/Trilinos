C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE INRENM (MSC, N23, CFLAG, RIN, IIN, IFOUND, NUMBER,
     &   NOROOM)
C***********************************************************************

C  SUBROUTINE INRENM = INPUTS A RENUMBERING CARD

C***********************************************************************

      DIMENSION NUMBER (MSC), RIN (IFOUND), IIN (IFOUND)

      CHARACTER * 80 NUMBER, CFLAG * 72

      LOGICAL NOROOM

      NOROOM = .TRUE.

      N23 = N23 + 1
      IF (N23 .GT. MSC)RETURN
      NUMBER (N23) = ' '
      NUMBER (N23) (1:5) = CFLAG (1:5)

C  INPUT A POINT - LINE - POINT CARD

      IF (CFLAG (1:5) .EQ. 'P-L-P') THEN
         IFOUND = MIN0 (IFOUND, 15)
         DO 100 IJ = 1, IFOUND
            I2 =  (IJ + 1) * 5
            I1 = I2 - 4
            WRITE (NUMBER (N23) (I1:I2), 10000)IIN (IJ)
  100    CONTINUE

C  INPUT AN X, Y LOCATION RENUMBERING CARD

      ELSEIF (CFLAG (1:3) .EQ. 'X-Y') THEN
         WRITE (NUMBER (N23) (11:20), 10010) RIN (1)
         WRITE (NUMBER (N23) (21:30), 10010) RIN (2)

C  INPUT A NODE UNIQUE ID RENUMBERING CARD

      ELSEIF (CFLAG (1:4) .EQ. 'NODE') THEN
         IFOUND = MIN0 (IFOUND, 7)
         DO 110 IJ = 1, IFOUND
            I2 =  ( (IJ + 1) * 10)
            I1 = I2 - 9
            WRITE (NUMBER (N23) (I1:I2), 10020)IIN (IJ)
  110    CONTINUE

C  INDICATE ERROR IN RENUMBERING FLAG

      ELSE
         N23 = N23 - 1
         WRITE ( * , 10030) CFLAG (1:5)
      ENDIF

      NOROOM = .FALSE.
      RETURN

10000 FORMAT (I5)
10010 FORMAT (1PE10.3)
10020 FORMAT (I10)
10030 FORMAT (' RENUMBERING KEY WORD: ', A5, ' IS NOT ALLOWABLE',  / ,
     &   ' THIS RENUMBERING LIST WILL NOT BE INPUT INTO DATABASE')
      END
