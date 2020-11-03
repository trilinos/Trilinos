C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2,
     &   IFOUND)
C***********************************************************************

C  SUBROUTINE GETI12 = GETS TWO INTEGERS

C***********************************************************************

      DIMENSION IIN (MCOM), KIN (MCOM)
      CHARACTER*72 CIN (MCOM)

      IF ((ICOM .GT. JCOM) .OR. (CIN (ICOM) (1:1) .EQ. ' ')) THEN
         IFOUND = 0
         ICOM = ICOM + 1
         I1 = 0
         I2 = 0
      ELSE
         I1 = IIN (ICOM)
         ICOM = ICOM + 1
         IF ((ICOM .LE. JCOM) .AND. (KIN (ICOM) .GT. 0)) THEN
            I2 = IIN (ICOM)
            ICOM = ICOM + 1
            IFOUND = 2
         ELSE
            I2 = 0
            IFOUND = 1
         ENDIF
      ENDIF
      RETURN

      END
