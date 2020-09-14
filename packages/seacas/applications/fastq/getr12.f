C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE GETR12 (MCOM, ICOM, JCOM, CIN, RIN, KIN, R1, R2,
     &   IFOUND)
C***********************************************************************

C  SUBROUTINE GETR12 = GETS TWO REAL NUMBERS

C***********************************************************************

      DIMENSION RIN(MCOM), KIN(MCOM)
      CHARACTER*72 CIN(MCOM)

      IF ( (ICOM .GT. JCOM) .OR. (CIN(ICOM) (1:1) .EQ. ' ') ) THEN
         IFOUND = 0
         ICOM = ICOM+1
         R1 = 0.
         R2 = 0.
      ELSE
         R1 = RIN(ICOM)
         ICOM = ICOM+1
         IF ( (ICOM .LE. JCOM) .AND. (KIN(ICOM) .GT. 0) ) THEN
            R2 = RIN(ICOM)
            ICOM = ICOM+1
            IFOUND = 2
         ELSE
            R2 = 0.
            IFOUND = 1
         ENDIF
      ENDIF
      RETURN

      END
