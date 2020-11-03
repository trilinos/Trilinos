C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE GETI1R (MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN, I1, R1,
     &   IFOUND)
C***********************************************************************

C  SUBROUTINE GETI1R = GETS AN INTEGER AND A REAL INPUT NUMBER

C***********************************************************************

      DIMENSION IIN(MCOM), RIN(MCOM), KIN(MCOM)
      CHARACTER*72 CIN(MCOM)

      IF( (ICOM .GT. JCOM) .OR. (CIN(ICOM) (1:1) .EQ. ' ') ) THEN
         ICOM = ICOM+1
         IFOUND = 0
      ELSE
         I1 = IIN(ICOM)
         ICOM = ICOM+1
         IF ( (ICOM .LE. JCOM) .AND. (KIN(ICOM) .GT. 0) ) THEN
            R1 = RIN(ICOM)
            ICOM = ICOM + 1
            IFOUND = 2
         ELSE
            R1 = 0.
            IFOUND = 1
         ENDIF
      ENDIF
      RETURN

      END
