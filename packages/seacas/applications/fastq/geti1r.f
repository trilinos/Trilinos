C $Id: geti1r.f,v 1.1 1990/11/30 11:08:11 gdsjaar Exp $
C $Log: geti1r.f,v $
C Revision 1.1  1990/11/30 11:08:11  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]GETI1R.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE GETI1R (MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN, I1, R1,
     &   IFOUND)
C***********************************************************************
C
C  SUBROUTINE GETI1R = GETS AN INTEGER AND A REAL INPUT NUMBER
C
C***********************************************************************
C
      DIMENSION IIN(MCOM), RIN(MCOM), KIN(MCOM)
      CHARACTER*72 CIN(MCOM)
C
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
C
      END
