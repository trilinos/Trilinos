C $Id: getr12.f,v 1.1 1990/11/30 11:08:37 gdsjaar Exp $
C $Log: getr12.f,v $
C Revision 1.1  1990/11/30 11:08:37  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]GETR12.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE GETR12 (MCOM, ICOM, JCOM, CIN, RIN, KIN, R1, R2,
     &   IFOUND)
C***********************************************************************
C
C  SUBROUTINE GETR12 = GETS TWO REAL NUMBERS
C
C***********************************************************************
C
      DIMENSION RIN(MCOM), KIN(MCOM)
      CHARACTER*72 CIN(MCOM)
C
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
C
      END
