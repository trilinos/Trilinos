C $Id: geti12.f,v 1.1 1990/11/30 11:08:08 gdsjaar Exp $
C $Log: geti12.f,v $
C Revision 1.1  1990/11/30 11:08:08  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]GETI12.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2,
     &   IFOUND)
C***********************************************************************
C
C  SUBROUTINE GETI12 = GETS TWO INTEGERS
C
C***********************************************************************
C
      DIMENSION IIN (MCOM), KIN (MCOM)
      CHARACTER*72 CIN (MCOM)
C
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
C
      END
