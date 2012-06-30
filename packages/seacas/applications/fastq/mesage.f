C $Id: mesage.f,v 1.1 1990/11/30 11:11:59 gdsjaar Exp $
C $Log: mesage.f,v $
C Revision 1.1  1990/11/30 11:11:59  gdsjaar
C Initial revision
C
C
      SUBROUTINE MESAGE (PROMPT)
C***********************************************************************
C
C  SUBROUTINE MESAGE = PRINTS A MESSAGE ONTO THE SCREEN
C
C***********************************************************************
C
      CHARACTER * (*) PROMPT
C
      IF (PROMPT .EQ. ' ') THEN
         WRITE (*, 10000)
      ELSE
         WRITE (*, 10010)PROMPT
      ENDIF
      RETURN
C
10000 FORMAT ( / )
10010 FORMAT (' ', A)
      END
