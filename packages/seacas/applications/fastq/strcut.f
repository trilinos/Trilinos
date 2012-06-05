C $Id: strcut.f,v 1.1 1990/11/30 11:16:42 gdsjaar Exp $
C $Log: strcut.f,v $
C Revision 1.1  1990/11/30 11:16:42  gdsjaar
C Initial revision
C
C
      SUBROUTINE STRCUT (STRING)
C***********************************************************************
C
C  SUBROUTINE STRCUT = DELETES ALL PRECEDING BLANKS FROM A STRING
C
C***********************************************************************
C
      CHARACTER * (*) STRING, HOLD*80
      CALL STRIPB (STRING, ILEFT, IRIGHT)
      IF (IRIGHT .GT. ILEFT) THEN
         HOLD = STRING (ILEFT:)
         STRING = HOLD
      ELSE
         HOLD = ' '
      ENDIF
      RETURN
      END
