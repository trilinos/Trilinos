C $Id: erasec.f,v 1.1 1990/11/30 11:06:45 gdsjaar Exp $
C $Log: erasec.f,v $
C Revision 1.1  1990/11/30 11:06:45  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]ERASEC.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE ERASEC (OLDCUR)
C***********************************************************************
C
C  SUBROUTINE ERASEC = DEACTIVATES THE CROSSHAIRS
C
C***********************************************************************
C
      LOGICAL OLDCUR
C
      IF (OLDCUR) THEN
         WRITE (*,*) CHAR(27)//'G0'
         OLDCUR = .FALSE.
      ENDIF
      RETURN
C
      END
