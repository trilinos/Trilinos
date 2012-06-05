C $Id: lcolor.f,v 1.1 1990/11/30 11:10:59 gdsjaar Exp $
C $Log: lcolor.f,v $
C Revision 1.1  1990/11/30 11:10:59  gdsjaar
C Initial revision
C
C
      SUBROUTINE LCOLOR (COLOR)
C***********************************************************************
C
C  SUBROUTINE LCOLOR = SETS THE LINE COLOR
C
C***********************************************************************
C
      CHARACTER*5 COLOR
C
      IF (COLOR .EQ. 'WHITE') THEN
         CALL PLTSTD(1, 7.)
      ELSEIF (COLOR .EQ. 'BLACK') THEN
         CALL PLTSTD(1, 0.)
      ELSEIF (COLOR .EQ. 'YELOW') THEN
         CALL PLTSTD(1, 3.)
      ELSEIF (COLOR .EQ. 'PINK ') THEN
         CALL PLTSTD(1, 5.)
      ELSEIF (COLOR .EQ. 'RED  ') THEN
         CALL PLTSTD(1, 1.)
      ELSEIF (COLOR .EQ. 'BLUE ') THEN
         CALL PLTSTD(1, 4.)
      ENDIF
      RETURN
C
      END
