C $Id: initdg.f,v 1.1 1990/11/30 11:09:49 gdsjaar Exp $
C $Log: initdg.f,v $
C Revision 1.1  1990/11/30 11:09:49  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]INITDG.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE INITDG (MCOM, ICOM, JCOM, CIN, RIN, IIN, KIN, IDUMP,
     &   XX1, YY1, SCALE, CT, ST, X1, X2, Y1, Y2, DRWTAB, SNAP)
C***********************************************************************
C
C  SUBROUTINE INITDG = INITIALIZES THE DIGITIZING TABLET
C
C***********************************************************************
C
      DIMENSION KIN (MCOM), IIN (MCOM), RIN (MCOM)
C
      CHARACTER * 72 CIN (MCOM), BUTTON * 1
C
      LOGICAL DRWTAB, IANS, SNAP
C
      IZ = 0
C
C  CHECK TO MAKE SURE THAT THE DRAWING IS NOT BEING TOGGLED
C
      IF (DRWTAB) THEN
         CALL MESAGE ('DRAWING INITIALIZATION IS ALREADY ACTIVE')
         CALL INTRUP  ('TOGGLE ALL DRAWING INITIALIZATION OFF',
     &      IANS,  MCOM,  ICOM,  JCOM,  CIN,  IIN,  RIN,  KIN)
         IF (IANS) THEN
            DRWTAB = .FALSE.
            CALL TABINT (X1, X2, Y1, Y2, CT, ST, SCALE, XX1, YY1, XX2,
     &         YY2, DRWTAB)
            RETURN
         ENDIF
      ENDIF
C
C
C  GET THE ZOOM LIMITS
C
      CALL MESAGE (' ')
      IF (ICOM .GT. JCOM) THEN
         CALL FREFLD (IZ, IZ, 'ENTER DRAWING XMIN, XMAX, YMIN, YMAX:',
     &      MCOM, IOSTAT, JCOM, KIN, CIN, IIN, RIN)
         ICOM = 1
      ENDIF
      IF ( (JCOM - ICOM + 1) .GE. 4) THEN
         SNAP = .TRUE.
         X1 = RIN (ICOM)
         X2 = RIN (ICOM + 1)
         Y1 = RIN (ICOM + 2)
         Y2 = RIN (ICOM + 3)
         ICOM = ICOM + 4
      ELSE
         CALL MESAGE ('NOT ENOUGH INFORMATION DEFINED TO SPECIFY'//
     &      ' DRAWING LIMITS')
         CALL MESAGE ('INITIALIZATION ABORTED')
         CALL MESAGE (' ')
         RETURN
      ENDIF
C
C  GET THE DIGITIZING POINTS
C
      CALL MESAGE ('NOW DIGITIZE THOSE 2 POINTS')
      CALL MESAGE ('       PUSH "PUCK - 1" FOR LOWER LEFT')
      CALL MESAGE ('       PUSH "PUCK - 2" FOR UPPER RIGHT')
      CALL MESAGE ('       PUSH "PUCK - E" TO END')
  100 CONTINUE
      CALL DPREAD (X, Y, BUTTON)
      IF (BUTTON .EQ. '1') THEN
         XX1 = X
         YY1 = Y
         CALL MESAGE ('LOWER LEFT INPUT')
         GOTO 100
      ELSEIF (BUTTON .EQ. '2') THEN
         XX2 = X
         YY2 = Y
         CALL MESAGE ('UPPER RIGHT INPUT')
         GOTO 100
      ELSEIF (BUTTON .EQ. 'E') THEN
         CALL PLTBEL
         CALL PLTFLU
      ENDIF
      IF ( ( (YY2 - YY1 .EQ. 0.) .AND. (XX2 - XX1 .EQ. 0.))
     &   .OR. ( (Y2 - Y1 .EQ. 0.) .AND. (X2 - X1 .EQ. 0.))) THEN
         CALL MESAGE ('BAD INITIALIZATION  -  INITIALIZATION ABORTED')
         CALL MESAGE (' ')
         CALL PLTBEL
         CALL PLTFLU
         RETURN
      ENDIF
      DRWTAB = .TRUE.
      CALL TABINT (X1, X2, Y1, Y2, CT, ST, SCALE, XX1, YY1, XX2, YY2,
     &   DRWTAB)
      CALL MESAGE ('INITIALIZATION COMPLETE')
      CALL MESAGE (' ')
C
      RETURN
C
      END
