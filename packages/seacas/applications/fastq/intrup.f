C $Id: intrup.f,v 1.1 1990/11/30 11:10:17 gdsjaar Exp $
C $Log: intrup.f,v $
C Revision 1.1  1990/11/30 11:10:17  gdsjaar
C Initial revision
C
C
      SUBROUTINE INTRUP (PROMPT, IANS, MCOM, ICOM, JCOM, CIN, IIN, RIN,
     &   KIN)
C***********************************************************************
C
C  SUBROUTINE INTRUP = INPUTS A YES OR NO PLUS MORE IF NEEDED
C
C***********************************************************************
C
      DIMENSION IIN (MCOM), RIN (MCOM), KIN (MCOM)
      CHARACTER* (*) PROMPT
      CHARACTER*72 CIN (MCOM), ANS (4)*1, NEWPMT
      LOGICAL IANS
      DATA ANS / 'Y', 'y', 'N', 'n' /
C
      IZ = 0
      CALL STRLNG (PROMPT, LEN)
C
C  SEE IF A YES / NO ANSWER IS SITTING AS THE FIRST COMMAND IN THE LIST
C
      IF ( (ICOM .LE. JCOM) .AND. ( (CIN (ICOM) (1:1) .EQ. ANS (1)) .OR.
     &   (CIN (ICOM) (1:1) .EQ. ANS (2)))) THEN
         IANS = .TRUE.
         ICOM = ICOM + 1
      ELSEIF ( (ICOM .LE. JCOM) .AND. ( (CIN (ICOM) (1:1) .EQ. ANS (3))
     &   .OR. (CIN (ICOM) (1:1) .EQ. ANS (4)))) THEN
         IANS = .FALSE.
         ICOM = ICOM + 1
C
C  INPUT NEW COMMAND LISTS ONLY IF THE CURRENT ONES ARE USED UP
C  MAKE SURE THE FIRST ONE OF THESE COMMANDS IS EITHER YES OR NO
C
      ELSEIF (ICOM .GT. JCOM) THEN
         IF (LEN .LE. 71) THEN
            NEWPMT = PROMPT (1:LEN)
            NEWPMT (LEN + 1:LEN + 1) = '?'
         ELSE
            NEWPMT = PROMPT
         ENDIF
         CALL STRLNG (NEWPMT, NEWLEN)
         NEWLEN = MIN0 (72, NEWLEN + 1)
  100    CONTINUE
         CALL MESAGE (' ')
         CALL FREFLD (IZ, IZ, NEWPMT (1:NEWLEN), MCOM, IOSTAT, JCOM,
     &      KIN, CIN, IIN, RIN)
         ICOM = 1
         IF ( (CIN (ICOM) (1:1) .EQ. ANS (1)) .OR.
     &      (CIN (ICOM) (1:1) .EQ. ANS (2))) THEN
            IANS = .TRUE.
            ICOM = ICOM + 1
         ELSEIF ( (CIN (ICOM) (1:1) .EQ. ANS (3)) .OR.
     &      (CIN (ICOM) (1:1) .EQ. ANS (4))) THEN
            IANS = .FALSE.
            ICOM = ICOM + 1
         ELSE
            WRITE (*, 10000)
            GOTO 100
         ENDIF
C
C  OTHERWISE,  JUST GET A YES / NO RESPONSE AND RETURN
C
      ELSE
         CALL INQTRU (PROMPT, IANS)
      ENDIF
      RETURN
C
10000 FORMAT (' RESPONSE MUST BE EITHER YES OR NO  -  TRY AGAIN')
      END
