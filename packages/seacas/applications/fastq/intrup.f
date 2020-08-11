C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE INTRUP (PROMPT, IANS, MCOM, ICOM, JCOM, CIN, IIN, RIN,
     &   KIN)
C***********************************************************************

C  SUBROUTINE INTRUP = INPUTS A YES OR NO PLUS MORE IF NEEDED

C***********************************************************************

      DIMENSION IIN (MCOM), RIN (MCOM), KIN (MCOM)
      CHARACTER* (*) PROMPT
      CHARACTER*72 CIN (MCOM), ANS (4)*1, NEWPMT
      LOGICAL IANS
      DATA ANS / 'Y', 'y', 'N', 'n' /

      IZ = 0
      CALL STRLNG (PROMPT, LEN)

C  SEE IF A YES / NO ANSWER IS SITTING AS THE FIRST COMMAND IN THE LIST

      IF ( (ICOM .LE. JCOM) .AND. ( (CIN (ICOM) (1:1) .EQ. ANS (1)) .OR.
     &   (CIN (ICOM) (1:1) .EQ. ANS (2)))) THEN
         IANS = .TRUE.
         ICOM = ICOM + 1
      ELSEIF ( (ICOM .LE. JCOM) .AND. ( (CIN (ICOM) (1:1) .EQ. ANS (3))
     &   .OR. (CIN (ICOM) (1:1) .EQ. ANS (4)))) THEN
         IANS = .FALSE.
         ICOM = ICOM + 1

C  INPUT NEW COMMAND LISTS ONLY IF THE CURRENT ONES ARE USED UP
C  MAKE SURE THE FIRST ONE OF THESE COMMANDS IS EITHER YES OR NO

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

C  OTHERWISE,  JUST GET A YES / NO RESPONSE AND RETURN

      ELSE
         CALL INQTRU (PROMPT, IANS)
      ENDIF
      RETURN

10000 FORMAT (' RESPONSE MUST BE EITHER YES OR NO  -  TRY AGAIN')
      END
