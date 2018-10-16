C    Copyright(C) 2014-2017 National Technology & Engineering Solutions of
C    Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C    * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C    
C    * Redistributions in binary form must reproduce the above
C      copyright notice, this list of conditions and the following
C      disclaimer in the documentation and/or other materials provided
C      with the distribution.
C    
C    * Neither the name of NTESS nor the names of its
C      contributors may be used to endorse or promote products derived
C      from this software without specific prior written permission.
C    
C    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C    

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
