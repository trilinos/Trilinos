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

C $Id: inqtru.f,v 1.1 1990/11/30 11:10:04 gdsjaar Exp $
C $Log: inqtru.f,v $
C Revision 1.1  1990/11/30 11:10:04  gdsjaar
C Initial revision
C
C
      SUBROUTINE INQTRU (PROMPT, IANS)
C***********************************************************************
C
C  SUBROUTINE INQTRU = INPUTS A YES OR NO ANSWER
C
C***********************************************************************
C
      CHARACTER* (*) PROMPT
      CHARACTER*1 RESULT, ANS
      LOGICAL IANS
      DIMENSION ANS (4)
      DATA ANS / 'Y', 'y', 'N', 'n' /
C
  100 CONTINUE
      WRITE (*, 10000)PROMPT
      READ (*, 10010, END = 110, ERR = 120)RESULT
      IF ( (RESULT (1:1) .EQ. ANS (1)) .OR.
     &   (RESULT (1:1) .EQ. ANS (2))) THEN
         IANS = .TRUE.
      ELSEIF ( (RESULT (1:1) .EQ. ANS (3)) .OR.
     &   (RESULT (1:1) .EQ. ANS (4))) THEN
         IANS = .FALSE.
      ELSE
         WRITE (*, 10020)
         GOTO 100
      ENDIF
      RETURN
  110 CONTINUE
      WRITE (*, 10030)
      GOTO 100
  120 CONTINUE
      WRITE (*, 10040)
      GOTO 100
C
10000 FORMAT (' ', A, '? ')
10010 FORMAT (A1)
10020 FORMAT (' RESPONSE MUST BE EITHER YES OR NO  -  TRY AGAIN')
10030 FORMAT (' END OF DATA ENCOUNTERED  -  TRY AGAIN')
10040 FORMAT (' ERROR IN RESPONSE  -  TRY AGAIN')
      END
