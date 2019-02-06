C    Copyright(C) 2009-2017 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C
C        * Redistributions of source code must retain the above copyright
C          notice, this list of conditions and the following disclaimer.
C
C        * Redistributions in binary form must reproduce the above
C          copyright notice, this list of conditions and the following
C          disclaimer in the documentation and/or other materials provided
C          with the distribution.
C
C        * Neither the name of NTESS nor the names of its
C          contributors may be used to endorse or promote products derived
C          from this software without specific prior written permission.
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

      CHARACTER*16 FUNCTION ENGNOT (RNUM, IPREC)
      REAL RNUM

      CHARACTER*10 TFORM
      CHARACTER*32 SCRSTR

      NSIG = IPREC
      NSIG = MIN(8, MAX(2, NSIG)) - 1

         ENGNOT = ' '
         DO 10 I=1,3
            WRITE (TFORM,  30) I, NSIG+I+7, NSIG+I
            WRITE (SCRSTR, TFORM) RNUM
            READ  (SCRSTR(NSIG+I+6:NSIG+I+7), 40) IEXP
            IF (MOD(IEXP, 3) .EQ. 0) THEN
               ENGNOT(9-(NSIG+I):16) = SCRSTR(:LENSTR(SCRSTR))
               RETURN
            END IF
   10    CONTINUE
      RETURN
   30 FORMAT ('(',I1,'PE',I2.2,'.',I2.2,')')
   40 FORMAT (I2)
      END
