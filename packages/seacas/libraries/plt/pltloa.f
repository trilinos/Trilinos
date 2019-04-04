C Copyright (C) 2009-2017 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.
C
C     * Neither the name of NTESS nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C

C $Id: pltloa.f,v 1.1 1993/07/16 16:48:43 gdsjaar Exp $
C $Log: pltloa.f,v $
C Revision 1.1  1993/07/16 16:48:43  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      SUBROUTINE PLTLOA(LINE1,NUM,TYPE)
      CHARACTER*10 LINE
      CHARACTER*(*) LINE1
      INTEGER TYPE

      LINE1 = ' '
      LINE = ' '
      IF (TYPE.EQ.1) THEN
         IF (NUM.GE.0) THEN
            WRITE (LINE,10,ERR=20) INT(10.**NUM)

   10       FORMAT (I10)

   20       DO 2680 J = 1,10
               IF (LINE(J:J).NE.' ') THEN
                  GO TO 2690

               END IF

 2680       CONTINUE
 2690       CONTINUE
            LINE1 = LINE(J:)
            RETURN

         END IF

         LINE1(1:1) = '.'
         DO 2700 I = 1,ABS(NUM) - 1
            LINE1(I+1:I+1) = '0'
 2700    CONTINUE
         LINE1(I+1:I+1) = '1'

      ELSE
         WRITE (LINE,10,ERR=40) NUM
   40    DO 2720 J = 1,10
            IF (LINE(J:J).NE.' ') THEN
               GO TO 2730

            END IF

 2720    CONTINUE
 2730    CONTINUE
         LINE1 = LINE(J:)
      END IF

      RETURN

      END
