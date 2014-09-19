C    Copyright (c) 2014, Sandia Corporation.
C    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C    the U.S. Governement retains certain rights in this software.
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
C        * Neither the name of Sandia Corporation nor the names of its
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
C    

C $Id: range.f,v 1.2 1999/02/16 21:38:01 gdsjaar Exp $
C $Log: range.f,v $
C Revision 1.2  1999/02/16 21:38:01  gdsjaar
C Converted to read exodusII database format.  Somewhat tested, not
C ready for production yet.
C
C Revision 1.1.1.1  1991/02/21 15:45:22  gdsjaar
C NUMBERS: Greg Sjaardema, initial Unix release
C
c Revision 1.1  1991/02/21  15:45:21  gdsjaar
c Initial revision
c
      SUBROUTINE RANGE (LEN, LIST, IOMIN, IOMAX)
      LOGICAL LIST(*), INRNG
      INTEGER IRANGE(3)
      CHARACTER*80 LINE

      INRNG = .FALSE.
      IRINC = 0

      DO 10 I=1, LEN
         IF (LIST(I)) THEN
            IBEG = I+1
            LASTSL = I
            GO TO 20
         END IF
   10 CONTINUE

   20 CONTINUE
      DO 40 I = IBEG, LEN
         IF (LIST(I)) THEN
            IRINC1 = I - LASTSL
            IF (IRINC .EQ. 0) IRINC = IRINC1
            IF (.NOT. INRNG) THEN
               IRBEG = LASTSL
               IRINC = IRINC1
               INRNG = .TRUE.
            ELSE IF (INRNG .AND. IRINC1 .NE. IRINC) THEN
               IREND = I
               INRNG = .FALSE.
               IRANGE(1) = IRBEG
               IRANGE(2) = IREND - IRINC1
               IRANGE(3) = IRINC
               LINE = ' '
               CALL FFADDV (IRANGE, LINE)
               DO 30 IO = IOMIN, IOMAX
                  CALL LOGERR ('CMDSPEC', LINE(:LENSTR(LINE)), IO)
   30          CONTINUE
               IRINC = IRINC1
            END IF
            LASTSL = I
         END IF
   40 CONTINUE
      IRANGE(1) = IRBEG
      IRANGE(2) = LASTSL
      IRANGE(3) = IRINC
      LINE = ' '
      CALL FFADDV (IRANGE, LINE)
      DO 50 IO = IOMIN, IOMAX
         CALL LOGERR ('CMDSPEC', LINE(:LENSTR(LINE)), IO)
   50 CONTINUE
      RETURN
      END
