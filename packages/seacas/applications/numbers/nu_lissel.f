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

C $Id: lissel.f,v 1.1 1991/02/21 15:43:56 gdsjaar Exp $
C $Log: lissel.f,v $
C Revision 1.1  1991/02/21 15:43:56  gdsjaar
C Initial revision
C
C=======================================================================
      SUBROUTINE LISSEL (OPT, TYPE, IOMIN, IOMAX, LIST, SELECT, NUMLST)
C=======================================================================
C
C ... Output selected entities to Terminal and/or List file
C
C     OPT = IN = Option:
C                 'L' = Selected by logical list
C                 'A' = All selected are in list
C                 'R' = List in Range form
C
      LOGICAL SELECT(*)
      CHARACTER*(*) OPT, TYPE
      INTEGER LIST(*), ISCR(12)
      CHARACTER*80 STRING
      LOGICAL ISABRT

      WRITE (STRING, 10) TYPE(:LENSTR(TYPE))
   10 FORMAT ('List of Selected ',A)
      CALL SQZSTR(STRING,LSTR)
      DO 20 IO=IOMIN, IOMAX
         WRITE (IO, 30) STRING(:LSTR)
   20 CONTINUE
   30 FORMAT (/,1X,A/)

      II = 0
      IF (INDEX(OPT,'R') .GT. 0) THEN
         CALL RANGE(NUMLST, SELECT, IOMIN, IOMAX)
      ELSE
      DO 50 I = 1, NUMLST
         IF (SELECT(I)) THEN
            II = II + 1
            ISCR(II) = I
            IF (II .EQ. 12) THEN
               IF (ISABRT()) RETURN
               DO 40 IO = IOMIN, IOMAX
                  WRITE (IO, 70) (ISCR(ILST),ILST=1,II)
   40          CONTINUE
               II = 0
            END IF
         END IF
   50 CONTINUE
      IF (II .GT. 0) THEN
         DO 60 IO = IOMIN, IOMAX
            WRITE (IO, 70) (ISCR(ILST),ILST=1,II)
   60    CONTINUE
      END IF
   70 FORMAT ((1X,12I6))
      END IF
      RETURN
      END
