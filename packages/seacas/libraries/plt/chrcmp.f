C Copyright (C) 2009 Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software
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
C     * Neither the name of Sandia Corporation nor the names of its
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

C $Id: chrcmp.f,v 1.1 1993/07/16 16:46:17 gdsjaar Exp $ 
C $Log: chrcmp.f,v $
C Revision 1.1  1993/07/16 16:46:17  gdsjaar
C Changed plt to library rather than single source file.
C 
C=======================================================================
      LOGICAL FUNCTION CHRCMP(KWD,PART1,PART2)
      CHARACTER*(*) KWD
      CHARACTER*(*) PART1
      CHARACTER*(*) PART2

      CALL CHRTRM(KWD,LK)
      CALL CHRTRM(PART1,LF)
      CALL CHRTRM(PART2,LV)
      IF (LK.LT.LF .OR. LK.GT.LF+LV) THEN
         CHRCMP = .FALSE.
         RETURN

      END IF

      IF (KWD(1:LF).NE.PART1(1:LF)) THEN
         CHRCMP = .FALSE.
         RETURN

      END IF

      IF (LK.EQ.LF) THEN
         CHRCMP = .TRUE.
         RETURN

      END IF

      IF (KWD(LF+1:LK).NE.PART2(1:LK-LF)) THEN
         CHRCMP = .FALSE.
         RETURN

      END IF

      CHRCMP = .TRUE.
      RETURN

      END
