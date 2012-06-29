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

C $Id: chrstr.f,v 1.1 1993/07/16 16:46:22 gdsjaar Exp $ 
C $Log: chrstr.f,v $
C Revision 1.1  1993/07/16 16:46:22  gdsjaar
C Changed plt to library rather than single source file.
C 
C=======================================================================
      SUBROUTINE CHRSTR(LINE1,LINE2,L)
      CHARACTER*(*) LINE1,LINE2
      CHARACTER CH

      K = LEN(LINE1)
      J = 1
 2150 IF (.NOT. (J.LT.K)) GO TO 2170
      CH = LINE1(J:J)
      IF (CH.EQ.' ') THEN
         GO TO 2160

      END IF

      IF (CH.EQ.CHAR(9)) THEN
         GO TO 2160

      END IF

      GO TO 2170

 2160 J = J + 1
      GO TO 2150

 2170 CONTINUE
      IF (J.GT.K) THEN
         L = 0
         LINE2 = ' '
         RETURN

      END IF

      LINE2(1:K-J+1) = LINE1(J:K)
      CALL CHRTRM(LINE2(1:K-J+1),L)
      RETURN

      END
