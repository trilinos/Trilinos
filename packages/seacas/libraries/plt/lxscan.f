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

C $Id: lxscan.f,v 1.1 1993/07/16 16:46:49 gdsjaar Exp $ 
C $Log: lxscan.f,v $
C Revision 1.1  1993/07/16 16:46:49  gdsjaar
C Changed plt to library rather than single source file.
C 
C=======================================================================
      LOGICAL FUNCTION LXSCAN(DELIM,REM,L,CH)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT
      CHARACTER*(*) DELIM
      CHARACTER*(*) REM
      INTEGER L
      CHARACTER CH

      REM = ' '
      L = 0
      LXSCAN = .FALSE.
 2410 CONTINUE
      CH = ILINE(JLINE:JLINE)
      IF (CH.EQ.CHAR(0)) THEN
         RETURN

      END IF

      JLINE = JLINE + 1
      LXSCAN = (INDEX(DELIM,CH).NE.0)
      IF (LXSCAN) THEN
         RETURN

      END IF

      L = L + 1
      IF (L.GT.LEN(REM)) THEN
         GO TO 2420

      END IF

      REM(L:L) = CH
 2420 GO TO 2410

      END
