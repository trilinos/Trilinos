C Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software.
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

C=======================================================================
      INTEGER FUNCTION LENSTR (STRING)
C=======================================================================
C$Id: lenstr.f,v 1.2 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: lenstr.f,v $
CRevision 1.2  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:15:16  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:15:15  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:33  gdsjaar
c Initial revision
c 

C   --*** LENSTR *** (STRLIB) Return string length
C   --   Written by Amy Gilkey - revised 02/14/86
C   --
C   --LENSTR returns the length of a passed string.  There can be blanks
C   --embedded within the string.  To prevent problems, an empty string
C   --is returned with a length of 1.
C   --
C   --Parameters:
C   --   STRING - IN - the string to find the length of

      CHARACTER*(*) STRING

      LENSTR = LEN(STRING)
      IF (STRING(LENSTR:LENSTR) .EQ. ' ') THEN
         N = INDEX (STRING, ' ')
   10    CONTINUE
         IF (STRING(N:) .NE. ' ') THEN
            N = N + INDEX (STRING(N+1:), ' ')
            GOTO 10
         END IF
         LENSTR = MAX (1, N-1)
      END IF

      RETURN
      END
