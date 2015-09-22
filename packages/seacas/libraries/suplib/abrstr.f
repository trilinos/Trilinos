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
      SUBROUTINE ABRSTR (RETWRD, ABBR, STRTBL)
C=======================================================================
C$Id: abrstr.f,v 1.2 2009/03/25 12:46:01 gdsjaar Exp $
C$Log: abrstr.f,v $
CRevision 1.2  2009/03/25 12:46:01  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:11:58  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:11:57  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:03  gdsjaar
c Initial revision
c 

C   --*** ABRSTR *** (STRLIB) Find abbreviation for string
C   --   Written by Amy Gilkey - revised 02/14/86
C   --
C   --ABRSTR returns the non-abbreviated form of the given abbreviation
C   --from the list of possible strings.  The abbreviation must either
C   --be a complete string or it must only match one string.
C   --
C   --Parameters:
C   --   RETWRD - OUT - the string for the abbreviation; ' ' if none
C   --   ABBR - IN - the abbreviation
C   --   STRTBL - IN - the table of possible strings; ended by ' '

      CHARACTER*(*) RETWRD
      CHARACTER*(*) ABBR
      CHARACTER*(*) STRTBL(*)

      RETWRD = ' '

      IF (ABBR .EQ. ' ') RETURN

      L = INDEX (ABBR, ' ') - 1
      IF (L .LT. 0) L = LEN(ABBR)

      NFOUND = 0
      I = 1
  100 CONTINUE
      IF (STRTBL(I) .NE. ' ') THEN
         IF (ABBR .EQ. STRTBL(I)(1:L)) THEN
            RETWRD = STRTBL(I)
            IF (ABBR .EQ. STRTBL(I)) GOTO 110
            NFOUND = NFOUND + 1
         END IF
         I = I + 1
         GOTO 100
      END IF

      IF (NFOUND .GT. 1) RETWRD = ' '

  110 CONTINUE
      RETURN
      END
