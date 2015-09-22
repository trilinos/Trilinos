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
      LOGICAL FUNCTION FFMATC (IFLD, INTYP, CFIELD, MATCH, NLET)
C=======================================================================
C$Id: ffmatc.f,v 1.2 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: ffmatc.f,v $
CRevision 1.2  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:14:28  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:14:27  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:25  gdsjaar
c Initial revision
c 

C   --*** FFMATC *** (FFLIB) Parse free-field character string if match
C   --   Written by Amy Gilkey - revised 07/01/87
C   --
C   --FFMATC parses a character field and returns true (and increments
C   --IFLD) iff it is equal to the match string.  Only NLET letters must
C   --be in the input field to match, but if more letters are given, they
C   --must match the match string exactly.
C   --
C   --Parameters:
C   --   IFLD - IN/OUT - the index of the current field number, incremented
C   --      only if field matches
C   --   INTYP - IN - the input type from the free-field reader
C   --   CFIELD - IN - the character fields
C   --   MATCH - IN - the match string
C   --   NLET - IN - number of letters that must match

C   --Routines Called:
C   --   MATSTR - (STRLIB) Check if string matches

      INTEGER IFLD
      INTEGER INTYP(*)
      CHARACTER*(*) CFIELD(*)
      CHARACTER*(*) MATCH
      INTEGER NLET

      LOGICAL MATSTR

      IF (INTYP(IFLD) .GE. 0) THEN
         FFMATC = MATSTR (CFIELD(IFLD), MATCH, NLET)
         IF (FFMATC) IFLD = IFLD + 1
      ELSE
         FFMATC = .FALSE.
      END IF

      RETURN
      END
