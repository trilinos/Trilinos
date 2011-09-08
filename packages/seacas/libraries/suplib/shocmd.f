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
      SUBROUTINE SHOCMD (HEADER, LIST)
C=======================================================================
C$Id: shocmd.f,v 1.2 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: shocmd.f,v $
CRevision 1.2  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:16:19  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:16:18  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:45  gdsjaar
c Initial revision
c 

C   --*** SHOCMD *** (ETCLIB) Display list of strings
C   --   Written by Amy Gilkey - revised 10/27/86
C   --
C   --SHOCMD displays the list of strings.  The heading for the list is
C   --dependent on HEADER:
C   --   COMMANDS - "Valid Commands"
C   --   null - no heading
C   --   other - HEADER
C   --
C   --Parameters:
C   --   HEADER - IN - the show option (see above)
C   --   LIST - IN - the string list, last entry must be ' '

C   --Routines Called:
C   --   LOCSTR - (STRLIB) Find string

      CHARACTER*(*) HEADER
      CHARACTER*(*) LIST(*)

      IF (HEADER .EQ. 'COMMANDS') THEN
         WRITE (*, 10) 'Valid Commands:'
      ELSE IF (HEADER .NE. ' ') THEN
         WRITE (*, 10) HEADER
      END IF

      NCMD = LOCSTR (' ', 999, LIST) - 1
      WRITE (*, 20) (LIST(I), I=1,NCMD)

      RETURN
   10 FORMAT (/, 1X, 5A)
   20 FORMAT ((4X, 7(A8, :, 2X)))
      END
