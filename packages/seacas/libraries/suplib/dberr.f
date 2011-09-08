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
      SUBROUTINE DBERR (IOSTAT, ERRMSG)
C=======================================================================
C$Id: dberr.f,v 1.4 2009/03/25 12:46:01 gdsjaar Exp $
C$Log: dberr.f,v $
CRevision 1.4  2009/03/25 12:46:01  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.3  1991/09/30 20:08:33  gdsjaar
CIncreased error number format for Cray
C
c Revision 1.2  1991/02/04  08:34:36  gdsjaar
c Changed IOSTAT format to I3
c
c Revision 1.1.1.1  90/08/14  16:12:31  gdsjaar
c Testing
c 
c Revision 1.1  90/08/14  16:12:30  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:08  gdsjaar
c Initial revision
c 

C   --*** DBERR *** (EXOLIB) Display a database error message
C   --   Written by Amy Gilkey - revised 02/27/86
C   --
C   --DBERR displays a database read error message.
C   --
C   --Parameters:
C   --   IOSTAT - IN - the error status
C   --   ERRMSG - IN/OUT - the item being read when the error occurred;
C   --      compressed

C   --Routines Called:
C   --   SQZSTR - (STRLIB) Delete extra blanks from string

      CHARACTER*(*) ERRMSG

      CALL SQZSTR (ERRMSG, LSTR)
      WRITE (*, 10000) ERRMSG(:LSTR)
      IF (IOSTAT .LT. 0) THEN
         WRITE (*, 10010) IOSTAT, 'Unexpected end of file'
      ELSE IF (IOSTAT .EQ. 67) THEN
         WRITE (*, 10010) IOSTAT, 'Input record is too short'
      ELSE
         WRITE (*, 10010) IOSTAT
      END IF

      RETURN

10000  FORMAT (/, ' DATABASE ERROR - Reading ', A)
10010  FORMAT (3X, ' FORTRAN Error #', I6, :, ' - ', A)
      END
