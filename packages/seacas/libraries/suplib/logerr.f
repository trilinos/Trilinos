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
      SUBROUTINE LOGERR (ERRTYP, ERRMSG, IUN)
C=======================================================================
C$Id: logerr.f,v 1.2 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: logerr.f,v $
CRevision 1.2  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:15:27  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:15:25  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:35  gdsjaar
c Initial revision
c 

C   --*** LOGERR *** (ETCLIB) Print error message to output file
C   --   Written by Amy Gilkey - revised 02/19/88
C   --
C   --Parameters:
C   --   ERRTYP - IN - the type of error:
C   --      'FATAL', 'PROGRAM', 'ERROR', 'WARNING',
C   --      'CMDERR', 'CMDWARN', 'CMDREQ', 'CMDSPEC'
C   --   ERRMSG - IN - the error message
C   --   IUN - IN - the output unit number

      CHARACTER*(*) ERRTYP
      CHARACTER*(*) ERRMSG

      IF (ERRTYP .EQ. 'FATAL') THEN
         WRITE (IUN, *)
         WRITE (IUN, 10) ERRMSG
   10    FORMAT (' FATAL ERROR - ', A)
      ELSE IF (ERRTYP .EQ. 'PROGRAM') THEN
         WRITE (IUN, *)
         WRITE (IUN, 20) ERRMSG, ' - call code sponsor'
   20    FORMAT (' PROGRAM ERROR - ', A, A)
      ELSE IF (ERRTYP .EQ. 'ERROR') THEN
         WRITE (IUN, *)
         WRITE (IUN, 30) ERRMSG
   30    FORMAT (' ERROR - ', A)
      ELSE IF (ERRTYP .EQ. 'WARNING') THEN
         WRITE (IUN, *)
         WRITE (IUN, 40) ERRMSG
   40    FORMAT (' WARNING - ', A)
      ELSE IF (ERRTYP .EQ. 'CMDERR') THEN
         WRITE (IUN, 50) ERRMSG
   50    FORMAT (' *** ERROR - ', A)
      ELSE IF (ERRTYP .EQ. 'CMDWARN') THEN
         WRITE (IUN, 60) ERRMSG
   60    FORMAT (' *** WARNING - ', A)
      ELSE IF (ERRTYP .EQ. 'CMDREQ') THEN
         WRITE (IUN, 70) ERRMSG
   70    FORMAT (' *** ', A)
      ELSE IF (ERRTYP .EQ. 'CMDSPEC') THEN
         WRITE (IUN, 80) ERRMSG
   80    FORMAT (' *** ', A)
      END IF

      RETURN
      END
