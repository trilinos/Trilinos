C    Copyright(C) 2008 Sandia Corporation.  Under the terms of Contract
C    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C    certain rights in this software
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C    * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C              
C    * Redistributions in binary form must reproduce the above
C      copyright notice, this list of conditions and the following
C      disclaimer in the documentation and/or other materials provided
C      with the distribution.
C                            
C    * Neither the name of Sandia Corporation nor the names of its
C      contributors may be used to endorse or promote products derived
C      from this software without specific prior written permission.
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
C=======================================================================
      subroutine chkerr (routine, caller, ierr)
C=======================================================================
C     Modified 9/13/95 for EXODUSIIV2 API calls
C
C     This subroutine should be called after an EXODUSIIV2 subroutine has
C     been invoked.  The arguments of this subroutine are as follows:
C     routine - IN - The exodusIIv2 subroutine
C     caller  - IN - The subroutine invoking the exodusII call
C     ierr    - IN - The error code returned from the exodusII call
C     

      character*6 routine
      CHARACTER*8 caller
      character*(512) string, path, msg
      integer ierr

C     No error occurred
      if (ierr .eq. EXOK) RETURN

      write (path, '(A,A2,A)') CALLER(:LENSTR(CALLER)), '->', 
     &       ROUTINE(:LENSTR(ROUTINE))
      lp = lenstr(path)

      if (ierr .eq. EXWARN) then
         write (string, '(A,A)')'ExodusII V2 WARNING in ', path(:lp)
         call prterr('WARNING', string)
         RETURN
      else if (ierr .lt. EXOK) then
         write (string, '(A)')'ExodusII V2 ERROR'
         call prterr('ERROR', string)
      end if

      if (ierr .eq. EXFATL) then
         write (msg, '(A,A)')
     &   'A Fatal error occurred in ', path(:lp)
      else if (ierr .eq. EXMEMF) then
         write (msg, '(A,A)')
     &   'Memory allocation failure flag in ', path(:lp)
      else if (ierr .eq. EXBFMD) then
         write (msg, '(A,A)')
     &   'Bad file mode in ', path(:lp)
      else if (ierr .eq. EXBFID) then
         write (msg, '(A,A)') 
     &   'BAD file ID in ', path(:lp)
      else if (ierr .eq. EXBTID) then
         write (msg, '(A,A)') 
     &   'Property table lookup failed in ', path(:lp)
      else if (ierr .eq. EXBPRM) then
         write (msg, '(A,A)') 
     &   'Bad parameter passed in ', path(:lp)
      else
         write (msg, '(A,A,A,I4)')
     &   'Unknown Error in ', path(:lp), ': IERR =  ',ierr
      end if

      call prterr ('ERROR', msg)
      
      return
      end
