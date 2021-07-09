C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      subroutine chkerr (routine, caller, ierr)
C=======================================================================
C     Modified 9/13/95 for EXODUSIIV2 API calls

C     This subroutine should be called after an EXODUSIIV2 subroutine has
C     been invoked.  The arguments of this subroutine are as follows:
C     routine - IN - The exodusIIv2 subroutine
C     caller  - IN - The subroutine invoking the exodusII call
C     ierr    - IN - The error code returned from the exodusII call

      include 'exodusII.inc'

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
