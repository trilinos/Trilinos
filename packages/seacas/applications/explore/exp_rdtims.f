C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE RDTIMS (NDB, TIMES)
C=======================================================================

C   --*** RDTIMS *** (EXPLORE) Read database time step times
C   --
C   --RDTIMS reads all time steps of the database, storing the
C   --time for each step.
C   --
C   --Parameters:
C   --   NDB - IN - the database number
C   --   TIMES - OUT - the time step times

      include 'exodusII.inc'
      REAL TIMES(*)

      CHARACTER*80 ERRMSG

      call exgatm(ndb, times, ierr)
      if (ierr .ne. 0) go to 140

      return

 140  CONTINUE
      WRITE (ERRMSG, 10010) 'TIME STEP TIME'
      CALL WDBERR (IERR, ERRMSG)

      RETURN

10010 FORMAT (A)
      END

