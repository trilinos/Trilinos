C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      INTEGER FUNCTION IOWDSZ()
C=======================================================================
C convenience function for setting desired exodusII output word size
C via the EXT05 environment variable.

      character*80 ws
      character*8 cdum

      iows = 4
      call exname (-5, ws, llen)
      if (llen .lt. 1) goto 25
      read(ws,'(i1)',ERR=25)iows
      if (iows .ne. 4 .and. iows .ne. 8) then
        CALL PRTERR('WARNING', 'invalid output word size, set to 4')
        iows = 4
      endif
      goto 26
25    continue

      call exparm (cdum, cdum, idum, iows, idum, idum)
26    continue
      iowdsz = iows
      return
      end
