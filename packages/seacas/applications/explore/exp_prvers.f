C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE PRVERS (NDB, NOUT)
C=======================================================================
      include 'exodusII.inc'
      include 'exp_progqa.blk'

      character*1 cdum
      real libver

      call exinq(ndb, EXVERS, idum, apiver, cdum, ierr)
      call exinq(ndb, EXDBVR, idum, dbver,  cdum, ierr)
      call exinq(ndb, EXLBVR, idum, libver, cdum, ierr)

      IF (NOUT .GT. 0) WRITE (NOUT, 10000)

      CALL BANNER (NOUT, QAINFO,' ', ' ', ' ')
      IF (NOUT .GT. 0) THEN
         WRITE (NOUT, 10010) apiver, dbver, libver
      ELSE
         WRITE (*, 10010) apiver, dbver, libver
      END IF

10000 format(/, 1x, 'VERSION INFORMATION')
10010 format(/,
     &  1x, 'Database written with Exodus API version: ', F6.3,/,
     &  1x, 'Exodus database version:                  ', F6.3,/,
     &  1x, 'Exodus API Library version (explore using): ', F6.3)
      return
      end
