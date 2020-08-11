C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
*DECK,CLSFIL
      SUBROUTINE CLSFIL

C     ******************************************************************

C     SUBROUTINE TO CLOSE PREVIOUSLY OPENED DISK FILES

C     Called by ERROR, MAPVAR

C     ******************************************************************

      include 'ex2tp.blk'
      include 'tapes.blk'
      include 'ntpdat.blk'

C     ******************************************************************

      IF (IFILES(1).EQ.1) CLOSE (UNIT=NTPOUT, STATUS='keep')
      IF (IFILES(3).EQ.1) call exclos (ntp2ex,ierr)
      IF (IFILES(4).EQ.1) call exclos (ntp3ex,ierr)
      IF (IFILES(5).EQ.1) call exclos (ntp4ex,ierr)

      RETURN
      END
