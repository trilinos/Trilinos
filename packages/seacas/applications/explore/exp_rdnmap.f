C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE RDNMAP (NDB, NUMNP, MAPNO, ISEOF)
C=======================================================================

C   --*** RDNMAP *** (EXPLORE) Read database node order map
C   --
C   --RDMAP reads the node order map from the database.  An error
C   --message is displayed if the end of file is read.
C   --
C   --Parameters:
C   --   NDB - IN - the database file
C   --   NUMNP - IN - the number of nodes
C   --   MAPNO - OUT - the node order map
C   --   ISEOF - IN/OUT - set true if end of file read

      include 'exodusII.inc'

      INTEGER MAPNO(*)
      LOGICAL ISEOF

      CHARACTER*80 ERRMSG

      CALL INIINT (NUMNP, 0, MAPNO)

C ... Don't warn about no map stored in file
      call exopts (0, ierr1)
      call exgnnm(ndb, mapno, ierr)
      call exopts (EXVRBS, ierr1)
      if (ierr .lt. 0) then
        go to 100
      else if (ierr .ne. 0) then
        do 10 i=1, numnp
          mapno(i) = i
 10     continue
      end if

      RETURN

  100 CONTINUE
      WRITE (ERRMSG, 10000) 'NODE ORDER MAP'
      GOTO 110
  110 CONTINUE
      CALL WDBERR (IERR, ERRMSG)
      ISEOF = .TRUE.

      RETURN

10000  FORMAT (A)
      END
