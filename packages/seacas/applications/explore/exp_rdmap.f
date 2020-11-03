C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE RDMAP (NDB, NUMEL, MAPEL, ISEOF)
C=======================================================================

C   --*** RDMAP *** (EXPLORE) Read database element order map
C   --
C   --RDMAP reads the element order map from the database.  An error
C   --message is displayed if the end of file is read.
C   --
C   --Parameters:
C   --   NDB - IN - the database file
C   --   NUMEL - IN - the number of elements
C   --   MAPEL - OUT - the element order map
C   --   ISEOF - IN/OUT - set true if end of file read

      include 'exodusII.inc'

      INTEGER MAPEL(*)
      LOGICAL ISEOF

      CHARACTER*80 ERRMSG

      CALL INIINT (NUMEL, 0, MAPEL)

C ... Don't warn about no map stored in file
      call exopts (0, ierr1)
      call exgenm(ndb, mapel, ierr)
      call exopts (EXVRBS, ierr1)
      if (ierr .lt. 0) then
        go to 100
      else if (ierr .ne. 0) then
        do 10 i=1, numel
          mapel(i) = i
 10     continue
      end if

      RETURN

  100 CONTINUE
      WRITE (ERRMSG, 10000) 'ELEMENT ORDER MAP'
      GOTO 110
  110 CONTINUE
      CALL WDBERR (IERR, ERRMSG)
      ISEOF = .TRUE.

      RETURN

10000  FORMAT (A)
      END
