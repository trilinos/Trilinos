C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE RDNAM1 (NDB, TYPE, NBLK, NVAR, ISEVOK)
C=======================================================================

C   --*** RDNAM1 *** (EXOLIB) Internal to RDNAME
C   --
C   --RDNAM1 reads the element block variable truth table.
C   --An error message is displayed if the end of file is read.
C   --
C   --Parameters:
C   --   NDB - IN - the database number
C   --   NBLK - IN - the number of element blocks
C   --   NVAR - IN - the number of element variables
C   --   ISEVOK - OUT - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(i,j) is NOT 0

      include 'exodusII.inc'

      CHARACTER*1 TYPE
      INTEGER ISEVOK(NVAR,*)

      CHARACTER*80 ERRMSG

      CALL INIINT (NBLK*NVAR, 999, ISEVOK)

      if (nblk .gt. 0 .and. nvar .gt. 0) then
         if (type .eq. 'E') then
            call exgvtt (ndb, nblk, nvar, isevok, ierr)
         else if (type .eq. 'M') then
            call exgnstt (ndb, nblk, nvar, isevok, ierr)
         else if (type .eq. 'S') then
            call exgsstt (ndb, nblk, nvar, isevok, ierr)
         end if
        if (ierr .ne. 0) go to 100
      end if

      RETURN

 100  CONTINUE
      WRITE (ERRMSG, 10000) 'VARIABLE TRUTH TABLE'
      GOTO 110
 110  CONTINUE
      CALL WDBERR (IERR, ERRMSG)

      RETURN

10000 FORMAT (A)
      END
