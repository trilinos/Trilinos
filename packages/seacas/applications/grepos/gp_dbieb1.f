C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBIEBI (NDB, OPTION, IELB, NUMELB, NUMLNK, NUMATR,
     &                   LINK, ATRIB, NATRDM, NLNKDM, *)
C=======================================================================

C   --*** DBIEB1 *** (EXOLIB) Read database element block misc.
C   --
C   --DBIEB1 reads the element block connectivity and attribute information
C   --from the database.  An error message is displayed if the end of file
C   --is read.
C   --
C   --Parameters:
C   --   NDB    - IN  - the database file
C   --   OPTION - IN  - ' ' to not store, '*' to store all, else store options:
C   --                  'C' to store connectivity
C   --                  'A' to store attributes
C   --   IELB   - IN  - the element block number
C   --   NUMELB - IN  - the number of elements in the block
C   --   NUMLNK - IN  - the number of nodes per element;
C   --                  negate to not store connectivity
C   --   NUMATR - IN  - the number of attributes;
C   --                  negate to not store attributes
C   --   LINK   - OUT - the element connectivity for this block
C   --   ATRIB  - OUT - the attributes for this block
C   --   NATRDM - IN  - dimension of atrib array
C   --   NLNKDM - IN  - dimension of link array
C   --   *      - OUT - return statement if end of file or read error

      include 'exodusII.inc'
      INTEGER NDB
      CHARACTER*(*) OPTION
      INTEGER NUMELB, NUMLNK, NUMATR
      INTEGER LINK(NLNKDM, *)
      REAL ATRIB(NATRDM,*)

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'C') .GT. 0)) THEN
        if (numelb .gt. 0 .and. numlnk .gt. 0) then
          call exgelc(ndb, ielb, link, ierr)
        end if
      END IF

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'A') .GT. 0)) THEN
        if (numatr .gt. 0) then
           call exgeat(ndb, ielb, atrib, ierr)
        end if
      END IF

      RETURN

      END
