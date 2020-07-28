C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBOELB (NDB, NELBS, NELBE,
     &   IDELB, NUMELB, NUMLNK, NUMATR, LINK, BLKTYP, ATRIB)
C=======================================================================

C   --*** DBOELB *** (EXOLIB) Write database element blocks
C   --
C   --DBOELB writes the element block information to the database.
C   --Some dynamic dimensioning is done.
C   --
C   --Parameters:
C   --   NDB - IN - the database file
C   --   NELBS, NELBE - IN - the number of first and last element blocks
C   --      to write
C   --   IDELB - IN - the element block IDs for each block
C   --   NUMELB - IN - the number of elements in each block
C   --   NUMLNK - IN - the number of nodes per element in each block
C   --   NUMATR - IN - the number of attributes in each block
C   --   LINK - IN - the connectivity for each block
C   --   ATRIB - IN - the attributes for each block
C   --
      include 'exodusII.inc'

      INTEGER NDB
      INTEGER IDELB(*)
      INTEGER NUMELB(*)
      INTEGER NUMLNK(*)
      INTEGER NUMATR(*)
      INTEGER LINK(*)
      REAL ATRIB(*)
      CHARACTER*(mxstln) BLKTYP(*)

      IELNK = 0
      IEATR = 0

      call expclb(ndb, idelb, blktyp, numelb, numlnk, numatr,
     *  .FALSE., ierr)
      IF (IERR .NE. 0) THEN
        CALL PRTERR ('FATAL', 'Error defining element blocks')
      END IF

      DO 100 NELB = NELBS, NELBE
         IELB  = NELB-NELBS+1
         ISLNK = IELNK + 1
         IELNK = IELNK + NUMLNK(IELB) * NUMELB(IELB)
         ISATR = IEATR + 1
         IEATR = IEATR + NUMATR(IELB) * NUMELB(IELB)

         if (numelb(ielb) .gt. 0 .and. numlnk(ielb) .gt. 0) then
           call expelc(ndb, idelb(ielb), link(islnk), ierr)
           IF (IERR .NE. 0) THEN
             CALL PRTERR ('FATAL',
     *         'Error writing element block connectivity')
           END IF
         end if

         if (numatr(ielb) .gt. 0) then
           call expeat(ndb, idelb(ielb), atrib(isatr), ierr)
           IF (IERR .NE. 0) THEN
             CALL PRTERR ('FATAL',
     *         'Error writing element block attributes')
           END IF
         end if
 100  CONTINUE

      RETURN
      END
