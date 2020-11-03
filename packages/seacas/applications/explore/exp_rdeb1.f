C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE RDEB1 (NDB, IDELB, NUMELB, NUMLNK, NUMATR,
     &   LINK, ATRIB, ATRNM, NAMLEN)
C=======================================================================

C   --*** RDEB1 *** (EXPLORE) Read database element block misc.
C   --
C   --RDEB1 reads the element block connectivity and attribute information
C   --from the database.  An error message is displayed if the end of file
C   --is read.
C   --
C   --Parameters:
C   --   NDB - IN - the database file
C   --   IDELB - IN - the element block number (for errors)
C   --   NUMELB - IN - the number of elements in the block
C   --   NUMLNK - IN - the number of nodes per element
C   --   NUMATR - IN - the number of attributes
C   --   LINK - OUT - the element connectivity for this block
C   --   ATRIB - OUT - the attributes for this block
C   --
      include 'exodusII.inc'

      INTEGER LINK(*)
      REAL ATRIB(*)

      CHARACTER*80 ERRMSG
      CHARACTER*(NAMLEN) ATRNM(*)

      IF (NUMELB .GT. 0 .AND. NUMLNK .GT. 0) THEN
        CALL EXGELC (NDB, IDELB, LINK, IERR)
        IF (IERR .NE. 0) GO TO 100
      END IF

      IF (NUMELB .GT. 0 .AND. NUMATR .GT. 0) THEN
        CALL EXGEAT (NDB, IDELB, ATRIB, IERR)
        IF (IERR .NE. 0) GO TO 110

        CALL EXGEAN (NDB, IDELB, NUMATR, ATRNM, IERR)
        IF (IERR .NE. 0) GO TO 115
      END IF

      RETURN

  100 CONTINUE
      WRITE (ERRMSG, 10000, IOSTAT=IDUM)
     *  'CONNECTIVITY for block', IDELB
      GOTO 120
  110 CONTINUE
      WRITE (ERRMSG, 10000, IOSTAT=IDUM)
     *  'ATTRIBUTES for block', IDELB
      GOTO 120
 115  CONTINUE
      WRITE (ERRMSG, 10000, IOSTAT=IDUM)
     *  'ATTRIBUTE NAMES for block', IDELB
      GOTO 120
  120 CONTINUE
      CALL WDBERR (IERR, ERRMSG)

      RETURN

10000  FORMAT (5 (A, I12))
      END
