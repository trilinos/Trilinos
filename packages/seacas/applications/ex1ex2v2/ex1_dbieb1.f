C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBIEB1 (NDB, OPTION, IELB, NUMELB, NUMLNK, NUMATR,
     &   LINK, ATRIB, NATRDM, NLNKDM, *)
C=======================================================================

C   --*** DBIEB1 *** (EXOLIB) Read database element block misc.
C   --   Written by Amy Gilkey - revised 10/14/87
C   --   Modified by Greg Sjaardema - 8/8/90
C   --      ---Removed MAX from Dimension statements, Added NATRDM, NLNKDM
C   --
C   --DBIEB1 reads the element block connectivity and attribute information
C   --from the database.  An error message is displayed if the end of file
C   --is read.
C   --
C   --Parameters:
C   --   NDB - IN - the database file
C   --   OPTION - IN - ' ' to not store, '*' to store all, else store options:
C   --      'C' to store connectivity
C   --      'A' to store attributes
C   --   IELB - IN - the element block number
C   --   NUMELB - IN - the number of elements in the block
C   --   NUMLNK - IN - the number of nodes per element;
C   --      negate to not store connectivity
C   --   NUMATR - IN - the number of attributes;
C   --      negate to not store attributes
C   --   LINK - OUT - the element connectivity for this block
C   --   ATRIB - OUT - the attributes for this block
C   --   NATRDM - IN - dimension of atrib array
C   --   NLNKDM - IN - dimension of link array
C   --   * - return statement if end of file or read error
C   --
C   --Database must be positioned at start of element block misc. information
C   --upon entry; upon exit at end of element block misc. information.

      INTEGER NDB
      CHARACTER*(*) OPTION
      INTEGER NUMELB, NUMLNK, NUMATR
      INTEGER LINK(NLNKDM, *)
      REAL ATRIB(NATRDM,*)

      CHARACTER*80 ERRMSG

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'C') .GT. 0)) THEN
         READ (NDB, END=100, ERR=100, IOSTAT=IERR)
     &      ((LINK(I,NE), I=1,NUMLNK), NE=1,NUMELB)
      ELSE
         READ (NDB, END=100, ERR=100, IOSTAT=IERR)
      END IF

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'A') .GT. 0)) THEN
         READ (NDB, END=110, ERR=110, IOSTAT=IERR)
     &      ((ATRIB(I,NE), I=1,NUMATR), NE=1,NUMELB)
      ELSE
         READ (NDB, END=110, ERR=110, IOSTAT=IERR)
      END IF

      RETURN

  100 CONTINUE
      WRITE (ERRMSG, '(A, I8)', IOSTAT=IDUM)
     &   'CONNECTIVITY for ELEMENT BLOCK', IELB
      GOTO 120
  110 CONTINUE
      WRITE (ERRMSG, '(A, I8)', IOSTAT=IDUM)
     &   'ATTRIBUTES for ELEMENT BLOCK', IELB
      GOTO 120
  120 CONTINUE
      CALL DBERR (IERR, ERRMSG)
      RETURN 1
      END
