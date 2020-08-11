C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBIXYZ (NDB, OPTION, NDIM, NUMNP, XN, YN, ZN, *)
C=======================================================================

C   --*** DBIXYZ *** (EXOLIB) Read database coordinates
C   --   Written by Amy Gilkey - revised 10/14/87
C   --
C   --DBIXYZ reads the coordinate array from the database.  An error
C   --message is displayed if the end of file is read.
C   --
C   --Parameters:
C   --   NDB - IN - the database file
C   --   OPTION - IN - ' ' to not store, '*' to store all, else store options:
C   --      'C' to store coordinates
C   --   NDIM - IN - the number of coordinates per node
C   --   NUMNP - IN - the number of nodes
C   --   XN, YN, ZN - OUT - the coordinates (if OPTION)
C   --   * - return statement if end of file or read error
C   --
C   --Database must be positioned at start of coordinates upon entry;
C   --upon exit at end of coordinates.

      INTEGER NDB
      CHARACTER*(*) OPTION
      INTEGER NDIM, NUMNP
      REAL XN(*), YN(*), ZN(*)

      CHARACTER*80 ERRMSG

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'C') .GT. 0)) THEN
         IF (NDIM .EQ. 1) THEN
            READ (NDB, END=100, ERR=100, IOSTAT=IERR)
     &         (XN(INP), INP=1,NUMNP)
         ELSE IF (NDIM .EQ. 2) THEN
            READ (NDB, END=100, ERR=100, IOSTAT=IERR)
     &         (XN(INP), INP=1,NUMNP), (YN(INP), INP=1,NUMNP)
         ELSE IF (NDIM .EQ. 3) THEN
            READ (NDB, END=100, ERR=100, IOSTAT=IERR)
     &         (XN(INP), INP=1,NUMNP), (YN(INP), INP=1,NUMNP),
     &         (ZN(INP), INP=1,NUMNP)
         ELSE
            READ (NDB, END=100, ERR=100, IOSTAT=IERR)
         END IF
      ELSE
         READ (NDB, END=100, ERR=100, IOSTAT=IERR)
      END IF

      RETURN

  100 CONTINUE
      ERRMSG = 'NODAL COORDINATES'
      CALL DBERR (IERR, ERRMSG)
      RETURN 1
      END
