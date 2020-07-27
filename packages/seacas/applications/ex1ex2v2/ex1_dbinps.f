C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBINPS (NDB, OPTION, NUMNPS, LNPSNL,
     &   IDNPS, NNNPS, IXNNPS, LTNNPS, FACNPS, *)
C=======================================================================

C   --*** DBINPS *** (EXOLIB) Read database node sets
C   --   Written by Amy Gilkey - revised 10/14/87
C   --
C   --DBINPS reads the node set information from the database.
C   --An error message is displayed if the end of file is read.
C   --
C   --Parameters:
C   --   NDB - IN - the database file
C   --   OPTION - IN - ' ' to not store, '*' to store all, else store options:
C   --      'H' to store information about node sets
C   --      'N' to store node set nodes
C   --      'F' to store node set factors
C   --   NUMNPS - IN - the number of node sets
C   --   LNPSNL - IN - the length of the node sets node list
C   --   IDNPS - OUT - the node set ID for each set (if OPTION)
C   --   NNNPS - OUT - the number of nodes for each set (if OPTION)
C   --   IXNNPS - OUT - the index of the first node for each set (if OPTION)
C   --   LTNNPS - OUT - the nodes for all sets (if OPTION)
C   --   FACNPS - OUT - the distribution factors for all sets (if OPTION)
C   --   * - return statement if end of file or read error
C   --
C   --Database must be positioned at start of node set information
C   --upon entry; upon exit at end of node set information.

      INTEGER NDB
      CHARACTER*(*) OPTION
      INTEGER NUMNPS, LNPSNL
      INTEGER IDNPS(*)
      INTEGER NNNPS(*)
      INTEGER IXNNPS(*)
      INTEGER LTNNPS(*)
      REAL FACNPS(*)

      CHARACTER*80 ERRMSG

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'H') .GT. 0)) THEN
         READ (NDB, END=100, ERR=100, IOSTAT=IERR)
     &      (IDNPS(INPS), INPS=1,NUMNPS)
         READ (NDB, END=110, ERR=110, IOSTAT=IERR)
     &      (NNNPS(INPS), INPS=1,NUMNPS)
         READ (NDB, END=120, ERR=120, IOSTAT=IERR)
     &      (IXNNPS(INPS), INPS=1,NUMNPS)
      ELSE
         READ (NDB, END=100, ERR=100, IOSTAT=IERR)
         READ (NDB, END=110, ERR=110, IOSTAT=IERR)
         READ (NDB, END=120, ERR=120, IOSTAT=IERR)
      END IF

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'N') .GT. 0)) THEN
         READ (NDB, END=130, ERR=130, IOSTAT=IERR)
     &      (LTNNPS(NL), NL=1,LNPSNL)
      ELSE
         READ (NDB, END=130, ERR=130, IOSTAT=IERR)
      END IF

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'F') .GT. 0)) THEN
         READ (NDB, END=140, ERR=140, IOSTAT=IERR)
     &      (FACNPS(NL), NL=1,LNPSNL)
      ELSE
         READ (NDB, END=140, ERR=140, IOSTAT=IERR)
      END IF

      RETURN

  100 CONTINUE
      ERRMSG = 'NODE SET IDS'
      GOTO 150
  110 CONTINUE
      ERRMSG = 'NODE SET NUMBER OF NODES'
      GOTO 150
  120 CONTINUE
      ERRMSG = 'NODE SET INDICES'
      GOTO 150
  130 CONTINUE
      ERRMSG = 'NODE SET NODES'
      GOTO 150
  140 CONTINUE
      ERRMSG = 'NODE SET DISTRIBUTION FACTORS'
      GOTO 150
  150 CONTINUE
      CALL DBERR (IERR, ERRMSG)
      RETURN 1
      END
