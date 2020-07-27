C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBIINI (NDB, OPTION, NVERS, TITLE,
     &   NDIM, NUMNP, NUMEL, NELBLK,
     &   NUMNPS, LNPSNL, NUMESS, LESSEL, LESSNL, *)
C=======================================================================

C   --*** DBIINI *** (EXOLIB) Read database title and initial variables
C   --   Written by Amy Gilkey - revised 05/24/88
C   --
C   --DBIINI rewinds the database and reads the title and the initial
C   --variables from the database.  An error message is displayed if
C   --the end of file is read.
C   --
C   --Parameters:
C   --   NDB - IN - the database file
C   --   OPTION - IN - ' ' to not store, '*' to store all, else store options:
C   --      'T' to store title
C   --      'I' to store initial variables
C   --   NVERS - OUT - the version number (not read, always 1)
C   --   TITLE - OUT - the database title (if OPTION)
C   --   NDIM - OUT - the number of coordinates per node
C   --   NUMNP - OUT - the number of nodes
C   --   NUMEL - OUT - the number of elements
C   --   NELBLK - OUT - the number of element blocks
C   --   NUMNPS - OUT - the number of node sets
C   --   LNPSNL - OUT - the length of the node sets node list
C   --   NUMESS - OUT - the number of side sets
C   --   LESSEL - OUT - the length of the side sets element list
C   --   LESSNL - OUT - the length of the side sets node list
C   --   * - return statement if end of file or read error
C   --
C   --Database is rewound upon entry; upon exit positioned at end of
C   --initial variables.

      INTEGER NDB
      CHARACTER*(*) OPTION
      INTEGER NVERS
      CHARACTER*80 TITLE
      INTEGER NDIM, NUMNP, NUMEL, NELBLK,
     &   NUMNPS, LNPSNL, NUMESS, LESSEL, LESSNL

      CHARACTER*80 ERRMSG

      REWIND (NDB)

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'T') .GT. 0)) THEN
         READ (NDB, END=100, ERR=100, IOSTAT=IERR) TITLE
      ELSE
         READ (NDB, END=100, ERR=100, IOSTAT=IERR)
      END IF

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'I') .GT. 0)) THEN
         READ (NDB, END=110, ERR=110, IOSTAT=IERR)
     &      NUMNP, NDIM, NUMEL, NELBLK,
     &      NUMNPS, LNPSNL, NUMESS, LESSEL, LESSNL
         NVERS = 1
      ELSE
         READ (NDB, END=110, ERR=110, IOSTAT=IERR)
      END IF

      RETURN

  100 CONTINUE
      ERRMSG = 'TITLE'
      GOTO 120
  110 CONTINUE
      ERRMSG = 'INITIAL VARIABLES'
      GOTO 120
  120 CONTINUE
      CALL DBERR (IERR, ERRMSG)
      RETURN 1
      END
