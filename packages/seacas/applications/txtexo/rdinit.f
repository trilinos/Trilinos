C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE RDINIT (NTXT, VERS, TITLE, NDIM, NUMNP, NUMEL, NELBLK,
     &   NUMNPS, LNPSNL, LNPSDF, NUMESS, LESSEL, LESSNL, LESSDF, NAMLEN,
     *  *)
C=======================================================================

C   --*** RDINIT *** (TXTEXO) Read database title and initial variables
C   --   Written by Amy Gilkey - revised 12/16/87
C   --
C   --RDINIT reads the title and the initial variables from the text file.
C   --An error message is displayed if the end of file is read.
C   --
C   --Parameters:
C   --   NTXT - IN - the text file
C   --   NVERS - OUT - the version number
C   --   TITLE - OUT - the database title
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
C   --Database must be rewound before entry; upon exit positioned at end of
C   --initial variables.

      CHARACTER*80 TITLE
      CHARACTER*80 SCRATCH

      READ (NTXT, *, END=110, ERR=110)
      READ (NTXT, '(A)', END=110, ERR=110) TITLE

      READ (NTXT, *, END=120, ERR=120)

C ... Split scratch string into pieces.  The namlen does not exist on older
C     databases, so we need to see if it is there or not.  The others are always there.
      READ (NTXT, '(A)', END=110, ERR=110) SCRATCH
      READ(SCRATCH( 1:10),*) NDIM
      READ(SCRATCH(11:20),*) VERS
      if (scratch(21:30) .ne. '          ') then
        READ(SCRATCH(21:30),*) NAMLEN
      else
        namlen = 32
      end if

      READ (NTXT, *, END=130, ERR=130) NUMNP, NUMEL, NELBLK
      READ (NTXT, *, END=150, ERR=150) NUMNPS, NUMESS
      READ (NTXT, *, END=160, ERR=160) LNPSNL, LNPSDF
      READ (NTXT, *, END=170, ERR=170) LESSEL, LESSNL, LESSDF

      RETURN

  110 CONTINUE
      CALL PRTERR ('FATAL', 'Reading TITLE')
      GOTO 180
  120 CONTINUE
      CALL PRTERR ('FATAL', 'Reading NUMBER OF DIMENSIONS or VERSION')
      GOTO 180
  130 CONTINUE
      CALL PRTERR ('FATAL',
     &     'Reading NUMBER OF NODES, ELEMENTS, and ELEMENT BLOCKS')
      GOTO 180
  150 CONTINUE
      CALL PRTERR ('FATAL',
     &   'Reading NUMBER OF NODE SETS and SIDE SETS')
      GOTO 180
  160 CONTINUE
      CALL PRTERR ('FATAL', 'Reading NODE SET NODE LENGTHS')
      GOTO 180
  170 CONTINUE
      CALL PRTERR ('FATAL',
     &     'Reading SIDE SET ELEMENT, NODE, and FACTOR COUNTS')
      GOTO 180
  180 CONTINUE
      RETURN 1
      END
