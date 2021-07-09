C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE RWXYZ (NTXT, NDB, NDIM, NUMNP, XN, YN, ZN, NAMECO,
     *  NAMLEN,*)
C=======================================================================

C   --*** RDXYZ *** (TXTEXO) Read database coordinates
C   --   Written by Amy Gilkey - revised 02/27/86
C   --
C   --RDXYZ reads the coordinate array from the text file.  An error
C   --message is displayed if the end of file is read.
C   --
C   --Parameters:
C   --   NTXT - IN - the text file
C   --   NDIM - IN - the number of coordinates per node
C   --   NUMNP - IN - the number of nodes
C   --   XN, YN, ZN - OUT - the coordinates
C   --   * - return statement if end of file or read error
C   --
C   --Database must be positioned at start of coordinates upon entry;
C   --upon exit at end of coordinates.

      REAL XN(*), YN(*), ZN(*)
      CHARACTER*(NAMLEN) NAMECO(*)
      integer kval(3)
      real rval(3)
      integer ival(3)

      character*512 scratch
      CHARACTER*32 STRA

      INP = 0
      nfield = 3
      READ (NTXT, *, END=110, ERR=110)
      READ (ntxt, '(A)', END=110, ERR=110) scratch
      idcont = 0
      call ffistr (scratch, ndim, idcont, nfield, kval, nameco,
     *  ival, rval)

      READ (NTXT, *, END=120, ERR=120)
      DO 100 INP = 1, NUMNP
         IF (NDIM .EQ. 1) THEN
            READ (NTXT, *, END=120, ERR=120) XN(INP)
         ELSE IF (NDIM .EQ. 2) THEN
            READ (NTXT, *, END=120, ERR=120) XN(INP), YN(INP)
         ELSE IF (NDIM .EQ. 3) THEN
            READ (NTXT, *, END=120, ERR=120) XN(INP), YN(INP), ZN(INP)
         END IF
  100 CONTINUE

      call expcon (ndb, nameco, ierr)
      call expcor (ndb, xn, yn, zn, ierr)

      RETURN

 110  continue
      CALL PRTERR ('FATAL', 'Reading COORDINATE NAMES')
      RETURN 1

  120 CONTINUE
      CALL INTSTR (1, 0, INP, STRA, LSTRA)
      CALL PRTERR ('FATAL',
     &   'Reading COORDINATES for node ' // STRA(:LSTRA))
      RETURN 1
      END
