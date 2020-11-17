C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBOXYZ (NDB, NDIM, NUMNP, XN, YN, ZN)
C=======================================================================

C   --*** DBOXYZ *** (EXOLIB) Write database coordinates
C   --   Written by Amy Gilkey - revised 02/27/86
C   --
C   --DBOXYZ writes the coordinate array to the database.
C   --
C   --Parameters:
C   --   NDB - IN - the database file
C   --   NDIM - IN - the number of coordinates per node
C   --   NUMNP - IN - the number of nodes
C   --   XN, YN, ZN - IN - the coordinates
C   --
C   --Database must be positioned at start of coordinates upon entry;
C   --upon exit at end of coordinates.

      INTEGER NDB
      INTEGER NDIM, NUMNP
      REAL XN(*), YN(*), ZN(*)

      IF ((NUMNP .GT. 0) .AND. (NDIM .GT. 0)) THEN
         IF (NDIM .EQ. 1) THEN
            WRITE (NDB) (XN(INP), INP=1,NUMNP)
         ELSE IF (NDIM .EQ. 2) THEN
            WRITE (NDB) (XN(INP), INP=1,NUMNP), (YN(INP), INP=1,NUMNP)
         ELSE IF (NDIM .GE. 3) THEN
            WRITE (NDB) (XN(INP), INP=1,NUMNP), (YN(INP), INP=1,NUMNP),
     &         (ZN(INP), INP=1,NUMNP)
         ELSE
            WRITE (NDB) 0
         END IF
      ELSE
         WRITE (NDB) 0
      END IF

      RETURN
      END
