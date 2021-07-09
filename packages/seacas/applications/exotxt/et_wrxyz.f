C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE WRXYZ (NTXT, NDIM, NUMNP, XN, YN, ZN, nameco, namlen)
C=======================================================================

C   --*** WRXYZ *** (EXOTXT) Write database coordinates
C   --   Written by Amy Gilkey - revised 02/27/86
C   --
C   --WRXYZ writes the coordinate array from the database.
C   --
C   --Parameters:
C   --   NTXT - IN - the text file
C   --   NDIM - IN - the number of coordinates per node
C   --   NUMNP - IN - the number of nodes
C   --   XN, YN, ZN - IN - the coordinates
C   --
C   --Database must be positioned at start of coordinates upon entry;
C   --upon exit at end of coordinates.

      REAL XN(*), YN(*), ZN(*)
      character*(namlen) nameco(*)

      WRITE (NTXT, '(A)')
     &   '! Coordinate names'
      write (ntxt, 10020) (nameco(i),i=1,ndim)

      WRITE (NTXT, '(A)') '! Coordinates'

      DO 100 INP = 1, NUMNP
         IF (NDIM .EQ. 1) THEN
            WRITE (NTXT, 10000) XN(INP)
         ELSE IF (NDIM .EQ. 2) THEN
            WRITE (NTXT, 10000) XN(INP), YN(INP)
         ELSE IF (NDIM .EQ. 3) THEN
            WRITE (NTXT, 10000) XN(INP), YN(INP), ZN(INP)
         END IF
  100 CONTINUE

      RETURN
10000  FORMAT (5(1pE16.7))
10020  FORMAT (3(A,1x))
       END
