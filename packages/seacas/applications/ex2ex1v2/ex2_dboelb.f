C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBOELB (NDB, NELBS, NELBE,
     &   IDELB, NUMELB, NUMLNK, NUMATR, LINK, ATRIB)
C=======================================================================

C   --*** DBOELB *** (EXOLIB) Write database element blocks
C   --   Written by Amy Gilkey - revised 10/12/87
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
C   --Database must be positioned at start of element block information
C   --upon entry; upon exit at end of element block information.

      INTEGER NDB
      INTEGER IDELB(*)
      INTEGER NUMELB(*)
      INTEGER NUMLNK(*)
      INTEGER NUMATR(*)
      INTEGER LINK(*)
      REAL ATRIB(*)

      IELNK = 0
      IEATR = 0

      DO 100 NELB = NELBS, NELBE
         IELB = NELB-NELBS+1

         WRITE (NDB)
     &      IDELB(IELB), NUMELB(IELB), NUMLNK(IELB), NUMATR(IELB)

         ISLNK = IELNK + 1
         IELNK = IELNK + NUMLNK(IELB) * NUMELB(IELB)
         ISATR = IEATR + 1
         IEATR = IEATR + NUMATR(IELB) * NUMELB(IELB)

         CALL DBOEB1 (NDB, NELB,
     &      NUMELB(IELB), NUMLNK(IELB), NUMATR(IELB),
     &      LINK(ISLNK), ATRIB(ISATR),
     $      MAX(NUMLNK(IELB),1), MAX(NUMATR(IELB),1))
  100 CONTINUE

      RETURN
      END
