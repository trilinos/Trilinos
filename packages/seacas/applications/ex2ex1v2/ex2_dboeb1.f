C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBOEB1 (NDB, IELB, NUMELB, NUMLNK, NUMATR, LINK, ATRIB,
     $     NLNKDM, NATRDM)
C=======================================================================

C   --*** DBOEB1 *** (EXOLIB) Write database element block misc.
C   --   Written by Amy Gilkey - revised 10/12/87
C   --
C   --DBOEB1 writes the element block connectivity and attribute information
C   --to the database.
C   --
C   --Parameters:
C   --   NDB - IN - the database file
C   --   IELB - IN - the element block number
C   --   NUMELB - IN - the number of elements in the block
C   --   NUMLNK - IN - the number of nodes per element
C   --   NUMATR - IN - the number of attributes
C   --   LINK - IN - the element connectivity for this block
C   --   ATRIB - IN - the attributes for this block
C   --
C   --Database must be positioned at start of element block misc. information
C   --upon entry; upon exit at end of element block misc. information.

      INTEGER NDB
      INTEGER NUMELB, NUMLNK, NUMATR
      INTEGER LINK(NLNKDM,*)
      REAL ATRIB(NATRDM,*)

      IF ((NUMLNK .GT. 0) .AND. (NUMELB .GT. 0)) THEN
         WRITE (NDB) ((LINK(I,NE), I=1,NUMLNK), NE=1,NUMELB)
      ELSE
         WRITE (NDB) 0
      END IF

      IF ((NUMATR .GT. 0) .AND. (NUMELB .GT. 0)) THEN
         WRITE (NDB) ((ATRIB(I,NE), I=1,NUMATR), NE=1,NUMELB)
      ELSE
         WRITE (NDB) 0
      END IF

      RETURN
      END
