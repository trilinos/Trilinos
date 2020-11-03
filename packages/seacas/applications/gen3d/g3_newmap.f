C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE NEWMAP (MAPEL, MAPEL3, IXEL, INCEL, NREL, IELCOL)
C=======================================================================

C   --*** NEWMAP *** (GEN3D) Calculate 3D element order map
C   --   Written by Amy Gilkey - revised 05/05/86
C   --
C   --NEWMAP calculates the element order map for the 3D database.
C   --
C   --Parameters:
C   --   MAPEL - IN - the 2D element order map
C   --   MAPEL3 - OUT - the 3D element order map
C   --   IXEL - IN - the new index for each element
C   --   INCEL - IN - the increment for each element, needed for blocks
C   --      that become multiple blocks
C   --   NREL - IN - the number of new elements generated for each element
C   --   IELCOL - IN - the row number for each element, 0 if not needed
C   --
C   --Common Variables:
C   --   Uses NUMEL of /DBNUMS/
C   --   Uses NUMNP3 of /DBNUM3/
C   --   Uses NEREPL of /PARAMS/

      INCLUDE 'g3_dbnums.blk'
      INCLUDE 'g3_dbnum3.blk'
      INCLUDE 'g3_params.blk'

      INTEGER MAPEL(NUMEL), MAPEL3(NUMEL3)
      INTEGER IXEL(*), INCEL(*), NREL(*), IELCOL(*)

C   --Element order map - add on elements for each plate/slice
C ... To avoid potential numbering conflicts, we throw away the old map and just create a
C     dummy map from 1..numel3

      DO 20 I = 1, NUMEL3
        MAPEL3(I) = I
   20 CONTINUE

      RETURN
      END
