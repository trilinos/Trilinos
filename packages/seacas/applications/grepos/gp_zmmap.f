C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE ZMMAP (NUMEL, MAPEL)
C=======================================================================

C   --*** ZMMAP *** (GJOIN) Compress element and node order maps
C   --   Written by Amy Gilkey - revised 01/20/88
C   --
C   --ZMMAP compresses the element/node order map by removing deleted
C     elements/nodes.
C   --
C   --Parameters:
C   --   NUMEL - IN - the number of elements/nodes
C   --   MAPEL - IN/OUT - the element/node order map

      INTEGER MAPEL(*)

      JEL = 0
      DO 100 IEL = 1, NUMEL
         IF (MAPEL(IEL) .GT. 0) THEN
            JEL = JEL + 1
            MAPEL(JEL) = MAPEL(IEL)
         END IF
  100 CONTINUE
      RETURN
      END
