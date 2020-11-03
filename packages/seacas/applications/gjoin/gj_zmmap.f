C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE ZMMAP (NUMEL, MAPEL)
C=======================================================================

C   --*** ZMMAP *** (GJOIN) Compress element order map
C   --   Written by Amy Gilkey - revised 01/20/88
C   --
C   --ZMMAP compresses the element order map by removing deleted elements.
C   --
C   --Parameters:
C   --   NUMEL - IN/OUT - the number of elements
C   --   MAPEL - IN/OUT - the element order map

      INTEGER MAPEL(*)

      JEL = 0
      DO 100 IEL = 1, NUMEL
         IF (MAPEL(IEL) .GT. 0) THEN
            JEL = JEL + 1
            MAPEL(JEL) = MAPEL(IEL)
         END IF
  100 CONTINUE

      NUMEL = JEL

      RETURN
      END
