C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE REFND1 (IELB, NUMELB, NUMLNK, LINK, MAP)
C=======================================================================

C   --*** REFND1 *** (GREPOS) Reference nodes on all elements
C   --   Written by Greg Sjaardema - revised 03/02/90
C   --
C   --Parameters:
C   --   IELB - IN - the element block number
C   --   NUMELB - IN - the number of elements in the block
C   --   NUMLNK - IN - the number of nodes per element
C   --   LINK - IN/OUT - the element connectivity for this block
C   --   MAP - OUT - list of active nodes
C   --

      INTEGER NUMELB, NUMLNK
      INTEGER LINK(NUMLNK,*)
      INTEGER MAP(*)

      DO 20 ILNK = 1, NUMLNK
         DO 10 NE = 1, NUMELB
            MAP(LINK(ILNK, NE)) = 1
   10    CONTINUE
   20 CONTINUE

      RETURN
      END
