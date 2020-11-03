C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE RENELB (NELBLK, IOFFNP, IXNP, NUMELB, NUMLNK, LINK)
C=======================================================================

C   --*** RENELB *** (GJOIN) Renumber connectivity in element blocks
C   --   Written by Amy Gilkey - revised 09/29/87
C   --
C   --RENELB renumbers the nodes in the connectivity arrays.
C   --
C   --Parameters:
C   --   NELBLK - IN - the number of element blocks
C   --   IOFFNP - IN - the nodal offset: if positive, add to node number;
C   --      if negative, use IXNP
C   --   IXNP - IN - the new node number for each node
C   --   NUMELB - IN - the number of elements in each block
C   --   NUMLNK - IN - the number of nodes per element in each block
C   --   LINK - IN - the connectivity for each block

      INTEGER IXNP(*)
      INTEGER NUMELB(*)
      INTEGER NUMLNK(*)
      INTEGER LINK(*)

      ILNK = 1

      DO 100 IELB = 1, NELBLK
C ... FALSE in call to RENIX means map old node to new instead of
C     deleting the node if IXNP(I) < 0
         CALL RENIX (NUMELB(IELB) * NUMLNK(IELB),
     &      IOFFNP, IXNP, LINK(ILNK), .FALSE.)
         ILNK = ILNK + NUMLNK(IELB) * NUMELB(IELB)
  100 CONTINUE

      RETURN
      END
