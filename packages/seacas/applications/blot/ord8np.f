C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE ORD8NP (NELBLK, LENE, NLNKE, LINKE)
C=======================================================================

C   --*** ORD8NP *** (BLOT) Order 8-node 2D elements
C   --   Modified by John H. Glick -- 10/25/88
C   --   Written by Amy Gilkey, revised 10/29/87
C   --
C   --ORD8NP orders the nodes of 8 and 9 node 2D elements so that they are
C   --connected.  Normally nodes are ordered 1-2-3-4-5-6-7-8-9, where nodes
C   --1 to 4 are corner nodes, 5 to 8 are side nodes, and 9 is the center
C   --node.  The nodes are returned as 1-5-2-6-3-7-4-8-9.
C   --
C   --Parameters:
C   --   NELBLK - IN - the number of element blocks
C   --   LENE - IN - the cumulative element counts by element block
C   --   NLNKE - IN - the number of nodes per element
C   --   LINK - IN/OUT - the element connectivity

      INTEGER LENE(0:*), LINKE(*)
      INTEGER NLNKE(*)

      INTEGER L(9)

      DO 140 IELB = 1,  NELBLK
         IF ((NLNKE(IELB) .EQ. 8)) THEN
            IXL0 = IDBLNK (IELB, 0, LENE, NLNKE) - 1
            DO 110 IEL = LENE(IELB-1)+1, LENE(IELB)
               L(1) = LINKE(IXL0+1)
               L(2) = LINKE(IXL0+5)
               L(3) = LINKE(IXL0+2)
               L(4) = LINKE(IXL0+6)
               L(5) = LINKE(IXL0+3)
               L(6) = LINKE(IXL0+7)
               L(7) = LINKE(IXL0+4)
               L(8) = LINKE(IXL0+8)
               DO 100 K = 1, NLNKE(IELB)
                  LINKE(IXL0+K) = L(K)
  100          CONTINUE
               IXL0 = IXL0 + NLNKE(IELB)
  110       CONTINUE
         ELSE IF ((NLNKE(IELB) .EQ. 9)) THEN
            IXL0 = IDBLNK (IELB, 0, LENE, NLNKE) - 1
            DO 130 IEL = LENE(IELB-1)+1, LENE(IELB)
               L(1) = LINKE(IXL0+1)
               L(2) = LINKE(IXL0+5)
               L(3) = LINKE(IXL0+2)
               L(4) = LINKE(IXL0+6)
               L(5) = LINKE(IXL0+3)
               L(6) = LINKE(IXL0+7)
               L(7) = LINKE(IXL0+4)
               L(8) = LINKE(IXL0+8)
               L(9) = LINKE(IXL0+9)
               DO 120 K = 1, NLNKE(IELB)
                  LINKE(IXL0+K) = L(K)
  120          CONTINUE
               IXL0 = IXL0 + NLNKE(IELB)
  130       CONTINUE
         END IF
  140 CONTINUE

      RETURN
      END
