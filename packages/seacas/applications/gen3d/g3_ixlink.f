C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE IXLINK (NCEN, IELROW, LINK, XN, YN, IROT)
C=======================================================================

C   --*** IXLINK *** (GEN3D) Get connectivity indices
C   --   Written by Amy Gilkey - revised 08/16/88
C   --
C   --IXLINK sets the connectivity indices IX1..4 for center elements.
C   --The indices are set so that:
C   --   LINK(IX1,i) holds the leftmost column, bottom row
C   --   LINK(IX2,i) holds the rightmost column, bottom row
C   --   LINK(IX3,i) holds the rightmost column, top row
C   --   LINK(IX4,i) holds the leftmost column, top row
C   --
C   --Parameters:
C   --   NCEN - IN - the number of center block elements
C   --   IELROW - IN - the indices of the center block elements
C   --   LINK - IN - the connectivity for the 2D elements, always 4 nodes
C   --   XN, YN - IN - the coordinates, needed to find node ordering
C   --   IROT - OUT - track link rotation for sset face update
C   --
C   --Common Variables:
C   --   Sets IX1, IX2, IX3, IX4 of /CENPAR/

      INCLUDE 'g3_cenpar.blk'

      INTEGER IELROW(*)
      INTEGER LINK(4,*)
      REAL XN(*), YN(*)
      INTEGER IROT(*)

      INTEGER LSV(4)

C   --Get the connectivity indices so that LINK(IX1,i) and LINK(IX4,i)
C   --are on the same column and LINK(IX2,i) and LINK(IX3,i) are on the
C   --next column

      DO 40 IXEL = 1, NCEN
         IEL = IELROW(IXEL)

         IX1 = 1
         DIFMAX = 0.0
         L2 = LINK(2,IEL)
         L3 = LINK(3,IEL)
         L4 = LINK(4,IEL)
         DO 10 ILINK = 1, 4
            L1 = LINK(ILINK,IEL)
            IF ((MAX (XN(L1), XN(L4)) .LT. MAX (XN(L2), XN(L3))) .AND.
     &         (MAX (YN(L1), YN(L2)) .LT. MAX (YN(L3), YN(L4)))) THEN
               DIF = (MAX (XN(L2), XN(L3)) - MAX (XN(L1), XN(L4)))
     &            + (MAX (YN(L3), YN(L4)) - MAX (YN(L1), YN(L2)))
               IF (DIF .GT. DIFMAX) THEN
                  IX1 = ILINK
                  DIFMAX = DIF
               END IF
            END IF
            L2 = L3
            L3 = L4
            L4 = LINK(ILINK,IEL)
   10    CONTINUE

         IROT(IEL) = IX1
         IF (IX1 .NE. 1) THEN
            DO 20 ILINK = 1, 4
               LSV(ILINK) = LINK(ILINK,IEL)
   20       CONTINUE
            IX = IX1
            DO 30 ILINK = 1, 4
               LINK(ILINK,IEL) = LSV(IX)
               IX = IX + 1
               IF (IX .GT. 4) IX = IX - 4
   30       CONTINUE
         END IF
   40 CONTINUE

      IX1 = 1
      IX2 = 2
      IX3 = 3
      IX4 = 4

      RETURN
      END
