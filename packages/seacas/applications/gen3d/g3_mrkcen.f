C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MRKCEN (LINK, XN, NUMCOL, NEROW, IELCEN,
     &   ICOL1, ALLCEN)
C=======================================================================

C   --*** MRKCEN *** (GEN3D) Find center column for each element row
C   --   Written by Amy Gilkey - revised 04/26/88
C   --
C   --MRKCEN identifies whether the first column of each row is on the
C   --center.  This is true iff the first node in the row is located to
C   --the left of the minimum location of the second node in all rows.
C   --
C   --Parameters:
C   --   LINK - IN - the connectivity for the 2D elements, always 4 nodes
C   --   XN - IN - the coordinates, needed to find node ordering
C   --   NUMCOL - IN - the number of columns in the center blocks
C   --   NEROW - IN - the number of element rows in the center blocks
C   --   IELCEN - IN - the element numbers of the center elements
C   --      by column and row (column 1 is not necessarily the center)
C   --   ICOL1 - OUT - set to 1 if the first element in the row is on the
C   --      center, else set to 0
C   --   ALLCEN - OUT - true iff all rows are on the center column
C   --
C   --Common Variables:
C   --   Uses IX1, IX2, IX3, IX4 of /CENPAR/

      INCLUDE 'g3_cenpar.blk'

      INTEGER LINK(4,*)
      REAL XN(*)
C ASSERT: NUMCOL > 0
      INTEGER IELCEN(NUMCOL,*)
      INTEGER ICOL1(*)
      LOGICAL ALLCEN

C   --Find the minimum X for the second node in all rows

      XMIN = 1E36
      DO 10 IROW = 1, NEROW
         IEL = IELCEN(1,IROW)
         INP = LINK(IX2,IEL)
         XMIN = MIN (XMIN, XN(INP))
   10 CONTINUE

C   --Mark as on center if the first node in the row is to the left of the
C   --minimum found above

      ALLCEN = .TRUE.

      DO 20 IROW = 1, NEROW
         IEL = IELCEN(1,IROW)
         INP = LINK(IX1,IEL)
         IF (XN(INP) .LT. XMIN) THEN
            ICOL1(IROW) = 1
         ELSE
            ICOL1(IROW) = 0
            ALLCEN = .FALSE.
         END IF
   20 CONTINUE

      RETURN
      END
