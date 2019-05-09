C Copyright(C) 2011-2017 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C
C * Redistributions of source code must retain the above copyright
C    notice, this list of conditions and the following disclaimer.
C
C * Redistributions in binary form must reproduce the above
C   copyright notice, this list of conditions and the following
C   disclaimer in the documentation and/or other materials provided
C   with the distribution.
C
C * Neither the name of NTESS nor the names of its
C   contributors may be used to endorse or promote products derived
C   from this software without specific prior written permission.
C
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
