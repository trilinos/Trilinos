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
      SUBROUTINE MAKROW (LINK, XN, NUMCEN, NUMCOL,
     &   NEROW, IELROW, IXROW)
C=======================================================================

C   --*** MAKROW *** (GEN3D) Order center elements in rows
C   --   Written by Amy Gilkey - revised 04/26/88
C   --
C   --MAKROW orders the elements in all center blocks into rows.  Only the
C   --elements in the needed number of columns for each row are stored.
C   --
C   --Parameters:
C   --   LINK - IN - the connectivity for the 2D elements, always 4 nodes
C   --   XN - IN - the coordinates, needed to find node ordering
C   --   NUMCEN - IN - the number of elements in the center blocks
C   --   NUMCOL - IN - the number of columns in the center block
C   --   NEROW - OUT - the number of element rows in the center blocks
C   --   IELROW - IN/OUT - the element numbers of the center elements;
C   --      returned as rows of elements
C   --   IXROW - OUT - the IELROW index of the starting column for each row
C   --
C   --Common Variables:
C   --   Uses IX1, IX2, IX3, IX4 of /CENPAR/

      INCLUDE 'g3_cenpar.blk'

      INTEGER LINK(4,*)
      REAL XN(*)
      INTEGER IELROW(*)
      INTEGER IXROW(*)

C   --Initialize the IELROW array
C   --   IELROW(1..ICOL1-1) holds the ordered elements for previous rows
C   --   IELROW(ICOL1..ICOL+NCOL-1) holds the ordered elements for this row
C   --   IELROW(IFILL..NUMCEN) holds the elements yet to be ordered

      IFILL = 1
      ICOL1 = 1

   10 CONTINUE
      IF (IFILL .LE. NUMCEN) THEN

C      --Pick leftmost unordered element for the next row

         ICHK = IFILL
         XMIN = XN(LINK(IX1,IELROW(IFILL)))
         DO 20 I = IFILL+1, NUMCEN
            IF (XN(LINK(IX1,IELROW(I))) .LT. XMIN) THEN
               ICHK = I
               XMIN = XN(LINK(IX1,IELROW(I)))
            END IF
   20    CONTINUE

         NEROW = NEROW + 1
         IC = IELROW(ICHK)
         IELROW(ICHK) = IELROW(IFILL)
         IFILL = IFILL + 1
         IELROW(ICOL1) = IC
         IXROW(NEROW) = ICOL1
         NCOL = 1

c      --Find elements to the left of the element until at leftmost column !#VAX

C      --Find elements to the right of the element until at rightmost column

   40    CONTINUE
         IEL = IELROW(ICOL1+NCOL-1)
         L2 = LINK(IX2,IEL)
         L3 = LINK(IX3,IEL)
         DO 50 ICHK = IFILL, NUMCEN
            IC = IELROW(ICHK)
            IF ((L2 .EQ. LINK(IX1,IC))
     &         .AND. (L3 .EQ. LINK(IX4,IC))) THEN
               IELROW(ICHK) = IELROW(IFILL)
               IFILL = IFILL + 1
               IELROW(ICOL1+NCOL) = IC
               NCOL = NCOL + 1
               GOTO 40
            END IF
   50    CONTINUE

         NCOL = MIN (NCOL, NUMCOL)
         ICOL1 = ICOL1 + NCOL

         GOTO 10
      END IF

      IXROW(NEROW+1) = ICOL1

      RETURN
      END
