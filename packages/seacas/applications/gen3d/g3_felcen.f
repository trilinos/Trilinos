C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE FELCEN (NUMCOL, NEROW, IELROW, IXROW, IELCEN)
C=======================================================================

C   --*** FELCEN *** (GEN3D) Put center elements into row by column array
C   --   Written by Amy Gilkey - revised 04/26/88
C   --
C   --It puts the center elements into a row by column array.
C   --
C   --Parameters:
C   --   NUMCOL - IN - the number of columns in the center blocks
C   --   NEROW - IN - the number of element rows in the center blocks
C   --   IELROW - IN - the element numbers of the rows of center elements
C   --   IXROW - IN - the IELROW index of the starting column for each row
C   --   IELCEN - OUT - the element numbers of the center elements
C   --      by column and row (column 1 is not necessarily the center)

      INTEGER IELROW(*)
      INTEGER IXROW(NEROW+1)
C...Assert NUMCOL > 0
      INTEGER IELCEN(NUMCOL,*)

C   --Put rows into row by column array

      DO 20 IROW = 1, NEROW
         IX = IXROW(IROW)
         DO 10 ICOL = 1, NUMCOL
            IF (IX .LE. IXROW(IROW+1)-1) THEN
               IELCEN(ICOL,IROW) = IELROW(IX)
               IX = IX + 1
            ELSE
               IELCEN(ICOL,IROW) = 0
            END IF
   10    CONTINUE
   20 CONTINUE

      RETURN
      END
