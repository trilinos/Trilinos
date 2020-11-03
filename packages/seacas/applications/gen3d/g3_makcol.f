C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MAKCOL (LINK, NUMNP, NUMCOL, NEROW, IELCEN, ICOL1,
     &   IRBOT, IRTOP, IRBDIF,
     &   NPBROW, NPBCOL, NPTROW, NPTCOL)
C=======================================================================

C   --*** MAKCOL *** (GEN3D) Compute center block row x column
C   --   Written by Amy Gilkey - revised 04/26/88
C   --
C   --MAKCOL orders the elements in a center block by row and column.
C   --The elements have already been ordered into rows, so this routine
C   --needs to determine which column the first element of a row is in.
C   --It does this by comparing the rows on top and bottom of a given row.
C   --
C   --Parameters:
C   --   LINK - IN - the connectivity for the 2D elements, always 4 nodes
C   --   NUMNP - IN - the number of nodes
C   --   NUMCOL - IN - the number of columns in the center blocks
C   --   NEROW - IN - the number of element rows in the center blocks
C   --   IELCEN - IN - the element numbers of the center elements
C   --      by column and row (column 1 is not necessarily the center)
C   --   ICOL1 - IN/OUT - the column number of the first element in a row;
C   --      upon entry, 1 iff center row, else 0; destroyed
C   --   IRBOT, IRTOP - OUT - the row number of the row on top and bottom
C   --   IRBDIF - SCRATCH - size = NEROW
C   --   NPBROW, NPBCOL, NPTROW, NPTCOL - SCRATCH - size = NUMNP
C   --
C   --Common Variables:
C   --   Uses IX1, IX2, IX3, IX4 of /CENPAR/

      INCLUDE 'g3_cenpar.blk'

      INTEGER LINK(4,*)
      INTEGER IELCEN(NUMCOL,*)
      INTEGER ICOL1(*)
      INTEGER IRBOT(*), IRTOP(*)
      INTEGER IRBDIF(*)
      INTEGER NPBROW(*), NPBCOL(*), NPTROW(*), NPTCOL(*)

      LOGICAL PROBLM

      PROBLM = .FALSE.

C   --Find out the corresponding rows and columns for each node

      CALL INIINT (NUMNP, 0, NPBROW)
      CALL INIINT (NUMNP, 0, NPTROW)

      DO 20 IROW = 1, NEROW
         DO 10 ICOL = 1, NUMCOL
            IEL = IELCEN(ICOL,IROW)
            IF (IEL .LE. 0) GOTO 10
            INP = LINK(IX1,IEL)
            NPTROW(INP) = IROW
            NPTCOL(INP) = ICOL
            INP = LINK(IX4,IEL)
            NPBROW(INP) = IROW
            NPBCOL(INP) = ICOL
   10    CONTINUE
   20 CONTINUE

C   --Find the row on top and on bottom of each row and check column differences

      CALL INIINT (NEROW, 0, IRBOT)
      CALL INIINT (NEROW, 0, IRTOP)
      CALL INIINT (NEROW, 0, IRBDIF)

      DO 40 IROW = 1, NEROW
         DO 30 ICOL = 1, NUMCOL
            IEL = IELCEN(ICOL,IROW)
            IF (IEL .LE. 0) GOTO 30
            INP = LINK(IX1,IEL)

            ITROW = NPTROW(INP)

            IF (NPBROW(INP) .NE. 0) THEN
C            --Mark the bottom row and point back to this row
               IBROW = NPBROW(INP)
               IRBOT(ITROW) = IBROW
               IRTOP(IBROW) = ITROW
C            --Mark column differences
               IRBDIF(ITROW) = NPTCOL(INP) - NPBCOL(INP)
               GOTO 40
            END IF
   30    CONTINUE
   40 CONTINUE

C   --Check columns of bottom and top row

      DO 60 IROW = 1, NEROW

         DO 50 ICOL = 1, NUMCOL
            IEL = IELCEN(ICOL,IROW)
            IF (IEL .LE. 0) GOTO 50
            IB = LINK(IX1,IEL)
            IT = LINK(IX4,IEL)

C         --Check that node is not connected to two top or bottom rows
            IF (NPTROW(IB) .NE. IROW) PROBLM = .TRUE.
            IF (NPBROW(IT) .NE. IROW) PROBLM = .TRUE.

            IF (NPBROW(IB) .NE. 0) THEN
C            --Check that same bottom row this node
               IF (NPBROW(IB) .NE. IRBOT(IROW)) PROBLM = .TRUE.
C            --Check column difference
               IF ((NPTCOL(IB) - NPBCOL(IB)) .NE. IRBDIF(IROW))
     &            PROBLM = .TRUE.
            ELSE IF (NPTROW(IT) .NE. 0) THEN
C            --Check that same top row this node
               IF (NPTROW(IT) .NE. IRTOP(IROW)) PROBLM = .TRUE.
C            --Check column difference
               ITROW = NPTROW(IT)
               IF ((NPTCOL(IT) - NPBCOL(IT)) .NE. IRBDIF(ITROW))
     &            PROBLM = .TRUE.
            END IF

   50    CONTINUE
   60 CONTINUE

C   --Get the column numbers of non-center rows by following a chain of rows
C   --through the top rows

      DO 90 IROW = 1, NEROW
C      --Start at a center row
         IF (ICOL1(IROW) .NE. 1) GOTO 90

         NROW = IROW
   70    CONTINUE
C      --Check the row if it has top-bottom row pairing
         IF (IRBOT(NROW) .GT. 0) THEN
            IBROW = IRBOT(NROW)

C         --Calculate the column from the paired row column and offset
            ICOL = ICOL1(NROW) + IRBDIF(NROW)

            IF ((ICOL1(IBROW) .GT. 1) .OR. (ICOL .LE. 0)) THEN
C            --If the paired row is already marked, stop chain
               IF (ICOL .LE. 0) PROBLM = .TRUE.
               IF (ICOL1(IBROW) .NE. ICOL) PROBLM = .TRUE.
               ICOL = 1
            END IF

            IF ((ICOL .GT. 1) .AND. (ICOL .LE. NUMCOL)) THEN
C            --Chain to the paired row
               ICOL1(IBROW) = ICOL
               NROW = IBROW

               GOTO 70
            END IF
         END IF

         NROW = IROW
   80    CONTINUE
C      --Check the row if it has top-bottom row pairing
         IF (IRTOP(NROW) .GT. 0) THEN
            ITROW = IRTOP(NROW)

C         --Calculate the column from the paired row column and offset
            ICOL = ICOL1(NROW) - IRBDIF(ITROW)

            IF ((ICOL1(ITROW) .GE. 1) .OR. (ICOL .LE. 0)) THEN
C            --If the paired row is already marked, stop chain
               IF (ICOL .LE. 0) PROBLM = .TRUE.
               IF (ICOL1(ITROW) .NE. ICOL) PROBLM = .TRUE.
               ICOL = 1
            END IF

            IF ((ICOL .GT. 1) .AND. (ICOL .LE. NUMCOL)) THEN
C            --Chain to the paired row
               ICOL1(ITROW) = ICOL
               NROW = ITROW

               GOTO 80
            END IF
         END IF
   90 CONTINUE

C   --Shift rows over if not on first column

      NROW = NEROW
      NEROW = 0
      DO 120 IROW = 1, NROW
C      --Eliminate row if its first column is too big
         IF ((ICOL1(IROW) .LT. 1) .OR. (ICOL1(IROW) .GT. NUMCOL)) THEN
            ICOL1(IROW) = 0

         ELSE
            NEROW = NEROW + 1

            IC = NUMCOL - ICOL1(IROW) + 1
            DO 100 ICOL = NUMCOL, ICOL1(IROW), -1
               IELCEN(ICOL,NEROW) = IELCEN(IC,IROW)
               IC = IC - 1
  100       CONTINUE
            DO 110 ICOL = 1, ICOL1(IROW)-1
               IELCEN(ICOL,NEROW) = 0
  110       CONTINUE

            ICOL1(IROW) = NEROW
         END IF
  120 CONTINUE

C   --Connect the row top/bottom pointers to the renumbered rows

      DO 130 IROW = 1, NEROW
         IR = IRTOP(IROW)
         IF (IR .GT. 0) IRTOP(IROW) = ICOL1(IR)
         IR = IRBOT(IROW)
         IF (IR .GT. 0) IRBOT(IROW) = ICOL1(IR)
  130 CONTINUE

      RETURN
      END
