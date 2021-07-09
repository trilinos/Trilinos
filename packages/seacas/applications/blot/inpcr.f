C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE INPCR (NUMNPF, HIDENP, XN, YN, ZN,
     &   COLMIN, ROWMIN, CRDELT, NUMCOL, NUMROW,
     &   IXNCRS, IXNCRE, NPCR, LNPCR, NPX, NPY)
C=======================================================================

C   --*** INPCR *** (MESH) Initialize the NPCR array
C   --   Written by Amy Gilkey - revised 04/04/88
C   --
C   --INPCR initializes the IXNCRS, IXNCRE, and NPCR arrays to point to the
C   --set of all visible nodes in a each block.
C   --
C   --Parameters:
C   --   NUMNPF - IN - the number of nodes
C   --   HIDENP - IN/OUT - node status (as in HIDDEN)
C   --   XN, YN, ZN - IN - the nodal coordinates
C   --   COLMIN, ROWMIN - IN - the minimum column and row value and 1 / interval
C   --   CRDELT - the column and row interval reciprical (1 / interval)
C   --   NUMCOL, NUMROW - IN - the number of columns and rows
C   --   IXNCRS, IXNCRE - OUT - the starting and ending NPCR index
C   --      for each column and row
C   --   NPCR - OUT - the nodes indexed by IXNCRS, IXNCRE
C   --   LNPCR - OUT - the length of the NPCR array
C   --   NPX, NPY - SCRATCH - size = NUMNPF

      PARAMETER (KNVIS=0, KNFOVR=10, KNHID=100)

      INTEGER HIDENP(NUMNPF)
      REAL XN(NUMNPF), YN(NUMNPF), ZN(NUMNPF)
      INTEGER IXNCRS(0:NUMCOL,0:NUMROW), IXNCRE(0:NUMCOL,0:NUMROW)
      INTEGER NPCR(NUMNPF)
      INTEGER NPX(NUMNPF), NPY(NUMNPF)

      CALL INIINT ((1+NUMCOL) * (1+NUMROW), 0, IXNCRS)

      DO 100 INP = 1, NUMNPF
         IF (HIDENP(INP) .LE. KNVIS) THEN
            ICOL = INT ((XN(INP) - COLMIN) * CRDELT)
            IROW = INT ((YN(INP) - ROWMIN) * CRDELT)
            NPX(INP) = ICOL
            NPY(INP) = IROW
            IXNCRS(ICOL,IROW) = IXNCRS(ICOL,IROW) + 1
         END IF
  100 CONTINUE

      NHEAP = 1
      DO 120 IROW = 0, NUMROW
         DO 110 ICOL = 0, NUMCOL
            L = IXNCRS(ICOL,IROW)
            IXNCRS(ICOL,IROW) = NHEAP
            IXNCRE(ICOL,IROW) = NHEAP - 1
            NHEAP = NHEAP + L
  110    CONTINUE
  120 CONTINUE

      LNPCR = NHEAP - 1

      DO 130 INP = 1, NUMNPF
         IF (HIDENP(INP) .LE. KNVIS) THEN
            ICOL = NPX(INP)
            IROW = NPY(INP)
            NHEAP = IXNCRE(ICOL,IROW) + 1
            IXNCRE(ICOL,IROW) = NHEAP
            NPCR(NHEAP) = INP
         END IF
  130 CONTINUE

      RETURN
      END
