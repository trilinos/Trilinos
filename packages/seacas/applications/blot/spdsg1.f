C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SPDSG1 (NNENUM, NENUM, IE2ELB, ISEVOK, NEL, ISEGEL)
C=======================================================================

C   --*** SPDSG1 *** (SPLOT) Find all defined element values
C   --   Written by Amy Gilkey - revised 11/05/87
C   --
C   --SPDSG1 returns a list of the defined element indices for an
C   --element variable.
C   --
C   --Parameters:
C   --   NNENUM - IN - the number of selected nodes/elements
C   --   NENUM - IN - the selected nodes/elements
C   --   IE2ELB - IN - the element block for each element
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable of block j exists iff ISEVOK(j)
C   --   NEL - OUT - the number of defined elements
C   --   ISEGEL - OUT - the NENUM indices of the defined elements;
C   --      ISEGEL(0) = the number of defined elements

      INTEGER NENUM(*)
      INTEGER IE2ELB(*)
      LOGICAL ISEVOK(*)
      INTEGER ISEGEL(0:*)

      NEL = 0
      DO 100 I = 1, NNENUM
         IF (ISEVOK(IE2ELB(NENUM(I)))) THEN
            NEL = NEL + 1
            ISEGEL(NEL) = I
         END IF
  100 CONTINUE
      ISEGEL(0) = NEL

      RETURN
      END
