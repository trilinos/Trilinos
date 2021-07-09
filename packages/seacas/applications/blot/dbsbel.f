C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBSBEL (NELBLK, NUMEL, LENE, INEL, NLISEL, LISEL)
C=======================================================================

C   --*** DBSBEL *** (BLOT) Select elements given list of elements
C   --   Written by Amy Gilkey - revised 01/05/88
C   --
C   --DBSBEL creates the element block selection array and the element
C   --selection array (by block) given a list of selected elements.
C   --
C   --Parameters:
C   --   NELBLK - IN - the number of element blocks
C   --   NUMEL - IN - the number of elements
C   --   LENE - IN - the cumulative element counts by element block
C   --   INEL - IN - the indices of the selected elements
C   --   NLISEL - IN/OUT - the number of selected elements for each block
C   --   LISEL - IN/OUT - the indices of the selected elements (by block)

      INTEGER LENE(0:*)
      INTEGER INEL(0:*)
      INTEGER NLISEL(0:*)
      INTEGER LISEL(0:*)

      LISEL(0) = 0
      CALL INIINT (NUMEL, 0, LISEL(1))

      NLISEL(0) = 0
      DO 110 IELB = 1, NELBLK
         IEL = LENE(IELB-1)+1
         LEL = LENE(IELB)
         NEL = 0
         DO 100 I = 1, INEL(0)
            IF ((INEL(I) .GE. IEL) .AND. (INEL(I) .LE. LEL)) THEN
               LISEL(IEL+NEL) = INEL(I)
               NEL = NEL + 1
            END IF
  100    CONTINUE
         LISEL(0) = LISEL(0) + NEL
         NLISEL(IELB) = NEL
         IF (NEL .GT. 0) NLISEL(0) = NLISEL(0) + 1
  110 CONTINUE

      RETURN
      END
