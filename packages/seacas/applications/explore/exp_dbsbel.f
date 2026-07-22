C     Copyright(C) 1999-2020, 2023, 2024 National Technology & Engineering Solutions
C     of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C     NTESS, the U.S. Government retains certain rights in this software.
C
C     See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE DBSBEL (NELBLK, NUMEL, LENE, INEL, NLISEL, LISEL, ADD)
C=======================================================================

C     --*** DBSBEL *** (BLOT) Select elements given list of elements
C     --
C     --DBSBEL creates the element block selection array and the element
C     --selection array (by block) given a list of selected elements.
C     --
C     --Parameters:
C     --   NELBLK - IN - the number of element blocks
C     --   NUMEL - IN - the number of elements
C     --   LENE - IN - the cumulative element counts by element block
C     --   INEL - IN - the indices of the selected elements
C     --   NLISEL - IN/OUT - the number of selected elements for each block

      INTEGER LENE(0:NELBLK)
      INTEGER INEL(0:*)
      INTEGER NLISEL(0:NELBLK)
      INTEGER LISEL(0:*)
      LOGICAL ADD

      if (.NOT. ADD) THEN
         do ielb = 1, nelblk
            nlisel(ielb) = 0
         end do
         nlisel(0) = 0
         do i = 1, numel
            lisel(i) = 0
         end do
      end if

      DO IELB = 1, NELBLK
         IEL = LENE(IELB-1)+1
         LEL = LENE(IELB)
         NEL = NLISEL(IELB)
         DO I = 1, INEL(0)
            IF ((INEL(I) .GE. IEL) .AND. (INEL(I) .LE. LEL)) THEN
               ISEL = INEL(I)
               IF (lisel(isel) .eq. 0) then
                  NEL = NEL + 1
                  lisel(isel) = isel
               end if
            END IF
         end do
         NLISEL(IELB) = NEL
         IF (NEL .GT. 0) NLISEL(0) = NLISEL(0) + 1
      end do

      END
