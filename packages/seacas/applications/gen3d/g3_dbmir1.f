C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBMIR1 (IELB, NUMELB, NUMLNK, LINK)
C=======================================================================

C   --*** DBMIR1 *** (GEN3D) Fixup element connectivity for reflections
C   --   Written by Greg Sjaardema - revised 02/10/89
C   --   Modified from DBOEB1 Written by Amy Gilkey
C   --
C   --Parameters:
C   --   IELB - IN - the element block number
C   --   NUMELB - IN - the number of elements in the block
C   --   NUMLNK - IN - the number of nodes per element
C   --   LINK - IN/OUT - the element connectivity for this block
C   --

      INTEGER NUMELB, NUMLNK
      INTEGER LINK(NUMLNK,NUMELB)

      IF ((NUMLNK .EQ. 8) .AND. (NUMELB .GT. 0)) THEN
         DO NE = 1, NUMELB
            ILTMP = LINK (2,NE)
            LINK(2,NE) = LINK(4,NE)
            LINK(4,NE) = ILTMP

            ILTMP = LINK (6,NE)
            LINK(6,NE) = LINK(8,NE)
            LINK(8,NE) = ILTMP
          end do
        ELSE IF ((NUMLNK .EQ. 6) .AND. (NUMELB .GT. 0)) THEN
C ... UNTESTED
         DO NE = 1, NUMELB
            ILTMP = LINK (2,NE)
            LINK(2,NE) = LINK(3,NE)
            LINK(3,NE) = ILTMP

            ILTMP = LINK (5,NE)
            LINK(5,NE) = LINK(6,NE)
            LINK(6,NE) = ILTMP
          end do
      END IF

      RETURN
      END
