C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE ALLCUT (IE2ELB, LENF, IF2EL, IF2EL2, NEWELB)
C=======================================================================

C   --*** ALLCUT *** (MESH) Uncut 3D mesh
C   --   Written by Amy Gilkey - revised 03/08/88
C   --
C   --ALLCUT restores a cut 3D mesh to a whole mesh.
C   --
C   --Parameters:
C   --   IE2ELB - IN - the element block for each element
C   --   LENF - IN - the cumulative face counts by element block
C   --   IF2EL - IN - the element number of each face
C   --   IF2EL2 - IN - the secondary element number of each face
C   --   NEWELB - OUT - size = LENF(NELBLK+1)
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/

      include 'dbnums.blk'
      include 'd3nums.blk'

      INTEGER IE2ELB(NUMEL)
      INTEGER LENF(0:NELBLK+3)
      INTEGER IF2EL(*), IF2EL2(*)
      INTEGER NEWELB(*)

C   --Move the interior surface faces back to the interior set

      DO 110 IELB = 1, NELBLK+3
         DO 100 IFAC = LENF(IELB-1)+1, LENF(IELB)
            IF (IF2EL2(IFAC) .LE. 0) THEN
               NEWELB(IFAC) = IELB
            ELSE
               NEWELB(IFAC) = NELBLK+1
            END IF
  100    CONTINUE
  110 CONTINUE

C   --Move the faces that were OUT to their original set

      IELB = NELBLK+3
      DO 120 IFAC = LENF(IELB-1)+1, LENF(IELB)
         IF (IF2EL2(IFAC) .LE. 0) THEN
            NEWELB(IFAC) = IE2ELB(IF2EL(IFAC))
         ELSE
            NEWELB(IFAC) = NELBLK+1
         END IF
  120 CONTINUE

      RETURN
      END
