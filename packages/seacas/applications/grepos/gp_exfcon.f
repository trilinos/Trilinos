C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE EXFCON (IBLK, NUMELB, NUMLNK, LINK, ICONOD)
C=======================================================================

C   --*** DBMIR1 *** (GEN3D) Fixup element connectivity for reflections
C   --   Written by Greg Sjaardema - revised 02/10/89
C   --   Modified from DBOEB1 Written by Amy Gilkey
C   --
C   --Parameters:
C   --   NUMELB - IN - the number of elements in the block
C   --   NUMLNK - IN - the number of nodes per element
C   --   LINK   - IN - the element connectivity for this block
C   --

      INTEGER NUMELB, NUMLNK
      INTEGER LINK(NUMLNK,*)
      INTEGER ICONOD(*)

      DO 20 NLNK = 1, NUMLNK
         DO 10 NE = 1, NUMELB
            INODE = LINK (NLNK, NE)
            ICONOD(INODE) = IBLK
   10    CONTINUE
   20 CONTINUE
      RETURN
      END
