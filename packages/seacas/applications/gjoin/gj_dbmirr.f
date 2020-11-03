C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBMIRR (NELBS, NELBE, IDELB, NUMELB, NUMLNK, LINK,
     *  NAMELB, NDIM, NONQUD)
C=======================================================================

C   --*** DBMIRR *** (GEN3D) Fixup connectivity for reflections
C   --   Written by Greg Sjaardema - revised 02/10/89
C   --   Modified from DBOELB Written by Amy Gilkey
C   --
C   --Parameters:
C   --   NELBS, NELBE - IN - the number of first and last element blocks
C   --   IDELB - IN - the element block IDs for each block
C   --   NUMELB - IN - the number of elements in each block
C   --   NUMLNK - IN - the number of nodes per element in each block
C   --   LINK - IN/OUT - the connectivity for each block
C   --
C   --   NONQUAD - IN - true if model contains non-hex/non-quad elements

      include 'exodusII.inc'
      INTEGER IDELB(*)
      INTEGER NUMELB(*)
      INTEGER NUMLNK(*)
      INTEGER LINK(*)
      CHARACTER*(MXSTLN) NAMELB(*)
      LOGICAL NONQUD

      NONQUD = .FALSE.
      IELNK = 0

      DO 100 NELB = NELBS, NELBE
         IELB = NELB-NELBS+1

         ISLNK = IELNK + 1
         IELNK = IELNK + NUMLNK(IELB) * NUMELB(IELB)

         CALL DBMIR1 (IDELB(IELB), NUMELB(IELB),
     *     NUMLNK(IELB), LINK(ISLNK), NAMELB(IELB), NDIM, NONQUD)
  100 CONTINUE

      RETURN
      END
