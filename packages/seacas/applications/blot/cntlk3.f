C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CNTLK3 (NELBLK, LENE, NLNKE, LENLNK, NFACES, NAMELB,
     $     NLNKSC)
C=======================================================================

C   --*** CNTLK3 *** (MESH) Return link array length and number faces (3D)
C   --   Written by Amy Gilkey - revised 03/04/88
C   --
C   --CNTLK3 returns the maximum length of the connectivity array for
C   --the faces and the maximum number of faces that may be defined.
C   --
C   --Parameters:
C   --   NELBLK - IN - the number of element blocks
C   --   LENE - IN - the cumulative element counts by element block
C   --   NLNKE - IN - the number of nodes per element
C   --   LENLNK - OUT - the length of the connectivity array
C   --   NELEMS - OUT - the number of elements

      INTEGER LENE(0:*)
      INTEGER NLNKE(*)
      CHARACTER*(*) NAMELB(*)
      INTEGER NLNKSC(*)

      NFACES = 0
      LENLNK = 0
      DO 100 IELB = 1, NELBLK
         NUME = LENE(IELB) - LENE(IELB-1)
         IF (NAMELB(IELB)(:3) .EQ. 'HEX' .and. NLNKE(IELB) .GE. 8) THEN
            NF = 6
            NL = 4
         ELSE IF (NAMELB(IELB)(:3) .EQ. 'TET') THEN
            NF = 4

C ... Even though there are only 3 nodes per face, we fake that there
C     are 4 nodes for consistency with the rest of blot...
            NL = 4
         ELSE IF (NAMELB(IELB)(:3) .EQ. 'WED') THEN
            NF = 5
            NL = 4
         ELSE IF (NAMELB(IELB)(:3) .EQ. 'PYR') THEN
            NF = 5
            NL = 4
         ELSE
            NF = 1
            NL = NLNKE(IELB)
         END IF
         NLNKSC(IELB) = NL
         NFACES = NFACES + NUME * NF
         LENLNK = LENLNK + NUME * NF * NL
  100 CONTINUE

      RETURN
      END
