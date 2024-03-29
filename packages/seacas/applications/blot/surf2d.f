C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SURF2D (LENE, NLNKE, LINKE,
     &   LENF, NLNKF, LINKF, IF2EL, NFACES, LENLNK)
C=======================================================================

C   --*** SURF2D *** (MESH) Group 2D faces
C   --   Written by Amy Gilkey - revised 11/04/87
C   --
C   --SURF2D packs the elements by element block, eliminating empty elements.
C   --
C   --Parameters:
C   --   LENE - IN - the cumulative element counts by element block
C   --   NLNKE - IN - the number of nodes per element
C   --   LINKE - IN - the original connectivity; connectivity all zero
C   --      if element is undefined
C   --   LENF - OUT - the cumulative face counts by element block
C   --   NLNKF - OUT - the number of nodes per face
C   --   LINKF - OUT - the connectivity for all faces
C   --   IF2EL - OUT - the element number of each face
C   --   NFACES - OUT - the number of returned faces
C   --   LENLNK - OUT - the length of the returned LINKF array
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/

      common /debugc/ cdebug
      common /debugn/ idebug
      character*8 cdebug

      include 'dbnums.blk'
      include 'd3nums.blk'

      INTEGER LENE(0:NELBLK), LINKE(*)
      INTEGER NLNKE(NELBLK)
      INTEGER LENF(0:NELBLK+2)
      INTEGER NLNKF(NELBLK)
      INTEGER LINKF(*)
      INTEGER IF2EL(*)

      LENF(0) = 0
      NFACES = 0
      IXL0 = 0
      DO 110 IELB = 1, NELBLK
         NLNKF(IELB) = NLNKE(IELB)
         IXE0 = IDBLNK (IELB, 0, LENE, NLNKE) - 1
         DO 100 IEL = LENE(IELB-1)+1, LENE(IELB)
            IF (LINKE(IXE0+1) .NE. 0) THEN
               NFACES = NFACES + 1
               CALL CPYINT (NLNKF(IELB), LINKE(IXE0+1), LINKF(IXL0+1))
               IF2EL(NFACES) = IEL
               IXL0 = IXL0 + NLNKF(IELB)
            END IF
            IXE0 = IXE0 + NLNKE(IELB)
  100    CONTINUE
         LENF(IELB) = NFACES
  110 CONTINUE

      LENF(NELBLK+1) = NFACES
      LENF(NELBLK+2) = NFACES

      LENLNK = IXL0
      RETURN
      END
