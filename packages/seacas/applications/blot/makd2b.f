C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MAKD2B (LENF, NLNKF, LINKF, IELBST,
     &   IF2EL, IF2EL2, IE2ELB, IDN2B)
C=======================================================================

C   --*** MAKD2B *** (MESH) Create dead node to element block index
C   --   Written by Amy Gilkey - revised 03/09/88
C   --
C   --MAKD2B finds the element block to be associated with each dead node.
C   --A dead node is a node which is in a dead surface face (if 3D),
C   --but not in any alive faces.  Only selected element blocks are
C   --included in the calculation.
C   --
C   --Parameters:
C   --   LENF - IN - the cumulative face counts by element block
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF - IN - the connectivity for all faces
C   --   IELBST - IN - the element block status (>0 if selected)
C   --   IF2EL - IN - the element number of each face
C   --   IF2EL2 - IN - the secondary element number of each face
C   --   IE2ELB - IN - the element block for each element
C   --   IDN2B - OUT - the element block for each dead node; dead if >= 0
C   --
C   --Common Variables:
C   --   Uses NDIM, NELBLK of /DBNUMS/
C   --   Uses IS3DIM, NUMNPF of /D3NUMS/

      include 'dbnums.blk'
      include 'd3nums.blk'

      INTEGER LENF(0:NELBLK+2)
      INTEGER NLNKF(NELBLK)
      INTEGER LINKF(*)
      INTEGER IELBST(NELBLK)
      INTEGER IF2EL(*), IF2EL2(*)
      INTEGER IE2ELB(NUMEL)
      INTEGER IDN2B(NUMNPF)

C   --Mark all nodes with -999

      CALL INIINT (NUMNPF, -999, IDN2B)

C   --Calculate starting link index for dead faces

      IELB = NELBLK+1
      IXL0 = IDBLNK (IELB, 0, LENF, NLNKF) - 1
      DO 100 IFAC = LENF(IELB-1)+1, LENF(IELB)
         IB = IE2ELB(IF2EL(IFAC))
         IXL0 = IXL0 + NLNKF(IB)
  100 CONTINUE
      IELB = NELBLK+2

      IF (.NOT. IS3DIM) THEN

C      --Mark nodes in dead faces if not in any alive faces

         DO 120 IFAC = LENF(IELB-1)+1, LENF(IELB)
            IB = IE2ELB(IF2EL(IFAC))
            IF (IELBST(IB) .GT. 0) THEN
               DO 110 K = 1, NLNKF(IB)
                  INP = LINKF(IXL0+K)
                  IF (IDN2B(INP) .LT. 0) THEN
                     IDN2B(INP) = IB
                  ELSE IF (IDN2B(INP) .NE. IB) THEN
                     IDN2B(INP) = 0
                  END IF
  110          CONTINUE
            END IF
            IXL0 = IXL0 + NLNKF(IB)
  120    CONTINUE

      ELSE

C      --Mark nodes in dead surface faces if not in any alive faces

         DO 140 IFAC = LENF(IELB-1)+1, LENF(IELB)
            IB = IE2ELB(IF2EL(IFAC))
            IF (IELBST(IB) .GT. 0) THEN
               IF (IF2EL2(IFAC) .LE. 0) THEN
                  DO 130 K = 1, NLNKF(IB)
                     INP = LINKF(IXL0+K)
                     IF (IDN2B(INP) .LT. 0) THEN
                        IDN2B(INP) = IB
                     ELSE IF (IDN2B(INP) .NE. IB) THEN
                        IDN2B(INP) = 0
                     END IF
  130             CONTINUE
               END IF
            END IF
            IXL0 = IXL0 + NLNKF(IB)
  140    CONTINUE
      END IF

C   --Mark alive nodes (any element block) with -999

      DO 170 IELB = 1, NELBLK
         IXL0 = IDBLNK (IELB, 0, LENF, NLNKF) - 1
         DO 160 IFAC = LENF(IELB-1)+1, LENF(IELB)
            DO 150 K = 1, NLNKF(IELB)
               INP = LINKF(IXL0+K)
               IDN2B(INP) = -999
  150       CONTINUE
            IXL0 = IXL0 + NLNKF(IELB)
  160    CONTINUE
  170 CONTINUE

      RETURN
      END
