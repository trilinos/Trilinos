C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MAKN2B (LENF, NLNKF, LINKF, IELBST, IN2ELB)
C=======================================================================

C   --*** MAKN2B *** (MESH) Create node to element block index
C   --   Written by Amy Gilkey - revised 03/09/88
C   --
C   --MAKN2B finds the element block to be associated with each node.  The
C   --element block is determined by examining the faces which contain the
C   --node.  If they all are of the same element block, this is the node's
C   --element block.  If not, the node is assigned element block 0.
C   --Only selected element blocks are included in the calculation.
C   --
C   --Parameters:
C   --   LENF - IN - the cumulative face counts by element block
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF - IN - the connectivity for all faces
C   --   IELBST - IN - the element block status (>0 if selected)
C   --   IN2ELB - OUT - the element block for each node;
C   --      <0 if not in any selected element block
C   --      =0 if in more than one selected element block
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/
C   --   Uses IS3DIM, NUMNPF of /D3NUMS/

      include 'dbnums.blk'
      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM

      INTEGER LENF(0:NELBLK+2)
      INTEGER NLNKF(NELBLK)
      INTEGER LINKF(*)
      INTEGER IELBST(NELBLK)
      INTEGER IN2ELB(NUMNPF)

      CALL INIINT (NUMNPF, -999, IN2ELB)

      DO 120 IELB = 1, NELBLK
         IF (IELBST(IELB) .GT. 0) THEN
            IXL0 = IDBLNK (IELB, 0, LENF, NLNKF) - 1
            DO 110 IFAC = LENF(IELB-1)+1, LENF(IELB)
               DO 100 K = 1, NLNKF(IELB)
                  INP = LINKF(IXL0+K)
                  IF (IN2ELB(INP) .LT. 0) THEN
                     IN2ELB(INP) = IELB
                  ELSE IF (IN2ELB(INP) .NE. IELB) THEN
                     IN2ELB(INP) = 0
                  END IF
  100          CONTINUE
               IXL0 = IXL0 + NLNKF(IELB)
  110       CONTINUE
         END IF
  120 CONTINUE

      RETURN
      END
