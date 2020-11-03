C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MAKSUR (LENF, NLNKF, LINKF, DODEAD, IDN2B, NPSURF)
C=======================================================================

C   --*** MAKSUR *** (MESH) Determine which nodes are on surface
C   --   Written by Amy Gilkey - revised 03/10/88
C   --
C   --MAKSUR calculates the nodes-on-surface indicator (NPSURF).
C   --
C   --Parameters:
C   --   LENF - IN - the cumulative face counts by element block
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF - IN - the connectivity for all faces
C   --   DODEAD - IN - add in dead nodes from IDN2B iff true
C   --   IDN2B - IN - the element block for each dead node; dead if >= 0
C   --   NPSURF - OUT - the node numbers of the surface nodes
C   --
C   --Common Variables:
C   --   Uses NDIM, NELBLK of /DBNUMS/
C   --   Sets NNPSUR, NUMNPF of /D3NUMS/

      include 'dbnums.blk'
      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM

      INTEGER LENF(0:NELBLK+2)
      INTEGER NLNKF(NELBLK)
      INTEGER LINKF(*)
      LOGICAL DODEAD
      INTEGER IDN2B(NUMNPF)
      INTEGER NPSURF(*)

      CALL INIINT (NUMNPF, 0, NPSURF)

C   --Mark all nodes that are in surface faces

      DO 120 IELB = 1, NELBLK
         IXL0 = IDBLNK (IELB, 0, LENF, NLNKF) - 1
         DO 110 IFAC = LENF(IELB-1)+1, LENF(IELB)
            DO 100 K = 1, NLNKF(IELB)
               NPSURF(LINKF(IXL0+K)) = 1
  100       CONTINUE
            IXL0 = IXL0 + NLNKF(IELB)
  110    CONTINUE
  120 CONTINUE

C   --Add in dead nodes

      IF (DODEAD) THEN
         DO 130 INP = 1, NUMNPF
            IF (IDN2B(INP) .GE. 0) NPSURF(INP) = 1
  130    CONTINUE
      END IF

C   --Make up the list of surface nodes

      NNPSUR = 0
      DO 140 INP = 1, NUMNPF
         IF (NPSURF(INP) .GT. 0) THEN
            NNPSUR = NNPSUR + 1
            NPSURF(NNPSUR) = INP
         END IF
  140 CONTINUE

      RETURN
      END
