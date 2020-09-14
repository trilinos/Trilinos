C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MAKSU2 (LENL, LINSET, MAXMSH, DODEAD, IDN2B, NPSURF)
C=======================================================================

C   --*** MAKSU2 *** (MESH) Determine which nodes are on surface (2D)
C   --   Written by Amy Gilkey - revised 03/10/88
C   --
C   --MAKSU2 calculates the nodes-on-surface indicator (NPSURF).
C   --
C   --Parameters:
C   --   LENL - IN/OUT - the cumulative line counts by element block
C   --   LINSET - IN/OUT - the sorted line set
C   --   MAXMSH - IN - the maximum mesh line types to be displayed
C   --      (as in MSHLIN of /MSHOPT/)
C   --   DODEAD - IN - add in dead nodes from IDN2B iff true
C   --   IDN2B - IN - the element block for each dead node; dead if >= 0
C   --   NPSURF - OUT - the node numbers of the mesh boundary nodes
C   --
C   --Common Variables:
C   --   Uses NDIM, NELBLK of /DBNUMS/
C   --   Sets NNPSUR, NUMNPF of /D3NUMS/

      PARAMETER (MSHNON=0, MSHBOR=1, MSHDIV=2, MSHSEL=3, MSHALL=4)

      include 'dbnums.blk'
      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM

      INTEGER LENL(-2:NELBLK), LINSET(LLNSET,*)
      LOGICAL DODEAD
      INTEGER IDN2B(NUMNPF)
      INTEGER NPSURF(NUMNPF)

C   --Determine nodes to be included

      IF (MAXMSH .LE. MSHBOR) THEN
         NODSUR = -1
      ELSE IF (MAXMSH .LE. MSHDIV) THEN
         NODSUR = 0
      ELSE IF (MAXMSH .LE. MSHALL) THEN
         NODSUR = NELBLK
      END IF

      CALL INIINT (NUMNPF, 0, NPSURF)

C   --Mark all nodes in the selected line sets

      DO 110 IL = 1, LENL(NODSUR)
         DO 100 K = 1, 2
            NPSURF(LINSET(K,IL)) = 1
  100    CONTINUE
  110 CONTINUE

C   --Add in dead nodes

      IF (DODEAD) THEN
         DO 120 INP = 1, NUMNPF
            IF (IDN2B(INP) .GE. 0) NPSURF(INP) = 1
  120    CONTINUE
      END IF

C   --Make up the list of selected nodes

      NNPSUR = 0
      DO 130 INP = 1, NUMNPF
         IF (NPSURF(INP) .GT. 0) THEN
            NNPSUR = NNPSUR + 1
            NPSURF(NNPSUR) = INP
         END IF
  130 CONTINUE

      RETURN
      END
