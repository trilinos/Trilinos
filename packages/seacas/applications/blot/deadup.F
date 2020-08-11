C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DEADUP (HIDENP, XN, YN, ZN, IDN2B, *)
C=======================================================================

C   --*** DEADUP *** (MESH) Plot "dead" nodes with symbols
C   --   Written by Amy Gilkey - revised 03/09/88
C   --
C   --DEADUP draws a symbol at each dead node.  A dead node is a node
C   --which is not attached to any alive elements.
C   --Only nodes of selected element blocks are drawn.
C   --
C   --Parameters:
C   --   HIDENP(i) - IN - true iff node i is hidden (3D only)
C   --   XN, YN, ZN - IN - the nodal coordinates
C   --   IDN2B - IN - the element block for each dead node; dead if >= 0;
C   --      =0 if in more than one selected element block
C   --   * - return statement if the cancel function is active
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/
C   --   Uses IS3DIM, NUMNPF of /D3NUMS/

      include 'dbnums.blk'
      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM

      LOGICAL HIDENP(NUMNPF)
      REAL XN(NUMNPF), YN(NUMNPF), ZN(NUMNPF)
      INTEGER IDN2B(NUMNPF)

      LOGICAL GRABRT

      DO 110 IELB = 0, NELBLK

C      --Set the node element block color
         CALL GRCOLR (IELB)

         DO 100 INP = 1, NUMNPF

            IF (IS3DIM) THEN
               IF (HIDENP(INP)) GOTO 100
            END IF

C         --If node is dead and of the appropriate element block, draw symbol

            IF (IDN2B(INP) .EQ. IELB) THEN
               IF (GRABRT ()) RETURN 1
#if NeedsDoubleEscape
               CALL MPD2SY (1, XN(INP), YN(INP), '\\CDO')
#else
               CALL MPD2SY (1, XN(INP), YN(INP), '\CDO')
#endif
            END IF

  100    CONTINUE

         CALL PLTFLU
  110 CONTINUE

      RETURN
      END
