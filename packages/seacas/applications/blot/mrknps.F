C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MRKNPS (HIDENP,
     &   XN, YN, ZN, IX2NP, IN2ELB, DODEAD, IDN2B,
     &   NNPSET, ISSNPS, IDNPS, NNNPS, IXNNPS, LTNNPS, *)
C=======================================================================

C   --*** MRKNPS *** (SETS) Mark node set on mesh
C   --   Written by Amy Gilkey - revised 04/08/88
C   --
C   --MRKNPS marks all the nodes in selected node sets on the mesh.
C   --Each node is marked with an X in the color assigned to its node set.
C   --
C   --Parameters:
C   --   HIDENP(i)  - IN - true iff node i is hidden (3D only)
C   --   XN, YN, ZN - IN - the nodal coordinates
C   --   IX2NP      - IN - the node number for each mesh index
C   --   IN2ELB     - IN - the element block for each node;
C   --                     <0 if not in any selected element block
C   --                     =0 if in more than one selected element block
C   --   DODEAD     - IN - mark dead nodes iff true
C   --   IDN2B      - IN - the element block for each dead node; dead if >= 0
C   --   NNPSET     - IN - the number of selected node sets
C   --   ISSNPS     - IN - the indices of the selected node sets
C   --   IDNPS      - IN - the node set ID for each set
C   --   NNNPS      - IN - the number of nodes for each set
C   --   IXNNPS     - IN - the index of the first node for each set
C   --   LTNNPS     - IN - the nodes for all sets
C   --   * - return statement if the cancel function is active
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/
C   --   Uses IS3DIM, NUMNPF of /D3NUMS/

      PARAMETER (KHCHSZ=1, KSCHSZ=2)

      common /debugc/ cdebug
      common /debugn/ idebug
      character*8 cdebug

      include 'dbnums.blk'
      include 'dbnumgq.blk'
      include 'd3nums.blk'

      LOGICAL HIDENP(*)
      REAL XN(*), YN(*), ZN(*), NPFACT
      INTEGER IX2NP(NUMNPF)
      INTEGER IN2ELB(NUMNPF)
      LOGICAL DODEAD
      INTEGER IDN2B(NUMNPF)
      INTEGER ISSNPS(*)
      INTEGER IDNPS(*)
      INTEGER NNNPS(*)
      INTEGER IXNNPS(*)
      INTEGER LTNNPS(*)

      LOGICAL PLTGTT, PLTSTT, LDUM
      logical logt

C   --Set the symbol size, do not exit before resetting
      LDUM = PLTGTT (KSCHSZ, SZSYM)
C     If no sets are selected, return
      IF (NNPSET .EQ. 0) THEN
         RETURN
      END IF

C     sum the number of nodes for all selected node sets
      ISUM = 0
      DO 10 I = 1, NNPSET
         ISUM = ISUM + NNNPS(ISSNPS(I))
 10   CONTINUE

C     compute the npfact
      IF (ISUM .LT. 100) THEN
         NPFACT = 1.0
      ELSE IF (ISUM .GT. 2100) THEN
         NPFACT = 0.40
      ELSE
         NPFACT = -0.0003 * ISUM + 1.03
      END IF
      LDUM = PLTSTT (KSCHSZ, NPFACT*SZSYM)

      DO 110 IX = 1, NNPSET
         INPS = ISSNPS(IX)

         CALL GRCOLR (INPS)

         IX0 = IXNNPS(INPS) - 1
         DO 100 INE = 1, NNNPS(INPS)
            INP = LTNNPS(IX0+INE)
            IF (IS3DIM) THEN
               IF (HIDENP(INP)) GOTO 100
            END IF

            logt = .false.
            if (in2elb(inp) .ge. 0) then
               logt = .true.
            else if (DODEAD) then
               if (IDN2B(inp) .GE. 0) logt = .true.
            end if
            if (logt) then
#if NeedsDoubleEscape
               CALL MPD2SY (1, XN(INP), YN(INP), '\\CCI')
#else
               CALL MPD2SY (1, XN(INP), YN(INP), '\CCI')
#endif
            END IF
  100    CONTINUE

         CALL PLTFLU
  110 CONTINUE

C   --Reset the symbol size, do not exit before resetting
      LDUM = PLTSTT (KSCHSZ, SZSYM)
      RETURN

      END

