C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MRKESS (LENF, NLNKF, LINKF,
     &   HIDENP, HIDEF, XN, YN, ZN, XF, YF, ZF,
     &   IELBST, IX2NP, IN2ELB, DODEAD, IDN2B, IF2EL, IE2ELB,
     &   NESSET, ISSESS, IDESS, NEESS, NNESS, IXEESS, IXNESS,
     &   LTEESS, LTNESS, *)
C=======================================================================

C   --*** MRKESS *** (SETS) Mark side set on mesh
C   --   Written by Amy Gilkey - revised 04/08/88
C   --
C   --MRKESS marks all the nodes in selected side sets on the mesh.
C   --The nodes in each side are connected in the color assigned
C   --to the side set.
C   --
C   --Parameters:
C   --   LENF - IN - the cumulative face counts by element block
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF - IN - the connectivity for all faces
C   --   HIDENP(i) - IN - true iff node i is hidden (3D only)
C   --   HIDEF(i) - IN - true iff face i is hidden (3D only)
C   --   XN, YN, ZN - IN - the nodal coordinates
C   --   XF, YF, ZF - IN - the element center coordinates
C   --   IELBST - IN - the element block status (>0 if selected)
C   --   IX2NP - IN - the node number for each mesh index
C   --   IN2ELB - IN - the element block for each node;
C   --      <0 if not in any selected element block
C   --      =0 if in more than one selected element block
C   --   DODEAD - IN - mark dead nodes iff true
C   --   IDN2B - IN - the element block for each dead node; dead if >= 0
C   --   IF2EL - IN - the element number of each face
C   --   IE2ELB - IN - the element block for each element
C   --   NESSET - IN - the number of selected side sets
C   --   ISSESS - IN - the indices of the selected side sets
C   --   IDESS - IN - the side set ID for each set
C   --   NEESS - IN - the number of elements for each set
C   --   NNESS - IN - the number of nodes for each set
C   --   IXEESS - IN - the index of the first element for each set
C   --   IXNESS - IN - the index of the first node for each set
C   --   LTEESS - IN - the elements for all sets
C   --   LTNESS - IN - the nodes for all sets
C   --   * - return statement if the cancel function is active
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/
C   --   Uses IS3DIM of /D3NUMS/

      PARAMETER (KHCHSZ=1, KSCHSZ=2)

      common /debugc/ cdebug
      common /debugn/ idebug
      character*8 cdebug

      include 'dbnums.blk'
      include 'dbnumgq.blk'
      include 'd3nums.blk'

      INTEGER LENF(0:NELBLK)
      INTEGER NLNKF(NELBLK)
      INTEGER LINKF(*)
      LOGICAL HIDENP(*)
      LOGICAL HIDEF(*)
      REAL XN(*), YN(*), ZN(*)
      REAL XF(*), YF(*), ZF(*)
      INTEGER IELBST(NELBLK)
      INTEGER IX2NP(NUMNPF)
      INTEGER IN2ELB(NUMNPF)
      LOGICAL DODEAD
      INTEGER IDN2B(NUMNPF)
      INTEGER IF2EL(*)
      INTEGER IE2ELB(NUMEL)
      INTEGER ISSESS(*)
      INTEGER IDESS(*)
      INTEGER NEESS(*)
      INTEGER NNESS(*)
      INTEGER IXEESS(*)
      INTEGER IXNESS(*)
      INTEGER LTEESS(*)
      INTEGER LTNESS(*)

      LOGICAL GRABRT
      LOGICAL PLTGTT, PLTSTT, LDUM
      LOGICAL EXISTS
      INTEGER IFACES(10)
      LOGICAL NPOK
      INTEGER NOK(10)
      LOGICAL ARROW

      EXISTS (M) = (MOD(M,2) .NE. 0)

C   --Set the symbol size, do not exit before resetting
      LDUM = PLTGTT (KSCHSZ, SZSYM)
      LDUM = PLTSTT (KSCHSZ, 1.0 * SZSYM)

      DO 150 IX = 1, NESSET
         IESS = ISSESS(IX)
         IF (NEESS(IESS) .LE. 0) THEN
            GOTO 150
         END IF
         NPPSID = NNESS(IESS) / NEESS(IESS)
         NPOK = (NNESS(IESS) .EQ. (NEESS(IESS) * NPPSID))
         IF (.NOT. NPOK) THEN
            CALL PRTERR ('ERROR',
     &         'Number of elements and nodes do not match')
            GOTO 150
         END IF
         IF (NPPSID .GT. 10) THEN
            CALL PRTERR ('ERROR', 'Number of nodes per side > 10')
            GOTO 150
         END IF

         CALL GRCOLR (IESS)

         IXE0 = IXEESS(IESS) - 1

         DO 140 IEE = 1, NEESS(IESS)
            IEL = LTEESS(IXE0+IEE)
            IF (IELBST(IE2ELB(IEL)) .GT. 0) THEN
               IXN0 = IXNESS(IESS) + (IEE-1) * NPPSID - 1

               DO 100 INE = 1, NPPSID
                  NOK(INE) = 0
                  INP = LTNESS(IXN0+INE)
                  IF (IS3DIM) THEN
                     IF (HIDENP(INP)) GOTO 100
                  END IF

                  if (IN2ELB(inp) .GE. 0) then
                     nok(ine) = inp
                  else if (dodead) then
                     if (idn2b(inp) .ge. 0) nok(ine) = inp
                  end if
 100           CONTINUE

               IF (NPPSID .GT. 2) THEN
                  N2 = NOK(NPPSID)
               ELSE
                  N2 = 0
               END IF

               ARROW = .TRUE.
               DO 110 INE = 1, NPPSID
                  N1 = N2
                  N2 = NOK(INE)
                  IF ((N1 .GT. 0) .AND. (N2 .GT. 0)) THEN
                     IF (GRABRT ()) GOTO 160
                     IF (ARROW) THEN
                        CALL MP2VC (1, XN(N1), YN(N1), XN(N2), YN(N2),
     &                     DX1, DY1, DX2, DY2, MASK)
                        IF (EXISTS (MASK))
     &                    CALL PLTARR (DX1, DY1, DX2, DY2, 1.0, 0.015)
                        ARROW = .FALSE.
                     ELSE
                        CALL MPD2VC (1, XN(N1), YN(N1), XN(N2), YN(N2))
                     END IF
C????                  ELSE IF (N2 .GT. 0) THEN
C????                     I = INE + 1
C????                     IF (I .GT. NPPSID) I = 1
C????                     IF (NOK(I) .LE. 0) THEN
C????                        CALL MPD2SY (1, XN(N2), YN(N2), '\CX')
C????                     END IF
                  END IF
  110          CONTINUE

               if (cdebug .eq. 'SSETS') then
                  CALL FNDE2F (IEL, LENF, IF2EL, NFARY, IFACES, IELB)
                  DO 130 N = 1, NFARY
                     IFAC = IFACES(N)
                     IF (IS3DIM) THEN
                        IF (HIDEF(IFAC)) GOTO 130
                     END IF
                     IXL = IDBLNK (IELB, IFAC, LENF, NLNKF)
                     DO 120 INE = 1, NPPSID
                        INP = LTNESS(IXN0+INE)
                        IF (LOCINT (INP, NLNKF(IELB), LINKF(IXL))
     &                     .LE. 0) GOTO 130
  120                CONTINUE

                     IF (GRABRT ()) GOTO 160
#if NeedsDoubleEscape
                     CALL MPD2SY (1, XF(IFAC), YF(IFAC), '\\CX')
#else
                     CALL MPD2SY (1, XF(IFAC), YF(IFAC), '\CX')
#endif
  130             CONTINUE
               end if
            END IF
  140    CONTINUE

         CALL PLTFLU
  150 CONTINUE

C   --Reset the symbol size, do not exit before resetting
      LDUM = PLTSTT (KSCHSZ, SZSYM)
      RETURN

  160 CONTINUE
C   --Reset the symbol size, do not exit before resetting
      LDUM = PLTSTT (KSCHSZ, SZSYM)
      RETURN 1
      END
