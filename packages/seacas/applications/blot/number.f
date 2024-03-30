C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE NUMBER (NUMTYP, LENF, NLNKF, IX2NP, IF2EL,
     &   HIDENP, HIDEF, XN, YN, ZN, XF, YF, ZF,
     &   IELBST, IN2ELB, DODEAD, IDN2B,
     &   ISELTY, NNESEL, NESEL, BLKCOL, IDELB, MAPEL, MAPND, *)
C=======================================================================

C   --*** NUMBER *** (MESH) Number nodes and elements on mesh
C   --   Written by Amy Gilkey - revised 03/31/88
C   --
C   --NUMBER numbers nodes and/or elements.  Node numbers are in white
C   --to the top and left of the node.  Element numbers are in the
C   --element block color in the center of the element.
C   --
C   --Parameters:
C   --   NUMTYP - IN - the numbering type (as in MSHNUM of /MSHOPT/)
C   --   LENF - IN - the cumulative face counts by element block
C   --   NLNKF - IN - the number of nodes per face
C   --   IX2NP - IN - the node number for each mesh index
C   --   IF2EL - IN - the element number of each face
C   --   HIDENP(i) - IN - true iff node i is hidden (3D only)
C   --   HIDEF(i) - IN - true iff face i is hidden (3D only)
C   --   XN, YN, ZN - IN - the nodal coordinates
C   --   XF, YF, ZF - IN - the element center coordinates
C   --   IELBST - IN - the element block status (>0 if selected)
C   --   IN2ELB - IN - the element block for each node;
C   --      <0 if not in any selected element block
C   --      =0 if in more than one selected element block
C   --   DODEAD - IN - number dead nodes iff true
C   --   IDN2B - IN - the element block for each dead node; dead if >= 0
C   --   ISELTY - IN - the type of selected nodes/elements
C   --      (as in MSHNUM of /MSHOPT/)
C   --   NNESEL - IN - the number of selected nodes/elements
C   --   NESEL - IN - the selected nodes/elements
C   --   BLKCOL - IN/OUT - the user selected colors of the element blocks.
C   --                    BLKCOL(0) = 1 if the user defined material
C   --                                colors should be used in mesh plots.
C   --                              = -1 if program selected colors should
C   --                                be used.
C   --                    BLKCOL(i) = the user selected color of element
C   --                               block i:
C   --                                  -2 - no color selected by user.
C   --                                  -1 - black
C   --                                   0 - white
C   --                                   1 - red
C   --                                   2 - green
C   --                                   3 - yellow
C   --                                   4 - blue
C   --                                   5 - cyan
C   --                                   6 - magenta
C   --   * - return statement if the cancel function is active
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/
C   --   Uses IS3DIM, NUMNPF of /D3NUMS/

      PARAMETER (KHCHSZ=1, KSCHSZ=2)

      include 'debug.blk'
      include 'dbnums.blk'
      include 'd3nums.blk'

      CHARACTER*(*) NUMTYP
      INTEGER LENF(0:NELBLK)
      INTEGER NLNKF(NELBLK)
      INTEGER IX2NP(NUMNPF)
      INTEGER IF2EL(*)
      LOGICAL HIDENP(NUMNPF)
      LOGICAL HIDEF(*)
      REAL XN(NUMNPF), YN(NUMNPF), ZN(NUMNPF)
      REAL XF(*), YF(*), ZF(*)
      INTEGER IELBST(NELBLK)
      INTEGER IN2ELB(NUMNPF)
      LOGICAL DODEAD
      INTEGER IDN2B(NUMNPF)
      CHARACTER*(*) ISELTY
      INTEGER NESEL(*)
      INTEGER BLKCOL(0:NELBLK)
      INTEGER IDELB(*)
      INTEGER MAPEL(*), MAPND(*)

      LOGICAL PLTGTT, LDUM
      LOGICAL GRABRT
      LOGICAL EXISTS
      LOGICAL SOFTCH
      logical logt
      LOGICAL OK1, OK2
      CHARACTER*10 ISTR
      INTEGER IFACES(10)

      EXISTS (M) = (MOD(M,2) .NE. 0)

C   --Get software character flag for current device
      CALL GRGPARD ('SOFTCHAR', 0, SOFTCH, ISTR)

      IF ((NUMTYP .EQ. 'SELECTED') .AND. (ISELTY .EQ. 'NODE')) THEN

C      --Number selected nodes

         IF (SOFTCH) THEN
            DXO = .001
            DYO = .002
         ELSE
            DXO = 0.0
            DYO = 0.0
         END IF

         CALL UGRCOL (0, BLKCOL)

         DO 100 INE = 1, NNESEL
            INP = NESEL(INE)
            IF (IS3DIM) THEN
               IF (HIDENP(INP)) GOTO 100
            END IF

            if (IN2ELB(inp) .GE. 0) then
               logt = .true.
            else if (dodead) then
               if (IDN2B(inp) .GE. 0) logt = .true.
            end if

            if (logt) then
               IF (GRABRT ()) RETURN 1
               NNUM = MAPND(INP)
               CALL INTSTR (1, 0, NNUM, ISTR, LSTR)
               CALL MP2PT (1, XN(INP), YN(INP), DX0, DY0, MASK)
               IF (EXISTS (MASK))
     &            CALL GRTEXT (DX0+DXO, DY0+DYO, ISTR(:LSTR))
            END IF
  100    CONTINUE

         CALL PLTFLU

      ELSE IF ((NUMTYP .EQ. 'NODE') .OR. (NUMTYP .EQ. 'ALL')) THEN

C      --Number all nodes

         IF (SOFTCH) THEN
            DXO = .001
            DYO = .002
         ELSE
            DXO = 0.0
            DYO = 0.0
         END IF

         CALL UGRCOL (0, BLKCOL)

         DO 110 INP = 1, NUMNPF
            IF (IS3DIM) THEN
               IF (HIDENP(INP)) GOTO 110
            END IF

            if (IN2ELB(inp) .GE. 0) then
               logt = .true.
            else if (dodead) then
               if (IDN2B(inp) .GE. 0) logt = .true.
            end if

            if (logt) then
               IF (GRABRT ()) RETURN 1
               NNUM = MAPND(INP)
               CALL INTSTR (1, 0, NNUM, ISTR, LSTR)
               CALL MP2PT (1, XN(INP), YN(INP), DX0, DY0, MASK)
               IF (EXISTS (MASK))
     &            CALL GRTEXT (DX0+DXO, DY0+DYO, ISTR(:LSTR))
            END IF
  110    CONTINUE

         CALL PLTFLU
      END IF

      IF (ISELTY .EQ. 'NODE') THEN

C      --Draw arrows

         CALL GRCOLR (3)

         OK2 = .FALSE.
         DO 120 INE = 1, NNESEL
            OK1 = OK2
            OK2 = .FALSE.

            INP = NESEL(INE)
            IF (IS3DIM) THEN
               IF (HIDENP(INP)) GOTO 120
            END IF

            if (IN2ELB(inp) .GE. 0) then
               logt = .true.
            else if (dodead) then
               if (IDN2B(inp) .GE. 0) logt = .true.
            end if

            if (logt) then
               X2 = XN(INP)
               Y2 = YN(INP)
               OK2 = .TRUE.

               IF (OK1 .AND. OK2) THEN
                  IF (GRABRT ()) RETURN 1
                  CALL MP2VC (1, X1, Y1, X2, Y2,
     &               DX1, DY1, DX2, DY2, MASK)
                  IF (EXISTS (MASK))
     &               CALL PLTARR (DX1, DY1, DX2, DY2, .5, .0075)
               END IF

               X1 = X2
               Y1 = Y2
            END IF
  120    CONTINUE

         CALL PLTFLU
      END IF

      IF ((NUMTYP .EQ. 'SELECTED') .AND. (ISELTY .EQ. 'ELEMENT')) THEN

         IF (SOFTCH) THEN
            LDUM = PLTGTT (KSCHSZ, CHSIZ)
         ELSE
            LDUM = PLTGTT (KHCHSZ, CHSIZ)
         END IF
         DXO = 0.0
         DYO = -.5 * CHSIZ

C      --Number selected elements by element block color

         DO 140 IELB = 1, NELBLK
            IF (IELBST(IELB) .GT. 0) THEN

c               CALL UGRCOL (IDELB(IELB), BLKCOL)
               CALL UGRCOL (IELB, BLKCOL)

               DO 130 IFAC = LENF(IELB-1)+1, LENF(IELB)
                  IF (IS3DIM) THEN
                     IF (HIDEF(IFAC)) GOTO 130
                  END IF

                  IEL = IF2EL(IFAC)
                  IF (LOCINT (IEL, NNESEL, NESEL) .GT. 0) THEN
                     IF (GRABRT ()) RETURN 1
                     if ((cdebug .eq. 'HIDDEN')
     &                  .or. (cdebug .eq. 'NUMBER')) iel = ifac
                     NNUM = MAPEL(IEL)
                     CALL INTSTR (1, 0, NNUM, ISTR, LSTR)
                     CALL MP2PT (1, XF(IFAC), YF(IFAC), DX0, DY0, MASK)
                     IF (EXISTS (MASK))
     &                  CALL GRTEXC (DX0+DXO, DY0+DYO, ISTR(:LSTR))
                  END IF
  130          CONTINUE
            END IF
  140    CONTINUE

         CALL PLTFLU

      ELSE IF ((NUMTYP .EQ. 'ELEMENT') .OR. (NUMTYP .EQ. 'ALL')) THEN

         IF (SOFTCH) THEN
            LDUM = PLTGTT (KSCHSZ, CHSIZ)
         ELSE
            LDUM = PLTGTT (KHCHSZ, CHSIZ)
         END IF
         DXO = 0.0
         DYO = -.5 * CHSIZ

C      --Number elements by element block color

         DO 160 IELB = 1, NELBLK
            IF (IELBST(IELB) .GT. 0) THEN

c               CALL UGRCOL (IDELB(IELB), BLKCOL)
               CALL UGRCOL (IELB, BLKCOL)

               DO 150 IFAC = LENF(IELB-1)+1, LENF(IELB)
                  IF (IS3DIM) THEN
                     IF (HIDEF(IFAC)) GOTO 150
                  END IF

                  IEL = IF2EL(IFAC)
                  IF (GRABRT ()) RETURN 1
                  if ((cdebug .eq. 'HIDDEN')
     &               .or. (cdebug .eq. 'NUMBER')) iel = ifac
                  CALL INTSTR (1, 0, MAPEL(IEL), ISTR, LSTR)
                  CALL MP2PT (1, XF(IFAC), YF(IFAC), DX0, DY0, MASK)
                  IF (EXISTS (MASK))
     &               CALL GRTEXC (DX0+DXO, DY0+DYO, ISTR(:LSTR))
  150          CONTINUE
            END IF
  160    CONTINUE

         CALL PLTFLU
      END IF

      IF (ISELTY .EQ. 'ELEMENT') THEN

C      --Draw arrows

         DYO = -.5 * CHSIZ - .002

         CALL GRCOLR (3)

         OK2 = .FALSE.
         DO 180 INE = 1, NNESEL
            OK1 = OK2
            OK2 = .FALSE.

            IEL = NESEL(INE)
            NQARY = 10
            CALL FNDE2F (IEL, LENF, IF2EL, NQARY, IFACES, IELB)

            IF ((NQARY .GT. 0) .AND. (IELBST(IELB) .GT. 0)) THEN
               NOK = 0
               X2 = 0.0
               Y2 = 0.0
               DO 170 N = 1, NQARY
                  IFAC = IFACES(N)
                  IF (IS3DIM) THEN
                     IF (HIDEF(IFAC)) GOTO 170
                  END IF
                  X2 = X2 + XF(IFAC)
                  Y2 = Y2 + YF(IFAC)
                  NOK = NOK + 1
                  OK2 = .TRUE.
  170          CONTINUE
               IF (NOK .GT. 1) THEN
                  X2 = X2 / NOK
                  Y2 = Y2 / NOK
               END IF

               IF (OK1 .AND. OK2) THEN
                  IF (GRABRT ()) RETURN 1
                  CALL MP2VC (1, X1, Y1, X2, Y2,
     &               DX1, DY1, DX2, DY2, MASK)
                  IF (EXISTS (MASK))
     &               CALL PLTARR (DX1, DY1+DYO, DX2, DY2+DYO, .5, .0075)
               END IF

               X1 = X2
               Y1 = Y2
            END IF
  180    CONTINUE

         CALL PLTFLU
      END IF

      RETURN
      END
