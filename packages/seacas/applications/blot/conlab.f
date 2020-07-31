C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CONLAB (ICNTR, CNTR, NHIT, LABINC, LINSET,
     &   VARNP, XN, YN, ZN, IN2ELB, *)
C=======================================================================

C   --*** CONLAB *** (DETOUR) Label contour line
C   --   Written by Amy Gilkey - revised 06/01/87
C   --   D. P. Flanagan, 03/30/83
C   --
C   --CONLAB labels a contour line that intersects the given line with the
C   --contour letter.  If LABINC is greater than one, only every LABINCth
C   --crossing is labeled.
C   --
C   --Parameters:
C   --   ICNTR - IN - the contour label number
C   --   CNTR - IN - the contour value
C   --   NHIT - IN/OUT - the number of "hits" so far, initialize to zero
C   --   LABINC - IN - the number of "hits" to skip before labelling
C   --   LINSET - IN - the line set
C   --   VARNP - IN - the contour function values
C   --   XN, YN, ZN - IN - the nodal coordinates
C   --   IN2ELB - IN - the element block for each node;
C   --      <0 if not in any selected element block
C   --      =0 if in more than one selected element block
C   --   * - return statement if the cancel function is active
C   --
C   --Common Variables:
C   --   Uses IS3DIM of /D3NUMS/

      PARAMETER (DXO = .0015, DYO = .0015)

      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM

      INTEGER LINSET(3)
      REAL VARNP(*)
      REAL XN(*), YN(*), ZN(*)
      INTEGER IN2ELB(*)

      LOGICAL INTERP_BL
      LOGICAL GRABRT
      LOGICAL EXISTS

      EXISTS (M) = (MOD(M,2) .NE. 0)

      IF (IS3DIM) THEN
C      --Skip invisible line
         IF (LINSET(3) .EQ. 0) GOTO 100
      END IF

      N1 = LINSET(1)
      N2 = LINSET(2)
      IF ((IN2ELB(N1) .GE. 0) .AND. (IN2ELB(N2) .GE. 0)) THEN

         IF (INTERP_BL (CNTR, VARNP(N1), VARNP(N2), PSI)) THEN
            NHIT = NHIT + 1
            IF (NHIT .GE. LABINC) THEN
               IF (GRABRT ()) RETURN 1
               IF (IS3DIM) THEN
C               --Replace invisible node with partial line node
                  IF (ABS (LINSET(3)) .GT. 1)
     &               N2 = ABS (LINSET(3))
               END IF
               X0 = XN(N1) * (1.-PSI) + XN(N2) * PSI
               Y0 = YN(N1) * (1.-PSI) + YN(N2) * PSI
               CALL MP2PT (1, X0, Y0, DX0, DY0, MASK)
               IF (EXISTS (MASK)) CALL PLTXTS (DX0+DXO, DY0+DYO,
     &            CHAR (ICHAR('A')-1+ICNTR))
               NHIT = 0
            END IF
         END IF
      END IF

  100 CONTINUE
      RETURN
      END
