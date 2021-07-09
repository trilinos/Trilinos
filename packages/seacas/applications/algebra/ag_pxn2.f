C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PXN2 (RETMIN, PARTYP, LENARY, NELBLK, IXELB, ISEVOK,
     &           NUMIX, IXNODE, IXELEM, SIG11, SIG22, SIG12, RESULT)
C=======================================================================

C   --*** PXN2 *** (ALGEBRA) Calculate 2x2 principal values
C   --   Written by Amy Gilkey - revised 05/17/88
C   --
C   --PXN2 determines minimum or maximum principal values for a 2 by 2
C   --array.
C   --
C   --Parameters:
C   --   RETMIN - IN - true iff minimum versus maximum to be calculated
C   --   PARTYP - IN - the parameter type (element must be handled differently)
C   --   LENARY - IN - the length of the parameter arrays
C   --   NELBLK - IN - the number of element blocks
C   --   IXELB  - IN - the cumulative element counts for each element block
C   --   ISEVOK - IN - the element variable truth table for each element block
C   --   NUMIX  - IN - the number of selected values; <0 if all
C   --   IXNODE - IN - the indices of the selected nodes (only if NUMIX >= 0)
C   --   IXELEM - IN - the indices of the selected elements (only if NUMIX >= 0)
C   --   SIG11, SIG22, SIG12 - IN - the tensor components
C   --   RESULT - OUT - the returned principal values

      LOGICAL RETMIN
      CHARACTER PARTYP
      INTEGER IXELB(0:NELBLK)
      LOGICAL ISEVOK(*)
      INTEGER IXNODE(*)
      INTEGER IXELEM(*)
      REAL SIG11(*), SIG22(*), SIG12(*)
      REAL RESULT(*)

      IF (NUMIX .GE. 0) THEN
         IF (PARTYP .NE. 'E') THEN
            DO 100 I = 1, NUMIX
               J = IXNODE(I)
               SROOT = SQRT
     &            ((0.5 * (SIG11(J)-SIG22(J)))**2 + SIG12(J)**2)

               IF (RETMIN) THEN
                  RESULT(J) = 0.5 * (SIG11(J)+SIG22(J)) - SROOT
               ELSE
                  RESULT(J) = 0.5 * (SIG11(J)+SIG22(J)) + SROOT
               END IF
  100       CONTINUE
         ELSE
            DO 120 IELB = 1, NELBLK
               IF (ISEVOK(IELB)) THEN
                  DO 110 I = IXELB(IELB-1)+1, IXELB(IELB)
                     J = IXELEM(I)
                     SROOT = SQRT
     &                  ((0.5 * (SIG11(J)-SIG22(J)))**2 + SIG12(J)**2)

                     IF (RETMIN) THEN
                        RESULT(J) = 0.5 * (SIG11(J)+SIG22(J)) - SROOT
                     ELSE
                        RESULT(J) = 0.5 * (SIG11(J)+SIG22(J)) + SROOT
                     END IF
  110             CONTINUE
               END IF
  120       CONTINUE
         END IF

      ELSE
C        NUMIX < 0; ALL VALUES SELECTED
C        NUMEL=NUMELO; NUMNP=NUMNPO; IXELB(J)=IXELBO(J), J=1,NELBLK
         IF (PARTYP .NE. 'E') THEN
            DO 130 J = 1, LENARY
               SROOT = SQRT
     &            ((0.5 * (SIG11(J)-SIG22(J)))**2 + SIG12(J)**2)

               IF (RETMIN) THEN
                  RESULT(J) = 0.5 * (SIG11(J)+SIG22(J)) - SROOT
               ELSE
                  RESULT(J) = 0.5 * (SIG11(J)+SIG22(J)) + SROOT
               END IF
  130       CONTINUE

         ELSE
            DO 150 IELB = 1, NELBLK
               IF (ISEVOK(IELB)) THEN
                  DO 140 J = IXELB(IELB-1)+1, IXELB(IELB)
                     SROOT = SQRT
     &                  ((0.5 * (SIG11(J)-SIG22(J)))**2 + SIG12(J)**2)

                     IF (RETMIN) THEN
                        RESULT(J) = 0.5 * (SIG11(J)+SIG22(J)) - SROOT
                     ELSE
                        RESULT(J) = 0.5 * (SIG11(J)+SIG22(J)) + SROOT
                     END IF
  140             CONTINUE
               END IF
  150       CONTINUE
         END IF
      END IF

      RETURN
      END
