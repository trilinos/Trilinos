C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE TMAG (PARTYP, LENARY, NELBLK, IXELB, ISEVOK,
     &                 SIG1, SIG2, SIG3, SIG12, SIG23, SIG31,
     &                 NUMIX, IXNODE, IXELEM, RESULT)
C=======================================================================

C   --*** TMAG *** (ALGEBRA) Calculate tensor magnitude
C   --   Written by Amy Gilkey - revised 12/03/87
C   --
C   --TMAG calculates the tensor magnitudes.
C   --
C   --Parameters:
C   --   PARTYP - IN - the parameter type (element must be handled differently)
C   --   LENARY - IN - the length of the parameter arrays
C   --   NELBLK - IN - the number of element blocks
C   --   IXELB  - IN - the cumulative element counts for each element block
C   --   ISEVOK - IN - the element variable truth table for each element block
C   --   SIG1, SIG2, SIG3, SIG12, SIG23, SIG31 - IN - the tensor components
C   --   NUMIX  - IN  - the number of selected values; <0 if all
C   --   IXNODE - IN  - the indices of the selected nodes (only if NUMIX>=0)
C   --   IXELEM - IN  - the indices of the selected elements (only if NUMIX>=0)
C   --   RESULT - OUT - the returned tensor magnitudes

      CHARACTER PARTYP
      INTEGER IXELB(0:NELBLK)
      LOGICAL ISEVOK(*)
      REAL SIG1(*), SIG2(*), SIG3(*), SIG12(*), SIG23(*), SIG31(*)
      INTEGER NUMIX
      INTEGER IXNODE(*)
      INTEGER IXELEM(*)
      REAL RESULT(*)

      IF (NUMIX .GE. 0) THEN
C        COMPUTE TMAG ONLY FOR ELEMENTS/NODES THAT EXIST
         IF (PARTYP .NE. 'E') THEN
            DO 50 I = 1, NUMIX
               J = IXNODE(I)
               RESULT(J) = SQRT ((SIG1(J)-SIG2(J))**2
     &            + (SIG2(J)-SIG3(J))**2
     &            + (SIG3(J)-SIG1(J))**2
     &            + 6 * (SIG12(J)**2 + SIG23(J)**2 + SIG31(J)**2))
   50       CONTINUE
         ELSE
            DO 75 IELB = 1, NELBLK
               IF (ISEVOK(IELB)) THEN
                  DO 95 I = IXELB(IELB-1)+1, IXELB(IELB)
                     J = IXELEM(I)
                     RESULT(J) = SQRT ((SIG1(J)-SIG2(J))**2
     &                  + (SIG2(J)-SIG3(J))**2
     &                  + (SIG3(J)-SIG1(J))**2
     &                  + 6 * (SIG12(J)**2 + SIG23(J)**2 + SIG31(J)**2))
   95             CONTINUE
               END IF
   75       CONTINUE
         END IF
      ELSE
C        NUMIX < 0; All Values are selected
C        NUMEL = NUMELO; NUMNP = NUMNPO; IXELB(J) = IXELBO(J) J=1,NELBLK
         IF (PARTYP .NE. 'E') THEN
            DO 100 J = 1, LENARY
               RESULT(J) = SQRT ((SIG1(J)-SIG2(J))**2
     &            + (SIG2(J)-SIG3(J))**2
     &            + (SIG3(J)-SIG1(J))**2
     &            + 6 * (SIG12(J)**2 + SIG23(J)**2 + SIG31(J)**2))
  100       CONTINUE
         ELSE
            DO 120 IELB = 1, NELBLK
               IF (ISEVOK(IELB)) THEN
                  DO 110 J = IXELB(IELB-1)+1, IXELB(IELB)
                     RESULT(J) = SQRT ((SIG1(J)-SIG2(J))**2
     &                  + (SIG2(J)-SIG3(J))**2
     &                  + (SIG3(J)-SIG1(J))**2
     &                  + 6 * (SIG12(J)**2 + SIG23(J)**2 + SIG31(J)**2))
  110             CONTINUE
               END IF
  120       CONTINUE
         END IF
      END IF

      RETURN
      END
