C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DOFNC2 (FUNC, PARTYP, LENARY, NELBLK, IXELB, ISEVOK,
     &   NUMIX, IXNODE, IXELEM, PARM1, PARM2, RESULT)
C=======================================================================

C   --*** DOFNC2 *** (ALGEBRA) Perform function (2 parameters)
C   --   Written by Amy Gilkey - revised 05/17/88
C   --
C   --DOFNC2 performs an intrinsic or external function with two parameters
C   --on each array value.
C   --
C   --Parameters:
C   --   FUNC   - IN  - the intrinsic or external function
C   --   PARTYP - IN  - the parameter type (element must be handled differently)
C   --   LENARY - IN  - the length of the parameter arrays
C   --   NELBLK - IN  - the number of element blocks
C   --   IXELB  - IN  - the cumulative element counts for each element block
C   --   ISEVOK - IN  - the element variable truth table for each element block
C   --   NUMIX  - IN  - the number of selected values; <0 if all
C   --   IXNODE - IN  - the indices of the selected nodes (only if NUMIX>=0)
C   --   IXELEM - IN  - the indices of the selected elements (only if NUMIX>=0)
C   --   PARM1  - IN  - parameter 1 array
C   --   PARM2  - IN  - parameter 2 array
C   --   RESULT - OUT - the returned result array

      CHARACTER PARTYP
      INTEGER IXELB(0:NELBLK)
      LOGICAL ISEVOK(*)
      INTEGER IXNODE(*)
      INTEGER IXELEM(*)
      REAL PARM1(*), PARM2(*)
      REAL RESULT(*)

      IF (NUMIX .GE. 0) THEN
         IF (PARTYP .NE. 'E') THEN
            DO 100 I = 1, NUMIX
               J = IXNODE(I)
               RESULT(J) = FUNC (PARM1(J), PARM2(J))
  100       CONTINUE
         ELSE
            DO 120 IELB = 1, NELBLK
               IF (ISEVOK(IELB)) THEN
                  DO 110 I = IXELB(IELB-1)+1, IXELB(IELB)
                     J = IXELEM(I)
                     RESULT(J) = FUNC (PARM1(J), PARM2(J))
  110             CONTINUE
               END IF
  120       CONTINUE
         END IF

      ELSE
C        NUMIX < 0; All values selected
C        NUMEL = NUMELO; NUMNP = NUMNPO; IXELB(J) = IXELBO(J) J=1, NELBLK
         IF (PARTYP .NE. 'E') THEN
            DO 130 J = 1, LENARY
               RESULT(J) = FUNC (PARM1(J), PARM2(J))
  130       CONTINUE
         ELSE
            DO 150 IELB = 1, NELBLK
               IF (ISEVOK(IELB)) THEN
                  DO 140 J = IXELB(IELB-1)+1, IXELB(IELB)
                     RESULT(J) = FUNC (PARM1(J), PARM2(J))
  140             CONTINUE
               END IF
  150       CONTINUE
         END IF
      END IF

      RETURN
      END
