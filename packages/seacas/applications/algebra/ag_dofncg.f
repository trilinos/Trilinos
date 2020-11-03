C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DOFNCG (FUNC, PARTYP, LENARY, NELBLK, IXELB, ISEVOK,
     &                   NUMIX, IXNODE, IXELEM, PARM, RESULT)
C=======================================================================

C   --*** DOFNCG *** (ALGEBRA) Perform array-to-single-value function
C   --   Written by Amy Gilkey - revised 05/17/88
C   --
C   --DOFNCG performs an intrinsic or external function with one parameter
C   --that transforms the array into a single value (such as minimum of
C   --all array values).  Certain array values may be selected,
C   --meaning the function is only applied to those values.
C   --
C   --Parameters:
C   --   FUNC   - IN  - the intrinsic or external function
C   --   PARTYP - IN  - the parameter type (element must be handled differently)
C   --   LENARY - IN  - the length of the parameter array
C   --   NELBLK - IN  - the number of element blocks
C   --   IXELB  - IN  - the cumulative element counts for each element block
C   --   ISEVOK - IN  - the element variable truth table for each element block
C   --   NUMIX  - IN  - the number of selected values; <0 if all
C   --   IXNODE - IN  - the indices of the selected nodes (only if NUMIX>=0)
C   --   IXELEM - IN  - the indices of the selected elements (only if NUMIX>=0)
C   --   PARM   - IN  - the parameter array
C   --   RESULT - OUT - the returned result value

      CHARACTER PARTYP
      INTEGER IXELB(0:NELBLK)
      LOGICAL ISEVOK(*)
      INTEGER IXNODE(*)
      INTEGER IXELEM(*)
      REAL PARM(*)
      REAL RESULT

      LOGICAL INIT

      INIT = .TRUE.

      R = 0.0
      IF (NUMIX .GE. 0) THEN
         IF (PARTYP .NE. 'E') THEN
            DO 100 I = 1, NUMIX
               J = IXNODE(I)
               IF (INIT) THEN
                  R = PARM(J)
                  INIT = .FALSE.
               ELSE
                  R = FUNC (R, PARM(J))
               END IF
  100       CONTINUE
         ELSE
            DO 120 IELB = 1, NELBLK
               IF (ISEVOK(IELB)) THEN
                  DO 110 I = IXELB(IELB-1)+1, IXELB(IELB)
                     J = IXELEM(I)
                     IF (INIT) THEN
                        R = PARM(J)
                        INIT = .FALSE.
                     ELSE
                        R = FUNC (R, PARM(J))
                     END IF
  110             CONTINUE
               END IF
  120       CONTINUE
         END IF

      ELSE
C        NUMIX < 0; ALL VALUES SELECTED
C        NUMEL=NUMELO; NUMNP=NUMNPO; IXELB(J)=IXELBO(J), J=1,NELBLK
         IF (PARTYP .NE. 'E') THEN
            R = PARM(1)
            DO 130 J = 2, LENARY
               R = FUNC (R, PARM(J))
  130       CONTINUE
         ELSE
            DO 150 IELB = 1, NELBLK
               IF (ISEVOK(IELB)) THEN
                  DO 140 J = IXELB(IELB-1)+1, IXELB(IELB)
                     IF (INIT) THEN
                        R = PARM(J)
                        INIT = .FALSE.
                     ELSE
                        R = FUNC (R, PARM(J))
                     END IF
  140             CONTINUE
               END IF
  150       CONTINUE
         END IF
      END IF

      RESULT = R

      RETURN
      END
