C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DOOPER (OP, PARTYP, LENARY, NELBLK, IXELB, ISEVOK,
     &                   NUMIX, IXNODE, IXELEM, PARM1, PARM2, RESULT)
C=======================================================================

C   --*** DOOPER *** (ALGEBRA) Perform operation (2 parameters)
C   --   Written by Amy Gilkey - revised 05/17/88
C   --
C   --DOOPER performs an operation with up to two parameters on each array
C   --value.
C   --
C   --Parameters:
C   --   OP     - IN  - the operator character
C   --   PARTYP - IN  - the parameter type (element must be handled differently)
C   --   LENARY - IN  - the length of the parameter arrays
C   --   NELBLK - IN  - the number of element blocks
C   --   IXELB  - IN  - the cumulative element counts for each element block
C   --   ISEVOK - IN  - the element variable truth table for each element block
C   --   NUMIX  - IN  - the number of selected values; <0 if all
C   --   IXNODE - IN  - the indices of the selected nodes (only if NUMIX>=0)
C   --   IXELEM - IN  - the indices of the selected elements (only if NUMIX>=0)
C   --   PARM1  - IN  - the parameter 1 array
C   --   PARM2  - IN  - the parameter 2 array
C   --   RESULT - OUT - the returned result array

      CHARACTER OP
      CHARACTER PARTYP
      INTEGER IXELB(0:NELBLK)
      LOGICAL ISEVOK(*)
      INTEGER IXNODE(*)
      INTEGER IXELEM(*)
      REAL PARM1(*), PARM2(*)
      REAL RESULT(*)

      IF (OP .EQ. '~') THEN

C      --Unary minus

         IF (NUMIX .GE. 0) THEN
            IF (PARTYP .NE. 'E') THEN
               DO 100 I = 1, NUMIX
                  J = IXNODE(I)
                  RESULT(J) = - PARM1(J)
  100          CONTINUE
            ELSE
               DO 120 IELB = 1, NELBLK
                  IF (ISEVOK(IELB)) THEN
                     DO 110 I = IXELB(IELB-1)+1, IXELB(IELB)
                        J = IXELEM(I)
                        RESULT(J) = - PARM1(J)
  110                CONTINUE
                  END IF
  120          CONTINUE
            END IF

         ELSE
C           NUMIX < 0; All values are selected
C           numel = numelo numnp = numnpo
            IF (PARTYP .NE. 'E') THEN
               DO 130 J = 1, LENARY
                  RESULT(J) = - PARM1(J)
  130          CONTINUE
            ELSE
               DO 150 IELB = 1, NELBLK
                  IF (ISEVOK(IELB)) THEN
                     DO 140 J = IXELB(IELB-1)+1, IXELB(IELB)
                        RESULT(J) = - PARM1(J)
  140                CONTINUE
                  END IF
  150          CONTINUE
            END IF
         END IF

      ELSE IF (OP .EQ. '+') THEN

C      --Addition

         IF (NUMIX .GE. 0) THEN
            IF (PARTYP .NE. 'E') THEN
               DO 160 I = 1, NUMIX
                  J = IXNODE(I)
                  RESULT(J) = PARM1(J) + PARM2(J)
  160          CONTINUE
            ELSE
               DO 180 IELB = 1, NELBLK
                  IF (ISEVOK(IELB)) THEN
                     DO 170 I = IXELB(IELB-1)+1, IXELB(IELB)
                        J = IXELEM(I)
                        RESULT(J) = PARM1(J) + PARM2(J)
  170                CONTINUE
                  END IF
  180          CONTINUE
            END IF

         ELSE
C           NUMIX < 0; All values are selected
C           numel = numelo numnp = numnpo
            IF (PARTYP .NE. 'E') THEN
               DO 190 J = 1, LENARY
                  RESULT(J) = PARM1(J) + PARM2(J)
  190          CONTINUE
            ELSE
               DO 210 IELB = 1, NELBLK
                  IF (ISEVOK(IELB)) THEN
                     DO 200 J = IXELB(IELB-1)+1, IXELB(IELB)
                        RESULT(J) = PARM1(J) + PARM2(J)
  200                CONTINUE
                  END IF
  210          CONTINUE
            END IF
         END IF

      ELSE IF (OP .EQ. '-') THEN

C      --Subtraction

         IF (NUMIX .GE. 0) THEN
            IF (PARTYP .NE. 'E') THEN
               DO 220 I = 1, NUMIX
                  J = IXNODE(I)
                  RESULT(J) = PARM1(J) - PARM2(J)
  220          CONTINUE
            ELSE
               DO 240 IELB = 1, NELBLK
                  IF (ISEVOK(IELB)) THEN
                     DO 230 I = IXELB(IELB-1)+1, IXELB(IELB)
                        J = IXELEM(I)
                        RESULT(J) = PARM1(J) - PARM2(J)
  230                CONTINUE
                  END IF
  240          CONTINUE
            END IF

         ELSE
C           NUMIX < 0; All values are selected
C           numel = numelo numnp = numnpo
            IF (PARTYP .NE. 'E') THEN
               DO 250 J = 1, LENARY
                  RESULT(J) = PARM1(J) - PARM2(J)
  250          CONTINUE
            ELSE
               DO 270 IELB = 1, NELBLK
                  IF (ISEVOK(IELB)) THEN
                     DO 260 J = IXELB(IELB-1)+1, IXELB(IELB)
                        RESULT(J) = PARM1(J) - PARM2(J)
  260                CONTINUE
                  END IF
  270          CONTINUE
            END IF
         END IF

      ELSE IF (OP .EQ. '*') THEN

C      --Multiplication
C      --NUMIX - the number of selected values; NUMIX<0 if all selected
         IF (NUMIX .GE. 0) THEN
            IF (PARTYP .NE. 'E') THEN
               DO 280 I = 1, NUMIX
                  J = IXNODE(I)
                  RESULT(J) = PARM1(J) * PARM2(J)
  280          CONTINUE
            ELSE
               DO 300 IELB = 1, NELBLK
                  IF (ISEVOK(IELB)) THEN
                     DO 290 I = IXELB(IELB-1)+1, IXELB(IELB)
                        J = IXELEM(I)
                        RESULT(J) = PARM1(J) * PARM2(J)
  290                CONTINUE
                  END IF
  300          CONTINUE
            END IF

         ELSE
C           NUMIX < 0; All values are selected
C           numel = numelo numnp = numnpo
            IF (PARTYP .NE. 'E') THEN
               DO 310 J = 1, LENARY
                  RESULT(J) = PARM1(J) * PARM2(J)
  310          CONTINUE
            ELSE
               DO 330 IELB = 1, NELBLK
                  IF (ISEVOK(IELB)) THEN
                     DO 320 J = IXELB(IELB-1)+1, IXELB(IELB)
                        RESULT(J) = PARM1(J) * PARM2(J)
  320                CONTINUE
                  END IF
  330          CONTINUE
            END IF
         END IF

      ELSE IF (OP .EQ. '/') THEN

C      --Division

         IF (NUMIX .GE. 0) THEN
            IF (PARTYP .NE. 'E') THEN
               DO 340 I = 1, NUMIX
                  J = IXNODE(I)
                  RESULT(J) = PARM1(J) / PARM2(J)
  340          CONTINUE
            ELSE
               DO 360 IELB = 1, NELBLK
                  IF (ISEVOK(IELB)) THEN
                     DO 350 I = IXELB(IELB-1)+1, IXELB(IELB)
                        J = IXELEM(I)
                        RESULT(J) = PARM1(J) / PARM2(J)
  350                CONTINUE
                  END IF
  360          CONTINUE
            END IF

         ELSE
C           NUMIX < 0; All values are selected
C           numel = numelo numnp = numnpo
            IF (PARTYP .NE. 'E') THEN
               DO 370 J = 1, LENARY
                  RESULT(J) = PARM1(J) / PARM2(J)
  370          CONTINUE
            ELSE
               DO 390 IELB = 1, NELBLK
                  IF (ISEVOK(IELB)) THEN
                     DO 380 J = IXELB(IELB-1)+1, IXELB(IELB)
                        RESULT(J) = PARM1(J) / PARM2(J)
  380                CONTINUE
                  END IF
  390          CONTINUE
            END IF
         END IF

      ELSE IF (OP .EQ. '^') THEN

C      --Exponentiation

         IF (NUMIX .GE. 0) THEN
            IF (PARTYP .NE. 'E') THEN
               DO 400 I = 1, NUMIX
                  J = IXNODE(I)
                  III = int(PARM2(J))
                  IF (III .EQ. PARM2(J)) THEN
                     RESULT(J) = PARM1(J)**III
                  ELSE
                     RESULT(J) = PARM1(J)**PARM2(J)
                  END IF
  400          CONTINUE
            ELSE
               DO 420 IELB = 1, NELBLK
                  IF (ISEVOK(IELB)) THEN
                     DO 410 I = IXELB(IELB-1)+1, IXELB(IELB)
                        J = IXELEM(I)
                        III = int(PARM2(J))
                        IF (III .EQ. PARM2(J)) THEN
                           RESULT(J) = PARM1(J)**III
                        ELSE
                           RESULT(J) = PARM1(J)**PARM2(J)
                        END IF
  410                CONTINUE
                  END IF
  420          CONTINUE
            END IF

         ELSE
C           NUMIX < 0; All values are selected
C           numel = numelo numnp = numnpo
            IF (PARTYP .NE. 'E') THEN
               DO 430 J = 1, LENARY
                  III = int(PARM2(J))
                  IF (III .EQ. PARM2(J)) THEN
                     RESULT(J) = PARM1(J)**III
                  ELSE
                     RESULT(J) = PARM1(J)**PARM2(J)
                  END IF
  430          CONTINUE
            ELSE
               DO 450 IELB = 1, NELBLK
                  IF (ISEVOK(IELB)) THEN
                     DO 440 J = IXELB(IELB-1)+1, IXELB(IELB)
                        III = int(PARM2(J))
                        IF (III .EQ. PARM2(J)) THEN
                           RESULT(J) = PARM1(J)**III
                        ELSE
                           RESULT(J) = PARM1(J)**PARM2(J)
                        END IF
  440                CONTINUE
                  END IF
  450          CONTINUE
            END IF
         END IF
      END IF

      RETURN
      END
