C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DOIF (TEST, PARTYP, LENARY, NELBLK, IXELB, ISEVOK,
     &           NUMIX, IXNODE, IXELEM, PARMC, PARMT, PARMF, RESULT)
C=======================================================================

C   --*** DOIF *** (ALGEBRA) Perform operation (2 parameters)
C   --   Written by Amy Gilkey - revised 05/17/88
C   --
C   --DOIF performs an operation with up to two parameters on each array
C   --value.
C   --
C   --Parameters:
C   --   TEST   - IN  - the operator character
C   --   PARTYP - IN  - the parameter type (element must be handled differently)
C   --   LENARY - IN  - the length of the parameter arrays
C   --   NELBLK - IN  - the number of element blocks
C   --   IXELB  - IN  - the cumulative element counts for each element block
C   --   ISEVOK - IN  - the element variable truth table for each element block
C   --   NUMIX  - IN  - the number of selected values; <0 if all
C   --   IXNODE - IN  - the indices of the selected nodes (only if NUMIX>=0)
C   --   IXELEM - IN  - the indices of the selected elements (only if NUMIX>=0)
C   --   PARMC  - IN  - the condition arrays
C   --   PARMT, PARMF - IN - the true and false result arrays
C   --   RESULT - OUT - the returned result array

      CHARACTER*(*) TEST
      CHARACTER PARTYP
      INTEGER IXELB(0:NELBLK)
      LOGICAL ISEVOK(*)
      INTEGER IXNODE(*)
      INTEGER IXELEM(*)
      REAL PARMC(*)
      REAL PARMT(*), PARMF(*)
      REAL RESULT(*)

      IF (TEST .EQ. 'IFLZ') THEN

         IF (NUMIX .GE. 0) THEN
            IF (PARTYP .NE. 'E') THEN
               DO 100 I = 1, NUMIX
                  J = IXNODE(I)
                  IF (PARMC(J) .LT. 0.0) THEN
                     RESULT(J) = PARMT(J)
                  ELSE
                     RESULT(J) = PARMF(J)
                  END IF
  100          CONTINUE
            ELSE
               DO 120 IELB = 1, NELBLK
                  IF (ISEVOK(IELB)) THEN
                     DO 110 I = IXELB(IELB-1)+1, IXELB(IELB)
                        J = IXELEM(I)
                        IF (PARMC(J) .LT. 0.0) THEN
                           RESULT(J) = PARMT(J)
                        ELSE
                           RESULT(J) = PARMF(J)
                        END IF
  110                CONTINUE
                  END IF
  120          CONTINUE
            END IF

         ELSE
C           NUMIX < 0; ALL VALUES SELECTED
C           NUMEL=NUMELO; NUMNP=NUMNPO; IXELB(J)=IXELBO(J), J=1,NELBLK
            IF (PARTYP .NE. 'E') THEN
               DO 130 J = 1, LENARY
                  IF (PARMC(J) .LT. 0.0) THEN
                     RESULT(J) = PARMT(J)
                  ELSE
                     RESULT(J) = PARMF(J)
                  END IF
  130          CONTINUE
            ELSE
               DO 150 IELB = 1, NELBLK
                  IF (ISEVOK(IELB)) THEN
                     DO 140 J = IXELB(IELB-1)+1, IXELB(IELB)
                        IF (PARMC(J) .LT. 0.0) THEN
                           RESULT(J) = PARMT(J)
                        ELSE
                           RESULT(J) = PARMF(J)
                        END IF
  140                CONTINUE
                  END IF
  150          CONTINUE
            END IF
         END IF

      ELSE IF (TEST .EQ. 'IFEZ') THEN

         IF (NUMIX .GE. 0) THEN
            IF (PARTYP .NE. 'E') THEN
               DO 160 I = 1, NUMIX
                  J = IXNODE(I)
                  IF (PARMC(J) .EQ. 0.0) THEN
                     RESULT(J) = PARMT(J)
                  ELSE
                     RESULT(J) = PARMF(J)
                  END IF
  160          CONTINUE
            ELSE
               DO 180 IELB = 1, NELBLK
                  IF (ISEVOK(IELB)) THEN
                     DO 170 I = IXELB(IELB-1)+1, IXELB(IELB)
                        J = IXELEM(I)
                        IF (PARMC(J) .EQ. 0.0) THEN
                           RESULT(J) = PARMT(J)
                        ELSE
                           RESULT(J) = PARMF(J)
                        END IF
  170                CONTINUE
                  END IF
  180          CONTINUE
            END IF

         ELSE
C           NUMIX < 0; ALL VALUES SELECTED
C           NUMEL=NUMELO; NUMNP=NUMNPO; IXELB(J)=IXELBO(J), J=1,NELBLK
            IF (PARTYP .NE. 'E') THEN
               DO 190 J = 1, LENARY
                  IF (PARMC(J) .EQ. 0.0) THEN
                     RESULT(J) = PARMT(J)
                  ELSE
                     RESULT(J) = PARMF(J)
                  END IF
  190          CONTINUE
            ELSE
               DO 210 IELB = 1, NELBLK
                  IF (ISEVOK(IELB)) THEN
                     DO 200 J = IXELB(IELB-1)+1, IXELB(IELB)
                        IF (PARMC(J) .EQ. 0.0) THEN
                           RESULT(J) = PARMT(J)
                        ELSE
                           RESULT(J) = PARMF(J)
                        END IF
  200                CONTINUE
                  END IF
  210          CONTINUE
            END IF
         END IF

      ELSE IF (TEST .EQ. 'IFGZ') THEN

         IF (NUMIX .GE. 0) THEN
            IF (PARTYP .NE. 'E') THEN
               DO 220 I = 1, NUMIX
                  J = IXNODE(I)
                  IF (PARMC(J) .GT. 0.0) THEN
                     RESULT(J) = PARMT(J)
                  ELSE
                     RESULT(J) = PARMF(J)
                  END IF
  220          CONTINUE
            ELSE
               DO 240 IELB = 1, NELBLK
                  IF (ISEVOK(IELB)) THEN
                     DO 230 I = IXELB(IELB-1)+1, IXELB(IELB)
                        J = IXELEM(I)
                        IF (PARMC(J) .GT. 0.0) THEN
                           RESULT(J) = PARMT(J)
                        ELSE
                           RESULT(J) = PARMF(J)
                        END IF
  230                CONTINUE
                  END IF
  240          CONTINUE
            END IF

         ELSE
C           NUMIX < 0; ALL VALUES SELECTED
C           NUMEL=NUMELO; NUMNP=NUMNPO; IXELB(J)=IXELBO(J), J=1,NELBLK
            IF (PARTYP .NE. 'E') THEN
               DO 250 J = 1, LENARY
                  IF (PARMC(J) .GT. 0.0) THEN
                     RESULT(J) = PARMT(J)
                  ELSE
                     RESULT(J) = PARMF(J)
                  END IF
  250          CONTINUE
            ELSE
               DO 270 IELB = 1, NELBLK
                  IF (ISEVOK(IELB)) THEN
                     DO 260 J = IXELB(IELB-1)+1, IXELB(IELB)
                        IF (PARMC(J) .GT. 0.0) THEN
                           RESULT(J) = PARMT(J)
                        ELSE
                           RESULT(J) = PARMF(J)
                        END IF
  260                CONTINUE
                  END IF
  270          CONTINUE
            END IF
         END IF
      END IF

      RETURN
      END
