C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CHKPAR (TEST, MSG,
     &   PARTYP, LENARY, NELBLK, IXELB, ISEVOK,
     &   NUMIX, IXNODE, IXELEM,
     &   PARM1, PARM2, *)
C=======================================================================

C   --*** CHKPAR *** (ALGEBRA) Check parameter data
C   --   Written by Amy Gilkey - revised 05/17/88
C   --
C   --CHKPAR checks the data for certain conditions (such as equality to
C   --zero).  If the data is in error, an error message is printed and
C   --the alternate return is taken.
C   --
C   --Parameters:
C   --   TEST - IN - the test to be performed:
C   --      EQ0 - if any data is equal to zero
C   --      LT0 - if any data is less than zero
C   --      LE0 - if any data is less than or equal to zero
C   --      ABSGT1 - if the absolute value of any data is greater than one
C   --      2EQ0 - if both data values equal zero
C   --   MSG - IN - the message if the condition is met
C   --   PARTYP - IN - the parameter type (element must be handled differently)
C   --   LENARY - IN - the length of the parameter arrays
C   --   NELBLK - IN - the number of element blocks
C   --   IXELB - IN - the cumulative element counts for each element block
C   --   ISEVOK - IN - the element variable truth table for each element block
C   --   NUMIX - IN - the number of selected values; <0 if all
C   --   IXNODE - IN - the indices of the selected nodes (only if NUMIX >= 0)
C   --   IXELEM - IN - the indices of the selected elements (only if NUMIX >= 0)
C   --   PARM1, PARM2 - IN - the data arrays to be checked
C   --   * - the return statement if the condition is met

      CHARACTER*(*) TEST
      CHARACTER*(*) MSG
      CHARACTER PARTYP
      INTEGER IXELB(0:NELBLK)
      LOGICAL ISEVOK(*)
      INTEGER IXNODE(*)
      INTEGER IXELEM(*)
      REAL PARM1(*), PARM2(*)

      IF (TEST .EQ. 'EQ0') THEN

         IF (NUMIX .GE. 0) THEN
            IF (PARTYP .NE. 'E') THEN
               DO 100 I = 1, NUMIX
                  J = IXNODE(I)
                  IF (PARM1(J) .EQ. 0) GOTO 400
  100          CONTINUE
            ELSE
               DO 120 IELB = 1, NELBLK
                  IF (ISEVOK(IELB)) THEN
                     DO 110 I = IXELB(IELB-1)+1, IXELB(IELB)
                        J = IXELEM(I)
                        IF (PARM1(J) .EQ. 0) GOTO 400
  110                CONTINUE
                  END IF
  120          CONTINUE
            END IF

         ELSE
            IF (PARTYP .NE. 'E') THEN
               DO 130 J = 1, LENARY
                  IF (PARM1(J) .EQ. 0) GOTO 400
  130          CONTINUE
            ELSE
               DO 150 IELB = 1, NELBLK
                  IF (ISEVOK(IELB)) THEN
                     DO 140 J = IXELB(IELB-1)+1, IXELB(IELB)
                        IF (PARM1(J) .EQ. 0) GOTO 400
  140                CONTINUE
                  END IF
  150          CONTINUE
            END IF
         END IF

      ELSE IF (TEST .EQ. 'LT0') THEN

         IF (NUMIX .GE. 0) THEN
            IF (PARTYP .NE. 'E') THEN
               DO 160 I = 1, NUMIX
                  J = IXNODE(I)
                  IF (PARM1(J) .LT. 0) GOTO 400
  160          CONTINUE
            ELSE
               DO 180 IELB = 1, NELBLK
                  IF (ISEVOK(IELB)) THEN
                     DO 170 I = IXELB(IELB-1)+1, IXELB(IELB)
                        J = IXELEM(I)
                        IF (PARM1(J) .LT. 0) GOTO 400
  170                CONTINUE
                  END IF
  180          CONTINUE
            END IF

         ELSE
            IF (PARTYP .NE. 'E') THEN
               DO 190 J = 1, LENARY
                  IF (PARM1(J) .LT. 0) GOTO 400
  190          CONTINUE
            ELSE
               DO 210 IELB = 1, NELBLK
                  IF (ISEVOK(IELB)) THEN
                     DO 200 J = IXELB(IELB-1)+1, IXELB(IELB)
                        IF (PARM1(J) .LT. 0) GOTO 400
  200                CONTINUE
                  END IF
  210          CONTINUE
            END IF
         END IF

      ELSE IF (TEST .EQ. 'LE0') THEN

         IF (NUMIX .GE. 0) THEN
            IF (PARTYP .NE. 'E') THEN
               DO 220 I = 1, NUMIX
                  J = IXNODE(I)
                  IF (PARM1(J) .LE. 0) GOTO 400
  220          CONTINUE
            ELSE
               DO 240 IELB = 1, NELBLK
                  IF (ISEVOK(IELB)) THEN
                     DO 230 I = IXELB(IELB-1)+1, IXELB(IELB)
                        J = IXELEM(I)
                        IF (PARM1(J) .LE. 0) GOTO 400
  230                CONTINUE
                  END IF
  240          CONTINUE
            END IF

         ELSE
            IF (PARTYP .NE. 'E') THEN
               DO 250 J = 1, LENARY
                  IF (PARM1(J) .LE. 0) GOTO 400
  250          CONTINUE
            ELSE
               DO 270 IELB = 1, NELBLK
                  IF (ISEVOK(IELB)) THEN
                     DO 260 J = IXELB(IELB-1)+1, IXELB(IELB)
                        IF (PARM1(J) .LE. 0) GOTO 400
  260                CONTINUE
                  END IF
  270          CONTINUE
            END IF
         END IF

      ELSE IF (TEST .EQ. 'ABSGT1') THEN

         IF (NUMIX .GE. 0) THEN
            IF (PARTYP .NE. 'E') THEN
               DO 280 I = 1, NUMIX
                  J = IXNODE(I)
                  IF (ABS(PARM1(J)) .GE. 0) GOTO 400
  280          CONTINUE
            ELSE
               DO 300 IELB = 1, NELBLK
                  IF (ISEVOK(IELB)) THEN
                     DO 290 I = IXELB(IELB-1)+1, IXELB(IELB)
                        J = IXELEM(I)
                        IF (ABS(PARM1(J)) .GE. 0) GOTO 400
  290                CONTINUE
                  END IF
  300          CONTINUE
            END IF

         ELSE
            IF (PARTYP .NE. 'E') THEN
               DO 310 J = 1, LENARY
                  IF (ABS(PARM1(J)) .GE. 0) GOTO 400
  310          CONTINUE
            ELSE
               DO 330 IELB = 1, NELBLK
                  IF (ISEVOK(IELB)) THEN
                     DO 320 J = IXELB(IELB-1)+1, IXELB(IELB)
                        IF (ABS(PARM1(J)) .GE. 0) GOTO 400
  320                CONTINUE
                  END IF
  330          CONTINUE
            END IF
         END IF

      ELSE IF (TEST .EQ. '2EQ0') THEN

         IF (NUMIX .GE. 0) THEN
            IF (PARTYP .NE. 'E') THEN
               DO 340 I = 1, NUMIX
                  J = IXNODE(I)
                  IF ((PARM1(J) .EQ. 0) .AND. (PARM2(J) .EQ. 0))
     &               GOTO 400
  340          CONTINUE
            ELSE
               DO 360 IELB = 1, NELBLK
                  IF (ISEVOK(IELB)) THEN
                     DO 350 I = IXELB(IELB-1)+1, IXELB(IELB)
                        J = IXELEM(I)
                        IF ((PARM1(J) .EQ. 0) .AND. (PARM2(J) .EQ. 0))
     &                     GOTO 400
  350                CONTINUE
                  END IF
  360          CONTINUE
            END IF

         ELSE
            IF (PARTYP .NE. 'E') THEN
               DO 370 J = 1, LENARY
                  IF ((PARM1(J) .EQ. 0) .AND. (PARM2(J) .EQ. 0))
     &               GOTO 400
  370          CONTINUE
            ELSE
               DO 390 IELB = 1, NELBLK
                  IF (ISEVOK(IELB)) THEN
                     DO 380 J = IXELB(IELB-1)+1, IXELB(IELB)
                        IF ((PARM1(J) .EQ. 0) .AND. (PARM2(J) .EQ. 0))
     &                     GOTO 400
  380                CONTINUE
                  END IF
  390          CONTINUE
            END IF
         END IF
      END IF

      RETURN

  400 CONTINUE
      CALL PRTERR ('FATAL', MSG)
      RETURN 1
      END
