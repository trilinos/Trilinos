C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE SETN02 (MXND, NXL, LXK, KXL, LINE, NEND, NODE, N0, N2)
C***********************************************************************

C  SUBROUTINE SETN02 = PICKS THE NEXT LINE AROUND THE ELEMENTS ATTACHED
C                      TO LINE WITH ONE END AT NEND, AND THE OTHER END
C                      NOT AT NODE, AND FROM THE CONNECTIVITY OF THE
C                      ELEMENTS DETERMINES THE BOUNDING ANGULAR LINES
C                      AND NODES.

C***********************************************************************

      DIMENSION NXL(2, 3*MXND)
      DIMENSION LXK(4, MXND), KXL(2, 3*MXND)

      K1 = KXL (1, LINE)
      K2 = KXL (2, LINE)

C  FIND THE NEXT LINE IN K1

      DO 100 I = 1, 4
         IL = LXK (I, K1)
         IF ((NXL (1, IL) .EQ. NEND) .AND.
     &      (NXL (2, IL) .NE. NODE)) THEN
            L1 = IL
            NNEW1 = NXL (2, IL)
            GOTO 110
         ELSEIF ((NXL (2, IL) .EQ. NEND) .AND.
     &      (NXL (1, IL) .NE. NODE)) THEN
            L1 = IL
            NNEW1 = NXL (1, IL)
            GOTO 110
         ENDIF
  100 CONTINUE
      CALL MESAGE ('** PROBLEMS IN SETN02 FINDING NNEW1 **')
      RETURN

  110 CONTINUE

C  FIND THE NEXT LINE IN K2

      DO 120 I = 1, 4
         IL = LXK (I, K2)
         IF ((NXL (1, IL) .EQ. NEND) .AND.
     &      (NXL (2, IL) .NE. NODE)) THEN
            NNEW2 = NXL (2, IL)
            GOTO 130
         ELSEIF ((NXL (2, IL) .EQ. NEND) .AND.
     &      (NXL (1, IL) .NE. NODE)) THEN
            NNEW2 = NXL (1, IL)
            GOTO 130
         ENDIF
  120 CONTINUE
      CALL MESAGE ('** PROBLEMS IN SETN02 FINDING NNEW2 **')
      RETURN

  130 CONTINUE

C  NOW DETERMINE WHICH OF THESE NODES IS N0 AND WHICH IS N2 BASED
C  ON THE FACT THAT THE CONNECTIVITY OF THE ELEMENTS LINES IS ALWAYS IN
C  COUNTER-CLOCKWISE ORDER

      DO 140 I = 1, 4
         IF (LXK (I, K1) .EQ. LINE) THEN
            I0 = I - 1
            I2 = I + 1
            IF (I .EQ. 1) THEN
               I0 = 4
            ELSEIF (I .EQ. 4) THEN
               I2 = 1
            ENDIF
            L0 = LXK (I0, K1)
            L2 = LXK (I2, K1)
            IF (L0 .EQ. L1) THEN
               N0 = NNEW1
               N2 = NNEW2
            ELSEIF (L2 .EQ. L1) THEN
               N0 = NNEW2
               N2 = NNEW1
            ELSE
               CALL MESAGE ('** PROBLEMS IN SETN02 FINDING N0 '//
     &            'AND N2 **')
            ENDIF
            GOTO 150
         ENDIF
  140 CONTINUE
      CALL MESAGE ('** PROBLEMS IN SETN02 FINDING LINE AGAIN **')

  150 CONTINUE

      RETURN

      END
