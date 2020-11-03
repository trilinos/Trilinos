C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES, NODE,
     &   ERR)
C***********************************************************************

C  SUBROUTINE MARKSM = MARKS NODES WITHIN 2 LINE CONNECTIONS FROM NODE
C                      FOR SMOOTHING

C***********************************************************************

      DIMENSION LXK (4, MXND), KXL (2, 3*MXND)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)
      DIMENSION LNODES (MLN, MXND)
      DIMENSION L1LIST(20), L2LIST(20)

      LOGICAL ERR

      IF (LXN (1, NODE) .LE. 0) GOTO 120
      CALL GETLXN (MXND, LXN, NODE, L1LIST, NL1, ERR)
      IF (ERR) THEN
         CALL MESAGE ('** PROBLEMS IN MARKSM FINDING LXN **')
         GOTO 120
      ENDIF

      LNODES (4, NODE) = - IABS (LNODES (4, NODE))
      DO 110 I = 1, NL1
         NODE2 = NXL (1, L1LIST (I)) + NXL (2, L1LIST (I)) - NODE
         CALL GETLXN (MXND, LXN, NODE2, L2LIST, NL2, ERR)
         IF (ERR) THEN
            CALL MESAGE ('** PROBLEMS IN MARKSM FINDING LXN **')
            GOTO 120
         ENDIF
         LNODES (4, NODE2) = - IABS (LNODES (4, NODE2))
         DO 100 J = 1, NL2
            NODE1 = NXL (1, L2LIST (J)) + NXL (2, L2LIST (J)) - NODE2
            LNODES (4, NODE1) = - IABS (LNODES (4, NODE1))
  100    CONTINUE
  110 CONTINUE

  120 CONTINUE
      RETURN

      END
