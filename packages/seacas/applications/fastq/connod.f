C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE CONNOD (MXND, MLN, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &   ANGLE, LNODES, NNN, KKK, LLL, NNNOLD, LLLOLD, N0, N1, N2,
     &   NODE, I, NLOOP, IAVAIL, NAVAIL, GRAPH, VIDEO, NOROOM, ERR)
C***********************************************************************

C  SUBROUTINE CONNOD = ADDS A NEW ELEMENT TO AN EXISTING NODE

C***********************************************************************

      DIMENSION XN (MXND), YN (MXND), NUID (MXND)
      DIMENSION LXK (4, MXND), KXL (2, 3*MXND)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)
      DIMENSION ANGLE (MXND), LNODES (MLN, MXND)

      LOGICAL GRAPH, VIDEO, ERR, NOROOM

      ERR = .FALSE.

      NLOOP = NLOOP-2

C  TAKE CARE OF A COUNTERCLOCKWISE CONNECTION

      IF (LNODES (3, N2) .EQ. NODE) THEN
         LNODES (3, N0) = NODE
         LNODES (2, NODE) = N0
         LNODES (4, N1) = - 2
         LNODES (4, N2) = - 2
         ANGLE (N1) = 0.
         ANGLE (N2) = 0.
         LLL = LLL+1
         KKK = KKK+1
         NXL (1, LLL) = NODE
         NXL (2, LLL) = N0

C  MAKE THE NEW ELEMENT

         LXK (1, KKK) = LNODES (5, N0)
         LXK (2, KKK) = LNODES (5, N1)
         LXK (3, KKK) = LNODES (5, N2)
         LXK (4, KKK) = LLL

         CALL ADDKXL (MXND, KXL, KKK, LLL)
         CALL ADDKXL (MXND, KXL, KKK, LNODES (5, N0))
         CALL ADDKXL (MXND, KXL, KKK, LNODES (5, N1))
         CALL ADDKXL (MXND, KXL, KKK, LNODES (5, N2))

         LNODES (5, N0) = LLL
         IF ((GRAPH) .OR. (VIDEO)) THEN
            CALL D2NODE (MXND, XN, YN, N0, NODE)
            CALL SFLUSH
            IF (VIDEO) CALL SNAPIT (1)
         ENDIF
         I = NODE

C  TAKE CARE OF A CLOCKWISE CONNECTION

      ELSE IF (LNODES (2, N0) .EQ. NODE) THEN
         LNODES (3, NODE) = N2
         LNODES (2, N2) = NODE
         LNODES (4, N0) = - 2
         LNODES (4, N1) = - 2
         ANGLE (N0) = 0.
         ANGLE (N1) = 0.
         LLL = LLL+1
         KKK = KKK+1
         NXL (1, LLL) = NODE
         NXL (2, LLL) = N2

C  MAKE THE NEW ELEMENT

         LXK (1, KKK) = LNODES (5, NODE)
         LXK (2, KKK) = LNODES (5, N0)
         LXK (3, KKK) = LNODES (5, N1)
         LXK (4, KKK) = LLL

         CALL ADDKXL (MXND, KXL, KKK, LNODES (5, NODE))
         CALL ADDKXL (MXND, KXL, KKK, LNODES (5, N0))
         CALL ADDKXL (MXND, KXL, KKK, LNODES (5, N1))
         CALL ADDKXL (MXND, KXL, KKK, LLL)

         LNODES (5, NODE) = LLL

         IF ((GRAPH) .OR. (VIDEO)) THEN
            CALL D2NODE (MXND, XN, YN, NODE, N2)
            CALL SFLUSH
            IF (VIDEO) CALL SNAPIT (1)
         ENDIF
         I = N2

      ELSE
         CALL MESAGE ('PROBLEM HANDLING CONNECTION IN CONNOD')
         CALL MESAGE ('CONNECTION DOES NOT FORM A SIMPLE CORNER')
         ERR = .TRUE.
         GOTO 100
      ENDIF
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES, N1, ERR)
      IF (ERR) GOTO 100
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES, N2, ERR)
      IF (ERR) GOTO 100
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES, NODE, ERR)
      IF (ERR) GOTO 100

      CALL FIXLXN (MXND, LXN, NXL, NUID, NAVAIL, IAVAIL, NNN, LLL,
     &   NNNOLD, LLLOLD, ERR, NOROOM)
      IF (ERR) GOTO 100
      LLLOLD = LLL
      NNNOLD = NNN

  100 CONTINUE

      RETURN

      END
