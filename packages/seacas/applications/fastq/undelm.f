C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE UNDELM (MXND, MLN, LNODES, XN, YN, NUID, LXK, KXL, NXL,
     &   LXN, NNN, LLL, KKK, NAVAIL, IAVAIL, N0, N1, N2, N3, L1, L2, L3,
     &   K1, K2, NOROOM, ERR, GRAPH, VIDEO)
C***********************************************************************

C  SUBROUTINE UNDELM = UNDELETES AN ELEMENT BY EXPANDING N1 INTO A
C                      NEW ELEMENT

C***********************************************************************

      DIMENSION LXK(4, MXND), NXL(2, 3*MXND), KXL(2, 3*MXND)
      DIMENSION LXN(4, MXND), XN(MXND), YN(MXND), NUID(MXND)
      DIMENSION LNODES (MLN, MXND)

      LOGICAL ERR, NOROOM, GRAPH, VIDEO

      ERR = .FALSE.

C  MAKE SURE THAT N2 HAS AT LEAST FOUR LINES ATTACHED TO IT

      IF (LXN (4, N2) .EQ. 0) THEN
         ERR = .TRUE.
         CALL MESAGE ('** N2 IN UNDELM CANNOT BE USED'//
     &      ' TO EXPAND AN ELEMENT **')
         GOTO 140
      ENDIF

C  ERASE THE LINE L3

      IF ((GRAPH) .OR. (VIDEO)) THEN
         CALL LCOLOR ('BLACK')
         CALL D2NODE (MXND, XN, YN, N0, N2)
         IF (GRAPH) THEN
            CALL LCOLOR ('WHITE')
         ELSE
            CALL LCOLOR ('YELOW')
         ENDIF
         CALL SFLUSH
         IF (VIDEO) CALL SNAPIT (3)
      ENDIF

C  DEFINE THE NEW NODE AND THE TWO NEW LINES

      NNN = NNN + 1
      IF (NNN .GT. MXND) THEN
         NOROOM = .TRUE.
         GOTO 140
      ENDIF
      NNEW = NNN
      XN (NNEW) = (XN (N0) + XN (N2)) * .5
      YN (NNEW) = (YN (N0) + YN (N2)) * .5

      LLL = LLL + 2
      L4 = LLL -1
      L5 = LLL
      NXL (1, L4) = N1
      NXL (2, L4) = NNEW
      NXL (1, L5) = N3
      NXL (2, L5) = NNEW

C  NOW CHANGE LINE L3'S END POINT FROM N2 TO NNEW

      IF (NXL (1, L3) .EQ. N2) THEN
         NXL (1, L3) = NNEW
      ELSEIF (NXL (2, L3) .EQ. N2) THEN
         NXL (2, L3) = NNEW
      ELSE
         CALL MESAGE ('** PROBLEMS IN UNDLEM WITH L3''S END POINT **')
         ERR = .TRUE.
         GOTO 140
      ENDIF

C  NOW UPDATE THE LXN ARRAYS

      LXN (1, NNEW) = L3
      LXN (2, NNEW) = L4
      LXN (3, NNEW) = L5
      LXN (4, NNEW) = 0
      CALL FIXLXN (MXND, LXN, NXL, NUID, NAVAIL, IAVAIL, NNN, LLL,
     &   NNN, LLL, ERR, NOROOM)
      IF ((NOROOM) .OR. (ERR)) GOTO 140

C  REMOVE L3 FROM THE LIST OF LINES FOR N2

      CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, N2,
     &   L3, NNN, ERR, NOROOM)
      IF ((NOROOM) .OR. (ERR)) THEN
         CALL MESAGE ('** PROBLEMS IN UNDELM UNHOOKING L3 FROM N2 **')
         GOTO 140
      ENDIF

C  ADD LINE L4 TO N1

      CALL ADDLXN (MXND, LXN, NUID, NAVAIL, IAVAIL,
     &   N1, L4, NNN, ERR, NOROOM)
      IF ((NOROOM) .OR. (ERR)) THEN
         CALL MESAGE ('** PROBLEMS IN UNDELM HOOKING L4 TO N1 **')
         GOTO 140
      ENDIF

C  ADD LINE L5 TO N3

      CALL ADDLXN (MXND, LXN, NUID, NAVAIL, IAVAIL,
     &   N3, L5, NNN, ERR, NOROOM)
      IF ((NOROOM) .OR. (ERR)) THEN
         CALL MESAGE ('** PROBLEMS IN UNDELM HOOKING L5 TO N3 **')
         GOTO 140
      ENDIF

C  NOW ADD THE NEW ELEMENT

      KKK = KKK + 1
      LXK (1, KKK) = L1
      LXK (2, KKK) = L2
      LXK (3, KKK) = L5
      LXK (4, KKK) = L4

C  NOW FIX THE KXL ARRAY FOR LINE L1

      IF (KXL (1, L1) .EQ. K2) THEN
         KXL (1, L1) = KKK
      ELSEIF (KXL (2, L1) .EQ. K2) THEN
         KXL (2, L1) = KKK
      ELSE
         CALL MESAGE ('** PROBLEMS IN UNDELM REPLACING K2 FOR L1 **')
         ERR = .TRUE.
         GOTO 140
      ENDIF

C  NOW FIX THE KXL ARRAY FOR LINE L2

      IF (KXL (1, L2) .EQ. K1) THEN
         KXL (1, L2) = KKK
      ELSEIF (KXL (2, L2) .EQ. K1) THEN
         KXL (2, L2) = KKK
      ELSE
         CALL MESAGE ('** PROBLEMS IN UNDELM REPLACING K1 FOR L2 **')
         ERR = .TRUE.
         GOTO 140
      ENDIF

C  ADD THE KXL ENTRIES FOR THE NEW LINES

      KXL (1, L4) = K2
      KXL (2, L4) = KKK
      KXL (1, L5) = K1
      KXL (2, L5) = KKK

C  NOW FIX THE LXK ARRAY FOR THE ELEMENT K1

      DO 100 I = 1, 4
         IF (LXK (I, K1) .EQ. L2) THEN
            LXK (I, K1) = L5
            GOTO 110
         ENDIF
  100 CONTINUE
      CALL MESAGE ('** PROBLEMS IN UNDELM REPLACING L2 WITH L5 IN '//
     &   'K1 **')
      ERR = .TRUE.
      GOTO 140
  110 CONTINUE

C  NOW FIX THE LXK ARRAY FOR THE ELEMENT K2

      DO 120 I = 1, 4
         IF (LXK (I, K2) .EQ. L1) THEN
            LXK (I, K2) = L4
            GOTO 130
         ENDIF
  120 CONTINUE
      CALL MESAGE ('** PROBLEMS IN UNDELM REPLACING L1 WITH L4 IN '//
     &   'K2 **')
      ERR = .TRUE.
      GOTO 140
  130 CONTINUE

C  NOW REDRAW THE LINES

      IF ((GRAPH) .OR. (VIDEO)) THEN
         CALL D2NODE (MXND, XN, YN, N0, NNEW)
         CALL D2NODE (MXND, XN, YN, N1, NNEW)
         CALL D2NODE (MXND, XN, YN, N3, NNEW)
         CALL SFLUSH
         IF (VIDEO) CALL SNAPIT (3)
      ENDIF

      LNODES (4, NNEW) = 2
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   NNEW, ERR)
      IF (ERR) GOTO 140
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   N0, ERR)
      IF (ERR) GOTO 140
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   N1, ERR)
      IF (ERR) GOTO 140
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   N2, ERR)
      IF (ERR) GOTO 140
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   N3, ERR)
      IF (ERR) GOTO 140
  140 CONTINUE
      RETURN

      END
