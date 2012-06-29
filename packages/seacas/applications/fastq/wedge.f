C $Id: wedge.f,v 1.2 2004/01/21 05:18:40 gdsjaar Exp $
C $Log: wedge.f,v $
C Revision 1.2  2004/01/21 05:18:40  gdsjaar
C Initialized several variables identified by valgrind.
C
C Revision 1.1.1.1  1990/11/30 11:17:41  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:17:39  gdsjaar
c Initial revision
c 
C
CC* FILE: [.PAVING]WEDGE.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE WEDGE (MXND, MLN, NUID, LXK, KXL, NXL, LXN, XN, YN,
     &   LNODES, BNSIZE, IAVAIL, NAVAIL, LLL, KKK, NNN, LLLOLD, NNNOLD,
     &   N1, N6, NLOOP, PWEDGE, GRAPH, VIDEO, NOROOM, ERR)
C***********************************************************************
C
C  SUBROUTINE WEDGE = EXPANDS TWO SIDE LINES INTO A CORNER NODE
C                      THIS IS REFERRED TO AS A WEDGE.
C
C***********************************************************************
C
      DIMENSION NUID (MXND), XN (MXND), YN (MXND)
      DIMENSION LXK (4, MXND), KXL (2, 3*MXND)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)
      DIMENSION LNODES (MLN, MXND), BNSIZE (2, MXND)
      DIMENSION L1LIST(20)
C
      LOGICAL GRAPH, ERR, NOROOM, PWEDGE, VIDEO, TWOLIN
C
      ERR = .FALSE.
      NNNOLD = MIN (NNN, NNNOLD)
      IF (LXN (3, N1) .EQ. 0) THEN
         TWOLIN = .TRUE.
      ELSE
         TWOLIN = .FALSE.
      ENDIF
C
C  GET ALL THE DEFINITIONS IN ORDER
C
      l3 = 0
      N0 = LNODES (2, N1)
      N2 = LNODES (3, N1)
      L1 = LNODES (5, N0)
      L2 = LNODES (5, N1)
      KL1 = KXL (1, L1)
C
C  FIND L3, N4, AND KL3
C
      IF (.NOT. TWOLIN) THEN
         DO 100 I = 1, 4
            LTEST = LXK (I, KL1)
            IF (LTEST .NE. L1) THEN
               IF (NXL (1, LTEST) .EQ. N1) THEN
                  L3 = LTEST
                  N4 = NXL (2, LTEST)
                  GOTO 110
               ELSEIF (NXL (2, LTEST) .EQ. N1) THEN
                  L3 = LTEST
                  N4 = NXL (1, LTEST)
                  GOTO 110
               ENDIF
            ENDIF
  100    CONTINUE
         CALL MESAGE ('** PROBLEMS IN WEDGE FINDING L3 AND N4 **')
         ERR = .TRUE.
         GOTO 200
  110    CONTINUE
C
C  FIND THE ELEMENT KL3 - THE ELEMENT ON THE OTHER SIDE OF L3
C
         IF (KXL (1, L3) .EQ. KL1) THEN
            KL3 = KXL (2, L3)
         ELSEIF (KXL (2, L3) .EQ. KL1) THEN
            KL3 = KXL (1, L3)
         ELSE
            CALL MESAGE ('** PROBLEMS IN WEDGE FINDING KL3 **')
            ERR = .TRUE.
            GOTO 200
         ENDIF
      ENDIF
C
C  GET ALL THE N1 LINES
C
      CALL GETLXN (MXND, LXN, N1, L1LIST, NL, ERR)
      IF (ERR) THEN
         CALL MESAGE ('** PROBLEMS IN WEDGE GETTING N1 LINES **')
         GOTO 200
      ENDIF
C
C  ERASE THE LINES
C
      IF ((GRAPH) .OR. (VIDEO)) THEN
         CALL LCOLOR ('BLACK')
         DO 120 I = 1, NL
            NODE1 = NXL (1, L1LIST (I))
            NODE2 = NXL (2, L1LIST (I))
            CALL D2NODE (MXND, XN, YN, NODE1, NODE2)
  120    CONTINUE
         IF (GRAPH) THEN
            CALL LCOLOR ('WHITE')
         ELSE
            CALL LCOLOR ('YELOW')
         ENDIF
         CALL SFLUSH
      ENDIF
C
C  NOW THAT ALL THE NECESSARY VARAIBLES HAVE BEEN DEFINED,
C  START BY DEFINING THE LOCATION OF ALL THE NEW NODES
C
      NNN = NNN + 2
      IF (NNN .GT. MXND) THEN
         NOROOM = .TRUE.
         GOTO 200
      ENDIF
      N5 = NNN - 1
      N6 = NNN
      IF (PWEDGE) THEN
         XN (N5) = XN (N1)
         YN (N5) = YN (N1)
         XN (N1) = XN (N0) + ( (XN (N5) - XN (N0)) * .33333 )
         YN (N1) = YN (N0) + ( (YN (N5) - YN (N0)) * .33333 )
         XN (N6) = XN (N0) + ( (XN (N5) - XN (N0)) * .66667 )
         YN (N6) = YN (N0) + ( (YN (N5) - YN (N0)) * .66667 )
      ELSE
         XN (N6) = XN (N0) + XN (N2) - XN (N4)
         YN (N6) = YN (N0) + YN (N2) - YN (N4)
         XN (N5) = ( (XN (N1) + ( ( XN (N2) - XN (N1) ) / 3. ) ) +
     &      ((XN (N6) + XN (N2)) * .5) ) * .5
         YN (N5) = ( (YN (N1) + ( ( YN (N2) - YN (N1) ) / 3. ) ) +
     &      ((YN (N6) + YN (N2)) * .5) ) * .5
         XN (N1) = ( (XN (N0) + ( ( XN (N1) - XN (N0) ) * .666667 ) ) +
     &      ((XN (N6) + XN (N0)) * .5) ) * .5
         YN (N1) = ( (YN (N0) + ( ( YN (N1) - YN (N0) ) * .666667 ) ) +
     &      ((YN (N6) + YN (N0)) * .5) ) * .5
      ENDIF
C
C  NOW ADD LINES L4, L5, AND L6
C
      LLL = LLL + 3
      L4 = LLL - 2
      IF (TWOLIN) THEN
         NXL (1, L4) = N1
         NXL (2, L4) = N2
      ELSE
         NXL (1, L4) = N4
         NXL (2, L4) = N5
      ENDIF
      L5 = LLL - 1
      L6 = LLL
      NXL (1, L5) = N1
      NXL (2, L5) = N6
      NXL (1, L6) = N5
      NXL (2, L6) = N6
C
C  NOW UPDATE THE LXN ARRAYS
C
      DO 130 I = 1, 4
         LXN (I, N5) = 0
         LXN (I, N6) = 0
  130 CONTINUE
C
      CALL FIXLXN (MXND, LXN, NXL, NUID, NAVAIL, IAVAIL, NNN, LLL,
     &   NNNOLD, LLLOLD, NOROOM, ERR)
      LLLOLD = LLL
      NNNOLD = NNN
      IF ((NOROOM) .OR. (ERR)) GOTO 200
C
C  UPDATE EXISTING NODES AND THEIR LXN ARRAYS
C  - FIRST MAKE SURE THAT N5 TAKES N1'S BOUNDARY STATUS IF PWEDGE IS ON
C
      IF ((PWEDGE) .AND. (LXN (2, N1) .LT. 0))
     &   LXN (2, N5) = - IABS (LXN (2, N5))
C
C  MAKE SURE THAT ALL THE LINES THAT WERE ATTACHED TO N1
C  EXCLUDING L1, L3, (STILL HOOKED TO N1) L4 AND L6 (ALREADY
C  HOOKED TO N5) ARE NOW ATTACHED TO N5 AND HAVE N5 AS THE
C  CORRECT ENDPOINT
C
      DO 140 I = 1, NL
         LL = L1LIST (I)
         IF ((LL .NE. L1) .AND. (LL .NE. L3) .AND. (LL .NE. L4) .AND.
     &      (LL .NE. L6)) THEN
C
            CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, N1,
     &         LL, NNN, ERR, NOROOM)
            IF ((NOROOM) .OR. (ERR)) THEN
               CALL MESAGE ('** PROBLEMS IN WEDGE UNHOOKING LL FROM '//
     &            'N1 **')
               GOTO 200
            ENDIF
C
            CALL ADDLXN (MXND, LXN, NUID, NAVAIL, IAVAIL,
     &         N5, LL, NNN, ERR, NOROOM)
            IF ((NOROOM) .OR. (ERR)) THEN
               CALL MESAGE ('** PROBLEMS IN WEDGE HOOKING LL TO N5 **')
               GOTO 200
            ENDIF
C
            IF (NXL (1, LL) .EQ. N1) THEN
               NXL (1, LL) = N5
            ELSEIF (NXL (2, LL) .EQ. N1) THEN
               NXL (2, LL) = N5
            ELSE
               CALL MESAGE ('** PROBLEMS IN WEDGE CHANGING LL '//
     &            'ENDPOINT FROM N1 TO N5 **')
               ERR = .TRUE.
               GOTO 200
            ENDIF
C
         ENDIF
  140 CONTINUE
C
C  NOW ADD THE NEW ELEMENT
C
      KKK = KKK + 1
      IF (TWOLIN) THEN
         LXK (1, KKK) = L4
         LXK (2, KKK) = L2
      ELSE
         LXK (1, KKK) = L3
         LXK (2, KKK) = L4
      ENDIF
      LXK (3, KKK) = L6
      LXK (4, KKK) = L5
C
C  NOW FIX THE KXL ARRAY FOR LINE L3
C
      IF (.NOT. TWOLIN) THEN
         IF (KXL (1, L3) .EQ. KL3) THEN
            KXL (1, L3) = KKK
         ELSEIF (KXL (2, L3) .EQ. KL3) THEN
            KXL (2, L3) = KKK
         ELSE
            CALL MESAGE ('** PROBLEMS IN WEDGE REPLACING KL3 FOR L3 **')
            ERR = .TRUE.
            GOTO 200
         ENDIF
      ENDIF
C
C  ADD THE KXL ENTRIES FOR THE NEW LINES
C
      IF (TWOLIN) THEN
         KXL (1, L4) = KL1
         KXL (2, L4) = KKK
         KXL (1, L2) = KKK
      ELSE
         KXL (1, L4) = KKK
         KXL (2, L4) = KL3
      ENDIF
      KXL (1, L5) = KKK
      KXL (1, L6) = KKK
C
C  NOW FIX THE LXK ARRAY FOR THE ELEMENT KL1 IF TWOLIN
C
      IF (TWOLIN) THEN
         DO 150 I = 1, 4
            IF (LXK (I, KL1) .EQ. L2) THEN
               LXK (I, KL1) = L4
               GOTO 160
            ENDIF
  150    CONTINUE
         CALL MESAGE ('** PROBLEMS IN WEDGE REPLACING L2 WITH L4 IN '//
     &      'KL1 **')
         ERR = .TRUE.
         GOTO 200
  160    CONTINUE
C
C  OTHERWISE FIX THE LXK ARRAY FOR THE ELEMENT KL3
C
      ELSE
         DO 170 I = 1, 4
            IF (LXK (I, KL3) .EQ. L3) THEN
               LXK (I, KL3) = L4
               GOTO 180
            ENDIF
  170    CONTINUE
         CALL MESAGE ('** PROBLEMS IN WEDGE REPLACING L3 WITH L4 IN '//
     &      'KL3 **')
         ERR = .TRUE.
         GOTO 200
  180    CONTINUE
      ENDIF
C
C  NOW FIX THE LNODES ARRAY
C
      LNODES (1, N5) = 0
      LNODES (1, N6) = 0
C
      LNODES (2, N6) = N1
      LNODES (2, N5) = N6
      LNODES (2, N2) = N5
C
      LNODES (3, N1) = N6
      LNODES (3, N6) = N5
      LNODES (3, N5) = N2
C
      LNODES (4, N5) = - 1
      LNODES (4, N6) = - 1
C
      LNODES (5, N1) = L5
      LNODES (5, N6) = L6
      LNODES (5, N5) = L2
C
      LNODES (8, N5) = LNODES (8, N1)
      LNODES (8, N6) = LNODES (8, N1)
C
C  PUT THE BEGINNING BOUNDARY DISTANCE IN PLACE
C
      BNSIZE (1, N5) = BNSIZE (1, N1)
      BNSIZE (1, N6) = BNSIZE (1, N1)
      IF (BNSIZE (1, N1) .GT. 0.) THEN
         SIZMIN = AMIN1 (BNSIZE (1, N1) * BNSIZE (2, N1),
     &      BNSIZE (1, N0) * BNSIZE (2, N0),
     &      BNSIZE (1, N2) * BNSIZE (2, N2)) / BNSIZE (1, N1)
      ELSE
         SIZMIN = AMIN1 (BNSIZE (1, N1) * BNSIZE (2, N1),
     &      BNSIZE (1, N0) * BNSIZE (2, N0),
     &      BNSIZE (1, N2) * BNSIZE (2, N2)) /
     &      SQRT ( (XN (N1) - XN (N2)) ** 2 + (YN (N1) - YN (N2)) ** 2 )
      ENDIF
      IF (PWEDGE) THEN
         BNSIZE (2, N6) = SIZMIN
         BNSIZE (2, N5) = SIZMIN
         BNSIZE (2, N1) = SIZMIN
      ELSE
         BNSIZE (2, N5) = BNSIZE (2, N1) * 1.15
         BNSIZE (2, N6) = BNSIZE (2, N1)
         BNSIZE (2, N1) = BNSIZE (2, N1) * 1.15
      ENDIF
C
      NLOOP = NLOOP + 2
      ERR = .FALSE.
C
C  NOW REDRAW THE LINES
C
      IF ((GRAPH) .OR. (VIDEO)) THEN
         DO 190 I = 1, NL
            IF ((L1LIST (I) .NE. L1) .AND. (L1LIST (I) .NE. L3) .AND.
     &         (L1LIST (I) .NE. L2)) THEN
               NODE1 = NXL (1, L1LIST (I))
               NODE2 = NXL (2, L1LIST (I))
               CALL D2NODE (MXND, XN, YN, NODE1, NODE2)
            ENDIF
  190    CONTINUE
         CALL D2NODE (MXND, XN, YN, N0, N1)
         CALL D2NODE (MXND, XN, YN, N1, N6)
         CALL D2NODE (MXND, XN, YN, N2, N5)
         CALL D2NODE (MXND, XN, YN, N6, N5)
         IF (TWOLIN) THEN
            CALL D2NODE (MXND, XN, YN, N1, N2)
         ELSE
            CALL D2NODE (MXND, XN, YN, N1, N4)
            CALL D2NODE (MXND, XN, YN, N4, N5)
         ENDIF
         CALL SFLUSH
      ENDIF
C
  200 CONTINUE
C
      RETURN
C
      END
