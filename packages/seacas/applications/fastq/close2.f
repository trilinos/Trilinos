C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE CLOSE2 (MXND, MLN, NUID, XN, YN, ZN, LXK, KXL, NXL,
     &   LXN, LNODES, IAVAIL, NAVAIL, NNN, LLL, N1, XMIN, XMAX, YMIN,
     &   YMAX, ZMIN, ZMAX, PGRAPH, VIDEO, DEV1, KREG, NOROOM, ERR)
C***********************************************************************

C  SUBROUTINE CLOSE2 = SEALS OFF THE LAST 2 OPEN LINES WHILE CHECKING
C                      FOR FORMING A 2-LINE NODE ON THE INTERIOR
C                      (A 2-LINE NODE GENERATES 2 DEGENERATE QUADS)

C***********************************************************************

      DIMENSION NUID (MXND), XN (MXND), YN (MXND), ZN (MXND)
      DIMENSION LXK (4, MXND), KXL (2, 3*MXND)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)
      DIMENSION LNODES (MLN, MXND)

      LOGICAL ERR, NOROOM, FOUND, PGRAPH, DONE, CHECK, VIDEO

      CHARACTER*3 DEV1

      ERR = .FALSE.
      CHECK = .FALSE.

      N0 = LNODES (2, N1)
      LINE1 = LNODES (5, N0)
      LINE2 = LNODES (5, N1)

C  CHECK TO MAKE SURE THAT AT LEAST ONE OF THE LINES
C  IS NOT A BOUNDARY LINE AND GET THE NODE TO BE DELETED

  100 CONTINUE
      IF ((KXL (1, LINE1) .GT. 0) .OR.
     &   (KXL (1, LINE2) .GT. 0)) THEN

         FOUND = .TRUE.

         IF (KXL (1, LINE1) .GT. 0) THEN
            LNEW = LINE2
            LOLD = LINE1
         ELSE
            LNEW = LINE1
            LOLD = LINE2
         ENDIF
         KOLD = KXL (1, LOLD)
         KNEW = KXL (1, LNEW)

C  CHECK FOR ONE OF THE NODES BEING A TWO LINE NODE

         IF (KOLD. EQ. KNEW) THEN
            IF (LXN (3, N0) .EQ. 0) THEN
               NGONE = N0
               NTHERE = N1
            ELSEIF (LXN (3, N1) .EQ. 0) THEN
               NGONE = N1
               NTHERE = N0
            ELSE
               CALL MESAGE ('** PROBLEMS WITH NO TWO LINE NODE'//
     &            ' ATTACHED IN CLOSE2 **')
               ERR = .TRUE.
               GOTO 150
            ENDIF

C  DELETE THE TWO-LINE NODE, THE TWO LINES, AND THE ELEMENT

            KXL (1, LOLD) = 0
            KXL (2, LOLD) = 0
            NXL (1, LOLD) = 0
            NXL (2, LOLD) = 0
            KXL (1, LNEW) = 0
            KXL (2, LNEW) = 0
            NXL (1, LNEW) = 0
            NXL (2, LNEW) = 0

C  UNHOOK BOTH LINES FROM NTHERE

            CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, NTHERE,
     &         LOLD, NNN, ERR, NOROOM)
            IF (ERR) THEN
               CALL MESAGE ('** PROBLEMS IN CLOSE2 DELETING LOLD'//
     &            ' FROM NTHERE **')
               GOTO 150
            ENDIF
            CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, NTHERE,
     &         LNEW, NNN, ERR, NOROOM)
            IF (ERR) THEN
               CALL MESAGE ('** PROBLEMS IN CLOSE2 DELETING LNEW'//
     &            ' FROM NTHERE **')
               GOTO 150
            ENDIF

C  NOW DELETE NGONE AND THE ELEMENT

            DO 110 I = 1, 4
               LXN (I, NGONE) = 0
               IF ( (LXK (I, KOLD) .EQ. LNEW) .OR.
     &            (LXK (I, KOLD) .EQ. LOLD) ) LXK (I, KOLD) = 0
  110       CONTINUE
            LOLD = 0
            LNEW = 0
            DO 120 I = 1, 4
               IF (LXK (I, KOLD) .NE. 0) THEN
                  IF (LOLD .EQ. 0) THEN
                     LOLD = LXK (I, KOLD)
                  ELSE
                     LNEW = LXK (I, KOLD)
                  ENDIF
                  LXK (I, KOLD) = 0
               ENDIF
  120       CONTINUE
            KXL (1, LNEW) = KXL (1, LNEW) + KXL (2, LNEW) - KOLD
            KXL (2, LNEW) = 0
            KXL (1, LOLD) = KXL (1, LOLD) + KXL (2, LOLD) - KOLD
            KXL (2, LOLD) = 0

C  NOW RESET THE NECESSARY VARIABLES

            N1 = NXL (1, LNEW)
            N0 = NXL (2, LNEW)
            LINE1 = LOLD
            LINE2 = LNEW
            GOTO 100
         ENDIF

C  DELETE THE OLD LINE AND REDO LINK ARRAYS

         IF (KNEW .EQ. 0) THEN
            KXL (1, LNEW) = KOLD
            KXL (2, LNEW) = 0
         ELSE
            KXL (1, LNEW) = KNEW
            KXL (2, LNEW) = KOLD
         ENDIF
         KXL (1, LOLD) = 0
         KXL (2, LOLD) = 0
         NXL (1, LOLD) = 0
         NXL (2, LOLD) = 0

C  FIX THE LINES PER ELEMENT ARRAY FOR THE ONE ELEMENT CHANGING

         DO 130 II = 1, 4
            IF (LXK (II, KOLD) .EQ. LOLD) THEN
               LXK (II, KOLD) = LNEW
               GOTO 140
            ENDIF
  130    CONTINUE
         CALL MESAGE ('** PROBLEMS IN CLOSE2 WITH CHANGING ELEMENT **')
         ERR = .TRUE.
         GOTO 150
  140    CONTINUE

C  FIX LXN ARRAY
C  UNHOOK LOLD FROM N0 AND FROM N1

         CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, N0,
     &      LOLD, NNN, ERR, NOROOM)
         IF (ERR) THEN
            CALL MESAGE ('** PROBLEMS IN CLOSE2 DELETING NNN LINES **')
            GOTO 150
         ENDIF
         CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, N1,
     &      LOLD, NNN, ERR, NOROOM)
         IF (ERR) THEN
            CALL MESAGE ('** PROBLEMS IN CLOSE2 DELETING N1 LINES **')
            GOTO 150
         ENDIF

C NOW FIX THE LNODES ARRAY

         LNODES (4, N1) = - 2
         LNODES (4, N0) = - 2

      ELSE
         CALL MESAGE ('** PINCHED TOO FAR IN CLOSE2 **')
         GOTO 150
      ENDIF

C  NOW SEE IF THE CLOSURE HAS PRODUCED A 2-LINE NODE AND
C  THUS REQUIRES THAT ONE OF THE ELEMENTS MUST BE SQUASHED

      IF ((LXN (3, N0) .EQ. 0) .AND. (LXN (2, N0) .GT. 0)) THEN
         CALL DELEM (MXND, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &      NNN, NAVAIL, IAVAIL, N0, KXL (1, LNEW), IDUM1, IDUM2,
     &      DONE, CHECK, NOROOM, ERR)
         IF ((NOROOM) .OR. (ERR)) GOTO 150
      ELSEIF ((LXN (3, N1) .EQ. 0) .AND. (LXN (2, N1) .GT. 0)) THEN
         CALL DELEM (MXND, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &      NNN, NAVAIL, IAVAIL, N1, KXL (1, LNEW), IDUM1, IDUM2,
     &      DONE, CHECK, NOROOM, ERR)
         IF ((NOROOM) .OR. (ERR)) GOTO 150
      ENDIF

      IF ( (FOUND) .AND. ((PGRAPH) .OR. (VIDEO)) ) THEN
         CALL RPLOTL (MXND, XN, YN, ZN, NXL, XMIN, XMAX, YMIN, YMAX,
     &      ZMIN, ZMAX, LLL, DEV1, KREG)
         IF (VIDEO) CALL SNAPIT (1)
      ENDIF

  150 CONTINUE

      RETURN

      END
