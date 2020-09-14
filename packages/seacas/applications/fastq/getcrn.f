C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE GETCRN (MXND, MXCORN, MLN, LNODES, NCORN, LCORN,
     &   ANGLE, XN, YN, LXN, NLOOP, N1, ONLYC, PPOSBL, GRAPH, ERR)
C***********************************************************************

C  SUBROUTINE GETCRN = SETS UP ALL THE POSSIBLE CORNER (OR NON-SIDE)
C                      LOCATIONS IN THE MESH

C***********************************************************************

      DIMENSION LNODES (MLN, MXND), LCORN (MXCORN), ANGLE (MXND)
      DIMENSION XN (MXND), YN (MXND), LXN (4, MXND), X(1), Y(1)

      LOGICAL ONLYC, GRAPH, PPOSBL, ERR

      ERR = .FALSE.

C  COUNT THE CURRENT "POSSIBLE NON-SIDES" STARTING AT THE I COUNTER

      NODE = N1
      NOLD = N1
      KOUNTC = 0
      KOUNT = 0
      NCORN = 0
      PPOSBL = .TRUE.
  100 CONTINUE
      KOUNT = KOUNT + 1
      IF (KOUNT .GT. NLOOP) THEN
         CALL MESAGE ('PROBLEM IN GETCRN WITH UNCLOSED LOOP')
         ERR = .TRUE.
         GOTO 120
      ENDIF

C  CHECK IF A PRIMITIVE IS EVEN POSSIBLE HERE

      IF ( LNODES (6, NODE) .GT. 4) PPOSBL = .FALSE.

C  CHECK FOR A POSSIBLE "NON - SIDE" NODE

      IF (ONLYC) CALL NDSTAT (NODE, LXN (1, NODE), ANGLE (NODE), ISTAT)

C  A NEW "POSSIBLE NON-SIDE" NODE HAS BEEN FOUND

      IF ( ( (ONLYC) .AND. (ISTAT .EQ. 1) ) .OR.
     &   ( (.NOT. ONLYC) .AND. (LNODES (6, NODE) .NE. 3) ) ) THEN
         IF (GRAPH) THEN
            ISQR = NODE
            X(1) = XN (ISQR)
            Y(1) = YN (ISQR)
            IF ((ONLYC) .AND. (ISTAT .EQ. 1)) THEN
               CALL SYMBOL (1, X, Y, 'DIAMND')
            ELSEIF (.NOT. ONLYC) THEN
               IF (LNODES (6, NODE) .EQ. 1) THEN
                  CALL SYMBOL (1, X, Y, 'DIAMND')
               ELSEIF (LNODES (6, NODE) .EQ. 2) THEN
                  CALL SYMBOL (1, X, Y, 'CENTRX')
               ELSEIF (LNODES (6, NODE) .EQ. 4) THEN
                  CALL SYMBOL (1, X, Y, 'TRINGL')
               ELSE
                  CALL SYMBOL (1, X, Y, 'CIRCLE')
               ENDIF
            ENDIF
            CALL SFLUSH
         ENDIF

C  TOO MANY "NON-SIDES" HAVE BEEN FOUND - JUST PICK ONE AND GO

         IF (NCORN .EQ. MXCORN) THEN
            PPOSBL = .FALSE.
            DO 110 J = 1, NCORN
               IF (LNODES (6, LCORN (J)) .EQ. 1) THEN
                  NODE = J
                  IF (GRAPH) CALL LCOLOR ('WHITE')
                  GOTO 120
               ENDIF
  110       CONTINUE
            NODE = LCORN (1)
            IF (GRAPH) CALL LCOLOR ('WHITE')
            GOTO 120
         ENDIF

C  INPUT THIS "NON-SIDE"

         NCORN = NCORN + 1
         LCORN (NCORN) = NODE

C  ADD UP THE NUMBER OF NODES FROM THE LAST "NON-SIDE"

         IF (NCORN .GT. 1) THEN
            LNODES (7, LCORN (NCORN - 1)) = KOUNTC + 1
         ELSE
            KKC = KOUNTC + 1
         ENDIF
         KOUNTC = 0

C  THIS IS A SIDE - JUST CONTINUE

      ELSE
         KOUNTC = KOUNTC + 1
      ENDIF

C  CHECK FOR COMPLETION OF THE LOOP

      NODE = LNODES (3, NODE)
      IF (NODE .NE. NOLD) GOTO 100

C  GET THE FIRST CORNER'S DISTANCE FROM PREVIOUS CORNER CORRECT

      IF (NCORN .GE. 1) LNODES (7, LCORN (NCORN ) ) = KKC + KOUNTC

  120 CONTINUE

      RETURN

      END
