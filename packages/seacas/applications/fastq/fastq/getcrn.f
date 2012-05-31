C $Id: getcrn.f,v 1.1 1990/11/30 11:07:56 gdsjaar Exp $
C $Log: getcrn.f,v $
C Revision 1.1  1990/11/30 11:07:56  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]GETCRN.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE GETCRN (MXND, MXCORN, MLN, LNODES, NCORN, LCORN,
     &   ANGLE, XN, YN, LXN, NLOOP, N1, ONLYC, PPOSBL, GRAPH, ERR)
C***********************************************************************
C
C  SUBROUTINE GETCRN = SETS UP ALL THE POSSIBLE CORNER (OR NON-SIDE)
C                      LOCATIONS IN THE MESH
C
C***********************************************************************
C
      DIMENSION LNODES (MLN, MXND), LCORN (MXCORN), ANGLE (MXND)
      DIMENSION XN (MXND), YN (MXND), LXN (4, MXND), X(1), Y(1)
C
      LOGICAL ONLYC, GRAPH, PPOSBL, ERR
C
      ERR = .FALSE.
C
C  COUNT THE CURRENT "POSSIBLE NON-SIDES" STARTING AT THE I COUNTER
C
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
C
C  CHECK IF A PRIMITIVE IS EVEN POSSIBLE HERE
C
      IF ( LNODES (6, NODE) .GT. 4) PPOSBL = .FALSE.
C
C  CHECK FOR A POSSIBLE "NON - SIDE" NODE
C
      IF (ONLYC) CALL NDSTAT (NODE, LXN (1, NODE), ANGLE (NODE), ISTAT)
C
C  A NEW "POSSIBLE NON-SIDE" NODE HAS BEEN FOUND
C
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
C
C  TOO MANY "NON-SIDES" HAVE BEEN FOUND - JUST PICK ONE AND GO
C
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
C
C  INPUT THIS "NON-SIDE"
C
         NCORN = NCORN + 1
         LCORN (NCORN) = NODE
C
C  ADD UP THE NUMBER OF NODES FROM THE LAST "NON-SIDE"
C
         IF (NCORN .GT. 1) THEN
            LNODES (7, LCORN (NCORN - 1)) = KOUNTC + 1
         ELSE
            KKC = KOUNTC + 1
         ENDIF
         KOUNTC = 0
C
C  THIS IS A SIDE - JUST CONTINUE
C
      ELSE
         KOUNTC = KOUNTC + 1
      ENDIF
C
C  CHECK FOR COMPLETION OF THE LOOP
C
      NODE = LNODES (3, NODE)
      IF (NODE .NE. NOLD) GOTO 100
C
C  GET THE FIRST CORNER'S DISTANCE FROM PREVIOUS CORNER CORRECT
C
      IF (NCORN .GE. 1) LNODES (7, LCORN (NCORN ) ) = KKC + KOUNTC
C
  120 CONTINUE
C
      RETURN
C
      END
