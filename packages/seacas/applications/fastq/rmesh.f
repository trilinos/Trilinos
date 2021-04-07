C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE RMESH (NPER, MXND, X, Y, NID, XN, YN, NUID, LXK, KXL,
     &   NXL, LXN, M1, M2, KKK, KKKOLD, NNN, NNNOLD, LLL, LLLOLD,
     &   IAVAIL, NAVAIL, ERR)
C***********************************************************************

C  SUBROUTINE RMESH = GENERATES AN INITIAL LOGICAL RECTANGULAR MESH
C                     WHOSE PERIMETER IS  (X (I), Y (I), I=1, N).

C***********************************************************************

C  VARIABLES USED:
C     X     = X VALUES AROUND THE PERIMETER
C     Y     = Y VALUES AROUND THE PERIMETER
C     NID   = PERIMETER NODE UNIQUE ID'S
C     N     = NUMBER OF PERIMETER NODES
C     M1    = THE NUMBER OF INTERVALS ON THE FIRST SIDE OF THE RECTANGLE
C     IMAP  = CONTROLS THE DEFINITION OF THE VALUES OF THE COORDINATE
C             OF THE INTERIOR NODES.
C           = 1 MEANS SET ALL COORDINATES TO 0.
C           = 2 MEANS SET ALL COORDINATES TO THE CENTROID OF PERIMETER
C           = 3 MEANS APPLY THE UNIT SQUARE TRANSFORMATION FROM
C             W.A.COOK'S PAPER  (THIOKOL REPORT AFRDL - TR - 71 - 51)
C     ERR   = .TRUE. IF ERRORS WERE ENCOUNTERED
C     XN    = GLOBAL X VALUES OF NODES
C     YN    = GLOBAL Y VALUES OF NODES
C     NUID  = GLOBAL NODE UNIQUE IDENTIFIERS
C     LXK   = LINES PER ELEMENT
C     KXL   = ELEMENTS PER LINE
C     NXL   = NODES PER LINE
C     LXN   = LINES PER NODE
C  NOTE:
C     FOR *XN TABLES A NEGATIVE FLAG IN THE FOURTH COLUMN MEANS
C     GO TO THAT ROW FOR A CONTINUATION OF THE LIST.  IN THAT ROW
C     THE FIRST ELEMENT WILL BE NEGATED TO INDICATE THAT THIS IS
C     A CONTINUATION ROW.   (RMESH ITSELF GENERATES NO SUCH NEGATIVES.)
C     A NEGATIVE FLAG IN THE SECOND COLUMN OF THE LXN ARRAY MEANS
C     THAT THIS NODE IS A BOUNDARY NODE.

C***********************************************************************

      DIMENSION X (NPER), Y (NPER), NID (NPER)
      DIMENSION XN (MXND), YN (MXND), NUID (MXND)
      DIMENSION LXK (4, MXND), KXL (2, 3 * MXND)
      DIMENSION NXL (2, 3 * MXND), LXN (4, MXND)

      LOGICAL ERR, NOROOM

C  NOTE: NOROOM SHOULD NEVER BE TRUE WITH CROSS BEING CALLED IN RMESH.
C        THUS,  IT IS NEVER PASSED BACK TO QMESH.

      KKK = M1 * M2 + KKKOLD
      LLL =  (M1 *  (M2 + 1)) +  (M2 *  (M1 + 1)) + LLLOLD
      NNN =  (M1 + 1) *  (M2 + 1) + NNNOLD
      ERR = .TRUE.

C  CHECK INPUT

      IF (2 *  (NPER / 2) .NE. NPER) THEN
         CALL MESAGE ('IN RMESH,  NO. OF PERIMETER NODES IS ODD')
         RETURN
      ELSEIF ( (M1 .LT. 1) .OR. (M2 .LT. 1)) THEN
         WRITE ( * , 10000)NPER, M1
         RETURN
      ENDIF

C  COMPUTE CONSTANTS

      NLC = 2 * M1 + 1
      M1P1 = M1 + 1
      M2P1 = M2 + 1

C  PRODUCE LXK ARRAY

C  LINES FOR FIRST ELEMENT

      LXK (1, KKKOLD + 1) = 1 + LLLOLD
      LXK (2, KKKOLD + 1) = M1 + 1 + LLLOLD
      LXK (3, KKKOLD + 1) = M1 + 2 + LLLOLD
      LXK (4, KKKOLD + 1) = NLC + 1 + LLLOLD

C  FIRST ROW  (SHIFT FIRST ELEMENT TO SECOND,  ETC.)

      IF (M1 .GT. 1) THEN
         DO 110 K = 2, M1
            DO 100 I = 1, 4
               LXK (I, K + KKKOLD) = LXK (I, K + KKKOLD - 1) + 1
  100       CONTINUE
  110    CONTINUE
      ENDIF

C  SUCCEEDING ROWS  (SHIFT FIRST COLUMN TO SECOND, ETC.)

      IF (M2 .GT. 1) THEN
         K = M1 + KKKOLD
         DO 140 K2 = 2, M2
            DO 130 K1 = 1, M1
               K = K + 1
               KL = K - M1
               DO 120 I = 1, 4
                  LXK (I, K) = LXK (I, KL) + NLC
  120          CONTINUE
  130       CONTINUE
  140    CONTINUE
      ENDIF

C  PREPARE KXL TABLE BY USING SUBROUTINE CROSS ON THE LXK TABLE

      CALL CCROSS (4, KKK, 2, LLL, LXK, KXL, KKKOLD + 1, LLLOLD + 1,
     &   NOROOM, ERR)
      IF (ERR) RETURN
      ERR = .TRUE.

C  MAKE NXL TABLE

C  FIRST DO HORIZONTAL LINES

      DO 160 J = 1, M2P1
         NODE = 1 +  (J - 1) *  (M1 + 1) + NNNOLD
         L = 1 +  ( (J - 1) * NLC) + LLLOLD
         DO 150 I = 1, M1
            NXL (1, L) = NODE
            NXL (2, L) = NODE + 1
            NODE = NODE + 1
            L = L + 1
  150    CONTINUE
  160 CONTINUE

C  NEXT DO VERTICAL LINES

      DO 180 J = 1, M1P1
         NODE = J + NNNOLD
         L = J + M1 + LLLOLD
         DO 170 I = 1, M2
            NXL (1, L) = NODE
            NXL (2, L) = NODE + M1P1
            NODE = NODE + M1P1
            L = L + NLC
  170    CONTINUE
  180 CONTINUE

C  PREPARE LXN TABLE FROM NXL TABLE

      CALL CCROSS (2, LLL, 4, NNN, NXL, LXN, LLLOLD + 1, NNNOLD + 1,
     &   NOROOM, ERR)
      IF (ERR) RETURN
      ERR = .TRUE.

C  LINK - UP AVAILABLE LXN SPACE

      IAVAIL = NNN + 1
      NAVAIL = MXND - NNN
      DO 190 I = IAVAIL, MXND
         LXN (1, I) = 0
         LXN (2, I) = 0
         LXN (3, I) = 0
         LXN (4, I) = I + 1
  190 CONTINUE

C  LOGICAL CONNECTION TABLES ARE COMPLETE
C  FILL IN THE CO - ORDINATES OF THE INTERIOR POINTS
C  USE THE UNIT SQUARE TRANSFORMATION  (OF COOK / THIOKOL)

      IF ( (M1 .GT. 1) .AND. (M2 .GT. 1)) THEN

C  GET NODE NUMBERS FOR CORNERS

         I1Z = M1 + 1
         IZ1 = NPER -  (M2 - 1)
         I11 = M1 + M2 + 1

C  COLUMN LOOP

         DO 210 J = 2, M2
            KL = NPER + 2 - J
            KR = M1 + J
            ETA = DBLE(J - 1) / DBLE(M2)
            OMETA = 1.0 - ETA

C  ROW LOOP

            DO 200 I = 2, M1
               KB = I
               KT = IZ1 + 1 - I
               EPS = DBLE(I - 1) / DBLE(M1)
               OMEPS = 1.0 - EPS
               IM =  (J - 1) * M1P1 + I + NNNOLD
               XN (IM) =  (OMETA * X (KB)) + (ETA * X (KT)) +
     &            (OMEPS * X (KL)) + (EPS * X (KR)) -
     &            ((X (1) * OMETA * OMEPS) + (X (I1Z) * OMETA * EPS) +
     &            (X (IZ1) * ETA * OMEPS) + (X (I11) * ETA * EPS))
               YN (IM) =  (OMETA * Y (KB)) + (ETA * Y (KT)) +
     &            (OMEPS * Y (KL)) + (EPS * Y (KR)) -
     &            ((Y (1) * OMETA * OMEPS) + (Y (I1Z) * OMETA * EPS) +
     &            (Y (IZ1) * ETA * OMEPS) + (Y (I11) * ETA * EPS))
  200       CONTINUE
  210    CONTINUE
      ENDIF

C  DEFINE THE COORDINATES OF THE PERIMETER NODES.
C  ALSO FLAG SECOND ELEMENTS OF LXN ARRAY TO INDICATE
C  WHICH NODES ARE BOUNDARY NODES.
C  DEFINE UNIQUE NODE ID NUMBERS ALSO.

      DO 220 I = NNNOLD + 1, NNN
         NUID (I) = 0
  220 CONTINUE

C  BOTTOM

      IP = 0
      DO 230 I = 1, M1P1
         IM = I + NNNOLD
         IP = IP + 1
         LXN (2, IM) =  - LXN (2, IM)
         NUID (IM) = NID (IP)
         XN (IM) = X (IP)
         YN (IM) = Y (IP)
  230 CONTINUE

C  RIGHT

      IP = M1P1
      DO 240 I = 2, M2P1
         IP = IP + 1
         IM = IM + M1P1
         LXN (2, IM) =  - LXN (2, IM)
         NUID (IM) = NID (IP)
         XN (IM) = X (IP)
         YN (IM) = Y (IP)
  240 CONTINUE

C  TOP

      DO 250 I = 2, M1P1
         IP = IP + 1
         IM = IM - 1
         LXN (2, IM) =  - LXN (2, IM)
         NUID (IM) = NID (IP)
         XN (IM) = X (IP)
         YN (IM) = Y (IP)
  250 CONTINUE

C  LEFT

      DO 260 I = 2, M2
         IP = IP + 1
         IM = IM - M1P1
         LXN (2, IM) =  - LXN (2, IM)
         NUID (IM) = NID (IP)
         XN (IM) = X (IP)
         YN (IM) = Y (IP)
  260 CONTINUE

C  EXIT

      ERR = .FALSE.

      RETURN

10000 FORMAT (' IN RMESH, N = ', I5, ' AND M1 = ', I5,
     &   ' ARE INCOMPATIBLE')

      END
