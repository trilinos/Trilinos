C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE BPINCH (MXND, MLN, LNODES, XN, YN, LXN, NXL, ANGLE,
     &   N0, N1, N2, NLOOP, TOLER1, TOLER2, BOK, ERR)
C***********************************************************************

C  SUBROUTINE BPINCH = CHECKS THAT A PINCH IS ALLOWABLE AND THAT IT
C                      DOESN'T FORCE A  DEGENERATE BOUNDARY ELEMENT

C***********************************************************************

      DIMENSION XN (MXND), YN (MXND)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)
      DIMENSION ANGLE (MXND), LNODES (MLN, MXND)
      DIMENSION L1LIST(20)

      LOGICAL BOK, CORNP, PANGLE, ERR

      TWOPI = 2.0 * ATAN2(0.0, -1.0)

C  SEE IF THE ANGLE IS ELIGIBLE FOR PINCHING
C  FIRST CHECK A NONBOUNDARY NODE

      IF (LXN (2, N1) .GT. 0) THEN

C  CHECK A FOUR (OR LESS) LINE NODE

         IF (LXN (4, N1) .GE. 0) THEN
            IF (ANGLE (N1) .LT. TOLER1) THEN
               PANGLE = .TRUE.
            ELSE
               PANGLE = .FALSE.
            ENDIF

C  CHECK A FIVE (OR MORE) LINE NODE

         ELSE
            IF (ANGLE (N1) .LT. TOLER2) THEN
               PANGLE = .TRUE.
            ELSE
               PANGLE = .FALSE.
            ENDIF
         ENDIF

C  CHECK A BOUNDARY NODE

      ELSE
         IF ( (ANGLE (N1) .LT. TOLER1) .AND.
     &      (LXN (2, N0) * LXN (2, N2) .LT. 0) ) THEN
            PANGLE = .TRUE.
         ELSEIF ( (ANGLE (N1) .LT. TOLER1) .AND.
     &      (LXN (2, N0) .GT. 0) .AND.
     &      (LXN (2, N2) .GT. 0) ) THEN
            PANGLE = .TRUE.
         ELSE
            PANGLE = .FALSE.
         ENDIF
      ENDIF
      IF (PANGLE) THEN

C  ALL THREE ARE NOT ON THE BOUNDARY

         IF ( (LXN (2, N1) .GT. 0) .AND.
     &      (LXN (2, N0) .GT. 0) .AND.
     &      (LXN (2, N2) .GT. 0) ) THEN
            BOK = .TRUE.

C  N0 AND N2 ARE ON THE BOUNDARY

         ELSEIF ( (LXN (2, N0) .LT. 0) .AND.
     &      (LXN (2, N2) .LT. 0) ) THEN
            BOK = .FALSE.

C  N1 AND N0 ARE ON THE BOUNDARY - FIND THE ANGLE THAT
C  THE BOUNDARY AT N1 MAKES

         ELSEIF ( (LXN (2, N0) .LT. 0) .AND.
     &      (LXN (2, N1) .LT. 0) ) THEN
            CALL GETLXN (MXND, LXN, N1, L1LIST, NL, ERR)
            IF (ERR) THEN
               CALL MESAGE ('** PROBLEMS IN SEW2 FINDING LXN FOR J1 **')
               GOTO 140
            ENDIF
            DO 100 I = 1, NL
               LL = L1LIST (I)
               IF ( (LL .NE. LNODES (5, N0)) .AND.
     &            (LL .NE. LNODES (5, N1)) ) THEN
                  IP1 = NXL (1, LL)
                  IP2 = NXL (2, LL)
                  IF ((IP1 .EQ. N1) .AND. (LXN (2, IP2) .LT. 0)) THEN
                     NP = IP2
                     GOTO 110
                  ELSEIF ((IP2 .EQ. N1) .AND.
     &               (LXN (2, IP1) .LT. 0)) THEN
                     NP = IP1
                     GOTO 110
                  ENDIF
               ENDIF
  100       CONTINUE
            CALL MESAGE ('** PROBLEMS IN BPINCH FINDING N1 BOUNDARY'//
     &         ' ANGLE NODE **')
            GOTO 140
  110       CONTINUE
            ANG1 = ATAN2 (YN (N0) - YN (N1), XN (N0) - XN (N1))
            IF (ANG1 .LT. 0.) ANG1 = ANG1 + TWOPI
            ANG2 = ATAN2 (YN (NP) - YN (N1), XN (NP) - XN (N1))
            IF (ANG2 .LT. 0.) ANG2 = ANG2 + TWOPI
            ANG = ANG1 - ANG2
            IF (ANG .LT. 0.) ANG = ANG + TWOPI

C  NOW CHECK TO MAKE SURE THAT ANGLE IS NOT TOO LARGE

            IF (ANG .LT. 2.3561945) THEN
               IF (LXN (3, N1) .EQ. 0) THEN
                  BOK = .FALSE.
               ELSE
                  BOK = .TRUE.
               ENDIF
            ELSE
               IF (LXN (4, N1) .EQ. 0) THEN
                  BOK = .FALSE.
               ELSE
                  BOK = .TRUE.
               ENDIF
            ENDIF

C  N1 AND N2 ARE ON THE BOUNDARY

         ELSEIF ( (LXN (2, N1) .LT. 0) .AND.
     &      (LXN (2, N2) .LT. 0) ) THEN
            CALL GETLXN (MXND, LXN, N1, L1LIST, NL, ERR)
            IF (ERR) THEN
               CALL MESAGE ('** PROBLEMS IN SEW2 FINDING LXN FOR J1 **')
               GOTO 140
            ENDIF
            DO 120 I = 1, NL
               LL = L1LIST (I)
               IF ( (LL .NE. LNODES (5, N0)) .AND.
     &            (LL .NE. LNODES (5, N1)) ) THEN
                  IP1 = NXL (1, LL)
                  IP2 = NXL (2, LL)
                  IF ((IP1 .EQ. N1) .AND. (LXN (2, IP2) .LT. 0)) THEN
                     NP = IP2
                     GOTO 130
                  ELSEIF ((IP2 .EQ. N1) .AND.
     &               (LXN (2, IP1) .LT. 0)) THEN
                     NP = IP1
                     GOTO 130
                  ENDIF
               ENDIF
  120       CONTINUE
            CALL MESAGE ('** PROBLEMS IN BPINCH FINDING N1 BOUNDARY'//
     &         ' ANGLE NODE **')
            GOTO 140
  130       CONTINUE
            ANG1 = ATAN2 (YN (N2) - YN (N1), XN (N2) - XN (N1))
            IF (ANG1 .LT. 0.) ANG1 = ANG1 + TWOPI
            ANG2 = ATAN2 (YN (NP) - YN (N1), XN (NP) - XN (N1))
            IF (ANG2 .LT. 0.) ANG2 = ANG2 + TWOPI
            ANG = ANG2 - ANG1
            IF (ANG .LT. 0.) ANG = ANG + TWOPI

C  NOW CHECK THE ANGLE SIZE

            IF (ANG .LT. 2.3561945) THEN
               IF (LXN (3, N1) .EQ. 0) THEN
                  BOK = .FALSE.
               ELSE
                  BOK = .TRUE.
               ENDIF
            ELSE
               IF (LXN (4, N1) .EQ. 0) THEN
                  BOK = .FALSE.
               ELSE
                  BOK = .TRUE.
               ENDIF
            ENDIF

C  ONLY N0 IS ON THE BOUNDARY

         ELSEIF (LXN (2, N0) .LT. 0) THEN
            N0A = LNODES (2, N0)
            N0B = LNODES (2, N0A)
            IF ( (NLOOP .EQ. 6) .AND.
     &         (LXN (2, N0A) .LT. 0) .AND.
     &         (LXN (2, N0B) .LT. 0) .AND.
     &         (.NOT. CORNP (ANGLE (N0A)) ) ) THEN
               BOK = .FALSE.
            ELSE
               BOK = .TRUE.
            ENDIF

C  ONLY N1 IS ON THE BOUNDARY

         ELSEIF (LXN (2, N1) .LT. 0) THEN
            BOK = .TRUE.

C  ONLY N2 IS ON THE BOUNDARY

         ELSEIF (LXN (2, N2) .LT. 0) THEN
            N2A = LNODES (3, N2)
            N2B = LNODES (3, N2A)
            IF ( (NLOOP .EQ. 6) .AND.
     &         (LXN (3, N2A) .LT. 0) .AND.
     &         (LXN (3, N2B) .LT. 0) .AND.
     &         (.NOT. CORNP (ANGLE (N2A)) ) ) THEN
               BOK = .FALSE.
            ELSE
               BOK = .TRUE.
            ENDIF

         ENDIF

      ELSE
         BOK = .FALSE.
      ENDIF

  140 CONTINUE
      RETURN

      END
