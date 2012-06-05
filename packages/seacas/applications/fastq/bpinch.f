C $Id: bpinch.f,v 1.2 1991/03/21 15:44:21 gdsjaar Exp $
C $Log: bpinch.f,v $
C Revision 1.2  1991/03/21 15:44:21  gdsjaar
C Changed all 3.14159... to atan2(0.0, -1.0)
C
c Revision 1.1.1.1  1990/11/30  11:04:11  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:04:10  gdsjaar
c Initial revision
c 
C
CC* FILE: [.PAVING]BPINCH.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE BPINCH (MXND, MLN, LNODES, XN, YN, LXN, NXL, ANGLE,
     &   N0, N1, N2, NLOOP, TOLER1, TOLER2, BOK, ERR)
C***********************************************************************
C
C  SUBROUTINE BPINCH = CHECKS THAT A PINCH IS ALLOWABLE AND THAT IT
C                      DOESN'T FORCE A  DEGENERATE BOUNDARY ELEMENT
C
C***********************************************************************
C
      DIMENSION XN (MXND), YN (MXND)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)
      DIMENSION ANGLE (MXND), LNODES (MLN, MXND)
      DIMENSION L1LIST(20)
C
      LOGICAL BOK, CORNP, PANGLE, ERR
C
      TWOPI = 2.0 * ATAN2(0.0, -1.0)
C
C  SEE IF THE ANGLE IS ELIGIBLE FOR PINCHING
C  FIRST CHECK A NONBOUNDARY NODE
C
      IF (LXN (2, N1) .GT. 0) THEN
C
C  CHECK A FOUR (OR LESS) LINE NODE
C
         IF (LXN (4, N1) .GE. 0) THEN
            IF (ANGLE (N1) .LT. TOLER1) THEN
               PANGLE = .TRUE.
            ELSE
               PANGLE = .FALSE.
            ENDIF
C
C  CHECK A FIVE (OR MORE) LINE NODE
C
         ELSE
            IF (ANGLE (N1) .LT. TOLER2) THEN
               PANGLE = .TRUE.
            ELSE
               PANGLE = .FALSE.
            ENDIF
         ENDIF
C
C  CHECK A BOUNDARY NODE
C
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
C
C  ALL THREE ARE NOT ON THE BOUNDARY
C
         IF ( (LXN (2, N1) .GT. 0) .AND.
     &      (LXN (2, N0) .GT. 0) .AND.
     &      (LXN (2, N2) .GT. 0) ) THEN
            BOK = .TRUE.
C
C  N0 AND N2 ARE ON THE BOUNDARY
C
         ELSEIF ( (LXN (2, N0) .LT. 0) .AND.
     &      (LXN (2, N2) .LT. 0) ) THEN
            BOK = .FALSE.
C
C  N1 AND N0 ARE ON THE BOUNDARY - FIND THE ANGLE THAT
C  THE BOUNDARY AT N1 MAKES
C
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
C
C  NOW CHECK TO MAKE SURE THAT ANGLE IS NOT TOO LARGE
C
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
C
C  N1 AND N2 ARE ON THE BOUNDARY
C
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
C
C  NOW CHECK THE ANGLE SIZE
C
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
C
C  ONLY N0 IS ON THE BOUNDARY
C
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
C
C  ONLY N1 IS ON THE BOUNDARY
C
         ELSEIF (LXN (2, N1) .LT. 0) THEN
            BOK = .TRUE.
C
C  ONLY N2 IS ON THE BOUNDARY
C
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
C
         ENDIF
C
      ELSE
         BOK = .FALSE.
      ENDIF
C
  140 CONTINUE
      RETURN
C
      END
