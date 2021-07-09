C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE RESTRY (MXND, K, K2, LXK, NXL, KXL, LXN, XN, YN, NUID,
     &   NAVAIL, IAVAIL, NNN, DONE, ERR, NOROOM)
C***********************************************************************

C  SUBROUTINE RESTRY = TRY TO RESTRUCTURE K AND ONE OF ITS NEIGHBORS

C***********************************************************************

C  NOTE:
C     THE ELEMENT OPPOSITE THE LONGEST  -  OR IN SOME CASES THE SECOND
C     LONGEST  -  SIDE OF ELEMENT K WILL BE FOUND.
C     A RESTRUCTURE WILL ONLY BE DONE THEN IF
C          - THIS SAME SIDE IS THE LONGEST SIDE  (OR NEARLY SO)
C          OF THE SECOND ELEMENT ALSO.
C          - NO ANGLES LARGER THAN 180 DEGREES ARE CREATED
C          - THE MAX AND AVERAGE Q - NUMBERS BOTH DECREASE
C          - AN AREA BREAKDOWN CONSISTENT WITH THE IMPROVEMENT
C          IN THE AVERAGE Q - NUMBER OCCURS

C***********************************************************************

      DIMENSION LXK (4, MXND), NXL (2, 3*MXND)
      DIMENSION KXL (2, 3*MXND), XN (MXND)
      DIMENSION YN (MXND), NUID (MXND), LXN (4, MXND)
      DIMENSION NSA (4), NSB (4), NSC (4), NSD (4)
      DIMENSION NS1 (4), AL1 (4), ANG1 (4)
      DIMENSION NS2 (4), AL2 (4), ANG2 (4)

      LOGICAL ERR, CCW, CAREA, IOKB, IOKF, IBGOOD
      LOGICAL IFGOOD, DONE, LSIDE, NOROOM

      DONE = .FALSE.
      ERR = .FALSE.
      PI = ATAN2(0.0, -1.0)
      PITOL = PI* (195. / 180.)

C  GET DATA FOR ELEMENT K

      CCW = .TRUE.
      CAREA = .FALSE.
      CALL GNXKA (MXND, XN, YN, K, NS1, AREA1, LXK, NXL, CCW)
      CALL QAAVAL (MXND, NS1, ANG1, QRAT1, DUMMY, XN, YN, CAREA)
      CALL CONDNO (MXND, NS1, QRAT1, SRAT1, COND1, AL1, XN, YN, LSIDE)

C  FIND LONGEST AND SECOND LONGEST SIDES,  EXCLUDING BOUNDARY LINES

      S1MAX =  - 1.
      L1MAX = 0
      A1MAX = 0.
      S2MAX =  - 1.
      DO 130 I = 1, 4
         IF (AL1 (I) .GT. S2MAX) THEN
            N1 = NS1 (I)
            J = I + 1
            IF (I .EQ. 4)J = 1
            N2 = NS1 (J)
            DO 100 IL = 1, 4
               L = LXK (IL, K)
               NODE = NXL (1, L)
               IF ( (NODE  .EQ.  N1) .OR. (NODE  .EQ.  N2)) THEN
                  NODE = NXL (2, L)
                  IF ( (NODE  .EQ.  N1) .OR. (NODE  .EQ.  N2)) THEN
                     IF (KXL (2, L) .GT. 0) THEN
                        GOTO 110
                     ELSE
                        GOTO 120
                     ENDIF
                  ENDIF
               ENDIF
  100       CONTINUE
            WRITE (*, 10000)N1, N2, K
            ERR = .TRUE.
            RETURN

C  N1 TO N2 IS AN INTERIOR LINE

  110       CONTINUE

C  LONGEST INTERIOR LINE SO FAR

            IF (AL1 (I) .GT. S1MAX) THEN
               S2MAX = S1MAX
               L2MAX = L1MAX
               A2MAX = A1MAX
               S1MAX = AL1 (I)
               L1MAX = L
               A1MAX = ANG1 (I) + ANG1 (J)

C  SECOND LONGEST LINE SO FAR

            ELSE
               S2MAX = AL1 (I)
               L2MAX = L
               A2MAX = ANG1 (I) + ANG1 (J)
            ENDIF
         ENDIF
  120    CONTINUE
  130 CONTINUE

C***********************************************************************
C  NOTE:
C     IF LONGEST SIDE IS SUBSTANTIALLY LONGER THAN SECOND
C     LONGEST,  PAIR WITH ELEMENT OPPOSITE LONGEST SIDE.
C     IF LONGEST SIDE IS NOT SUBSTANTIALLY LONGER THAN SECOND
C     LONGEST SIDE,  PAIR WITH THE ELEMENT OPPOSITE EITHER THE
C     LONGEST OR SECOND LONGEST SIDE DEPENDING ON WHICH HAS THE
C     SMALLER SUM OF ADJACENT ANGLES.   (LOOK AT A TRAPEZOID TO
C     SEE WHY THIS IS REASONABLE.)

C***********************************************************************

      IF (L1MAX .LE. 0) RETURN
      SLEN  =  S1MAX
      LINT  =  L1MAX
      IF ( (L2MAX .GT. 0) .AND. (S1MAX .LE. 1.15 * S2MAX) .AND.
     &   (A1MAX .GT. A2MAX)) THEN
         SLEN  =  S2MAX
         LINT  =  L2MAX
      ENDIF
      K2  =  KXL (1, LINT) + KXL (2, LINT) - K

C  DOUBLE CHECK

      IF (K2 .LE. 0) RETURN

C  GET DATA FOR ELEMENT K2

      CALL GNXKA (MXND, XN, YN, K2, NS2, AREA2, LXK, NXL, CCW)
      CALL QAAVAL (MXND, NS2, ANG2, QRAT2, DUMMY, XN, YN, CAREA)
      CALL CONDNO (MXND, NS2, QRAT2, SRAT2, COND2, AL2, XN, YN, LSIDE)

C  FIND LONGEST SIDE IN SECOND ELEMENT

      SMAXB  =  AMAX1 (AL2 (1), AL2 (2), AL2 (3), AL2 (4))

C  IF THE INTERFACE SIDE IS SIGNIFICANTLY SHORTER THAN THE
C  LONGEST SIDE OF THE SECOND ELEMENT,  SKIP THE RESTRUCTURE.

      IF (SLEN .LT. 0.50 * SMAXB) RETURN

C  CIRCULARLY SHIFT THE TWO NODE LISTS TO CREATE CANONICAL ORDER.
C  IN CANONICAL ORDER THE FIRST NODE IS THE NODE IN BOTH ELEMENTS
C  WHOSE COUNTER - CLOCKWISE SUCCESSOR IN THE FIRST ELEMENT IS NOT
C  ALSO IN THE SECOND ELEMENT.
C  NOTE : ORDER OF SIDE LENGTH AND ANGLE DATA IS NO GOOD AFTER THIS

      N1  =  NXL (1, LINT)
      CALL NXKORD (NS1, N1)
      DO 140 I  =  1, 4
         IF (NS2 (I)  .EQ.  NS1 (2)) THEN
            N1  =  NXL (2, LINT)
            CALL NXKORD (NS1, N1)
            CALL NXKORD (NS2, N1)
            GOTO 150
         ENDIF
  140 CONTINUE
      CALL NXKORD (NS2, N1)
  150 CONTINUE

C  SEE IF THEY MATCH AS THEY SHOULD  (BUTTERFLY ELEMENTS MAY CAUSE
C  PROBLEMS WITH THE CCW ROUTINES

      IF (NS1 (4) .NE. NS2 (2)) THEN
         NSHOLD  =  NS2 (2)
         NS2 (2)  =  NS2 (4)
         NS2 (4)  =  NSHOLD
         IF (NS1 (4) .NE. NS2 (2)) THEN
            NSHOLD  =  NS1 (2)
            NS1 (2)  =  NS1 (4)
            NS1 (4)  =  NSHOLD
            IF (NS1 (4) .NE. NS2 (2)) THEN

C  ERROR MATCHING ELEMENTS ALONG A COMMON SIDE

               WRITE ( * , 10010)K, K2, NS1 (1), NS1 (4)
               ERR  =  .TRUE.
               RETURN
            ENDIF
         ENDIF
      ENDIF

C  COMPUTE ALL RELEVANT DATA FOR ALL THREE STRUCTURES

C  ORIGINAL STRUCTURE

      QMAX  =  AMAX1 (QRAT1, QRAT2)
      TOLQX  =  .95 * QMAX + .05
      QAVG  =  0.5 *  (QRAT1 + QRAT2)
      AMIN  =  AMIN1 (AREA1, AREA2)
      IF (AMIN .GT. 0.) THEN
         ARAT  =  AMAX1 (AREA1, AREA2) / AMIN
      ELSE
         ARAT  =  1.0E10
      ENDIF

C   * BACKWARDS *  STRUCTURE

      NSA (1)  =  NS2 (4)
      NSA (2)  =  NS1 (1)
      NSA (3)  =  NS1 (2)
      NSA (4)  =  NS1 (3)
      NSB (1)  =  NS2 (4)
      NSB (2)  =  NS1 (3)
      NSB (3)  =  NS2 (2)
      NSB (4)  =  NS2 (3)
      CAREA  =  .TRUE.
      CALL QAAVAL (MXND, NSA, ANG1, QRAT1B, AREA1B, XN, YN, CAREA)
      CALL QAAVAL (MXND, NSB, ANG2, QRAT2B, AREA2B, XN, YN, CAREA)
      IF ( (AMAX1 (ANG1 (1), ANG1 (2), ANG1 (3), ANG1 (4)) .GT. PITOL)
     &   .OR.
     &   (AMAX1 (ANG2 (1), ANG2 (2), ANG2 (3), ANG2 (4)) .GT. PITOL))
     &   THEN
         IOKB  =  .FALSE.
      ELSE
         IOKB  =  .TRUE.
      ENDIF
      QMAXB  =  AMAX1 (QRAT1B, QRAT2B)
      QAVGB  =  0.5 *  (QRAT1B + QRAT2B)
      AMIN  =  AMIN1 (AREA1B, AREA2B)
      IF (AMIN .GT. 0.) THEN
         ARATB  =  AMAX1 (AREA1B, AREA2B) / AMIN
      ELSE
         ARATB  =  1.0E10
      ENDIF

C   * FORWARDS *  STRUCTURE

      NSC (1)  =  NS1 (2)
      NSC (2)  =  NS1 (3)
      NSC (3)  =  NS1 (4)
      NSC (4)  =  NS2 (3)
      NSD (1)  =  NS1 (2)
      NSD (2)  =  NS2 (3)
      NSD (3)  =  NS2 (4)
      NSD (4)  =  NS2 (1)
      CALL QAAVAL (MXND, NSC, ANG1, QRAT1F, AREA1F, XN, YN, CAREA)
      CALL QAAVAL (MXND, NSD, ANG2, QRAT2F, AREA2F, XN, YN, CAREA)
      IF ( (AMAX1 (ANG1 (1), ANG1 (2), ANG1 (3), ANG1 (4)) .GT. PITOL)
     &   .OR.
     &   (AMAX1 (ANG2 (1), ANG2 (2), ANG2 (3), ANG2 (4)) .GT. PITOL))
     &   THEN
         IOKF  =  .FALSE.
      ELSE
         IOKF  =  .TRUE.
      ENDIF
      QMAXF  =  AMAX1 (QRAT1F, QRAT2F)
      QAVGF  =   0.5 *  (QRAT1F + QRAT2F)
      AMIN  =  AMIN1 (AREA1F, AREA2F)
      IF (AMIN .GT. 0.) THEN
         ARATF  =  AMAX1 (AREA1F, AREA2F) / AMIN
      ELSE
         ARATF  =  1.0E10
      ENDIF

C  SEE IF BACKWARD IS BETTER THAN ORIGINAL

      IF ( (IOKB) .AND. (QMAXB .LE. TOLQX) .AND. (QAVGB .LE. QAVG) .AND.
     &   (ARATB * QAVGB .LE. ARAT * QAVG)) THEN
         IBGOOD  =  .TRUE.
      ELSE
         IBGOOD  =  .FALSE.
      ENDIF

C  SEE IF FORWARD IS BETTER THAN ORIGINAL

      IF ( (IOKF) .AND. (QMAXF .LE. TOLQX) .AND. (QAVGF .LE. QAVG) .AND.
     &   (ARATF * QAVGF .LE. ARAT * QAVG)) THEN
         IFGOOD  =  .TRUE.
      ELSE
         IFGOOD  =  .FALSE.
      ENDIF

C  CHOOSE BEST ALTERNATIVE
C  IF BOTH FORWARD AND BACKWARD IS BETTER THAN ORIGINAL,  THEN
C  COMPUTE PAIR - VALUES TO CHOOSE BETWEEN FORWARD AND BACKWARD.
C  VALUE  =   (AVERAGE CONDITION NUMBER)  *  SQRT (AREA RATIO)

      IF ( (IFGOOD) .AND. (IBGOOD)) THEN
         LSIDE  =  .FALSE.
         CALL CONDNO (MXND, NSA, QRAT1B, SRAT1B, COND1B, AL1, XN, YN,
     &      LSIDE)
         CALL CONDNO (MXND, NSB, QRAT2B, SRAT2B, COND2B, AL1, XN, YN,
     &      LSIDE)
         VALUEB  =  ARATB *  (COND1B + COND2B) **2
         CALL CONDNO (MXND, NSC, QRAT1F, SRAT1F, COND1F, AL1, XN, YN,
     &      LSIDE)
         CALL CONDNO (MXND, NSD, QRAT2F, SRAT2F, COND2F, AL1, XN, YN,
     &      LSIDE)
         VALUEF  =  ARATF *  (COND1F + COND2F) **2
         IF (VALUEB .GT. VALUEF)IBGOOD  =  .FALSE.
      ENDIF

C  BACKWARD STRUCTURE IS BEST.  IMPLEMENT IT.

C  FIRST FIX LXK AND KXL ARRAYS

      IF (IBGOOD) THEN
         CALL FNDLNK (MXND, LXK, NXL, K, NS1 (3), NS1 (4), L1EE, ERR)
         IF (ERR) RETURN
         CALL FNDLNK (MXND, LXK, NXL, K2, NS2 (1), NS2 (4), L2EE, ERR)
         IF (ERR) RETURN
         CALL LSWAP (MXND, LXK, KXL, K, L1EE, K2, L2EE, ERR)
         IF (ERR) RETURN

C  FIX NXL ARRAY  (MOVE THE DIAGONAL)

         NXL (1, LINT)  =  NS2 (4)
         NXL (2, LINT)  =  NS1 (3)

C  FIX LXN ARRAY

         CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, NS1 (1), LINT,
     &      NNN, ERR, NOROOM)
         IF (ERR) RETURN
         CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, NS1 (4), LINT,
     &      NNN, ERR, NOROOM)
         IF (ERR) RETURN
         CALL ADDLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, NS1 (3), LINT,
     &      NNN, ERR, NOROOM)
         IF (ERR) RETURN
         CALL ADDLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, NS2 (4), LINT,
     &      NNN, ERR, NOROOM)
         IF (ERR) RETURN
         DONE  =  .TRUE.

C  FORWARD STRUCTURE IS BEST.  IMPLEMENT IT.

C  FIX LXK AND KXL ARRAYS

      ELSEIF (IFGOOD) THEN
         CALL FNDLNK (MXND, LXK, NXL, K, NS1 (1), NS1 (2), L1EE, ERR)
         IF (ERR) RETURN
         CALL FNDLNK (MXND, LXK, NXL, K2, NS2 (2), NS2 (3), L2EE, ERR)
         IF (ERR) RETURN
         CALL LSWAP (MXND, LXK, KXL, K, L1EE, K2, L2EE, ERR)
         IF (ERR) RETURN

C  FIX NXL ARRAY  (MOVE THE DIAGONAL)

         NXL (1, LINT)  =  NS1 (2)
         NXL (2, LINT)  =  NS2 (3)

C  FIX LXN ARRAY

         CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, NS1 (1), LINT,
     &      NNN, ERR, NOROOM)
         IF (ERR) RETURN
         CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, NS1 (4), LINT,
     &      NNN, ERR, NOROOM)
         IF (ERR) RETURN
         CALL ADDLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, NS1 (2), LINT,
     &      NNN, ERR, NOROOM)
         IF (ERR) RETURN
         CALL ADDLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, NS2 (3), LINT,
     &      NNN, ERR, NOROOM)
         IF (ERR) RETURN
         DONE  =  .TRUE.
      ENDIF

      RETURN

10000 FORMAT (' IN RESTRY, NODES', 2I5, /,
     &   ' DO NOT DEFINE A LINE IN ELEMENT', I5)
10010 FORMAT (' IN RESTRY,  ELEMENTS', 2I5, /,
     &   ' DO NOT CONTAIN A COMMON SIDE USING NODES', 2I5)

      END
