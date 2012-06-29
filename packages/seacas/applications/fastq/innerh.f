C $Id: innerh.f,v 1.2 1998/07/14 18:19:13 gdsjaar Exp $
C $Log: innerh.f,v $
C Revision 1.2  1998/07/14 18:19:13  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.1.1.1  1990/11/30 11:09:56  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:09:54  gdsjaar
c Initial revision
c 
C
CC* FILE: [.QMESH]INNERH.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE INNERH (MXND, NXH, NUID, LXK, KXL, NXL, LXN, KKK, LLL,
     &   NNN, NNNOLD, NH, ISTART, IAVAIL, NAVAIL, NOROOM, ERR)
C***********************************************************************
C
C  SUBROUTINE INNERH  =  INSERT A ROW OF ELEMENTS AROUND A HOLE
C
C***********************************************************************
C
      DIMENSION NUID(MXND), LXK(4, MXND)
      DIMENSION KXL(2, 3*MXND), NXL(2, 3*MXND), LXN(4, MXND)
      DIMENSION NXH(MXND)
C
      LOGICAL ERR, NOROOM
C
      ERR = .FALSE.
      IF ((KKK + NH .GT. MXND) .OR. (LLL + 2*NH .GT. 3*MXND) .OR.
     &   (NNN + NH .GT. MXND)) THEN
         NOROOM = .TRUE.
         GO TO 200
      ELSE
         NOROOM = .FALSE.
      END IF
C
C  GENERATE LINES
C
      LLLOLD = LLL
      I1 = ISTART
      DO 100 I = 1, NH
         LLL = LLL + 1
         NXL(1, LLL) = NXH(I1)
         CALL ADDLXN(MXND, LXN, NUID, NAVAIL, IAVAIL, NXL(1, LLL),
     &      LLL, NNN, ERR, NOROOM)
         IF (NOROOM .OR. ERR) GO TO 200
         NXL(2, LLL) = NNNOLD + I
         I1 = I1 + 1
         IF (I1 .GT. NH) I1 = 1
  100 CONTINUE
C
      DO 110 I = 1, NH - 1
         LLL = LLL + 1
         NXL(1, LLL) = NNNOLD + I
         NXL(2, LLL) = NNNOLD + I + 1
  110 CONTINUE
      LLL = LLL + 1
      NXL(1, LLL) = NNN
      NXL(2, LLL) = NNNOLD + 1
C
C  MAKE SPACE IN LXN LIST FOR NEW NODES
C
      IF (IAVAIL .LT. NNNOLD) THEN
         JJ = IAVAIL
         MAXLNK = 0
  120    CONTINUE
         IF (LXN(4, JJ) .LE. NNNOLD) THEN
            JJ = LXN(4, JJ)
            MAXLNK = MAXLNK + 1
            IF (MAXLNK .GT. NNNOLD) THEN
               CALL MESAGE ('INNERH - LXN LINKED LIST ERROR')
               ERR = .TRUE.
               GO TO 200
            ELSE
               GO TO 120
            END IF
         END IF
         LXN(4, JJ) = NNNOLD + NH + 1
         NAVAIL = NAVAIL - NH
      ELSE IF (IAVAIL .GT. NNNOLD) THEN
         DO 130 J = 1, IAVAIL - 1
            IF (LXN(4, J) .LT. 0 .AND. ABS(LXN(4, J)) .GT. NNNOLD)
     &         LXN(4, J) = LXN(4, J) - NH
  130    CONTINUE
         JJ = NNNOLD + NH
         DO 150 J = NNNOLD + 1, IAVAIL - 1
            JJ = JJ + 1
            DO 140 K = 1, 4
               LXN(K, JJ) = LXN(K, J)
  140       CONTINUE
  150    CONTINUE
         IAVAIL = JJ + 1
         NAVAIL = NAVAIL - NH
      ELSE
         IAVAIL = NNN + 1
         NAVAIL = MXND - NNN
      END IF
C
C  MARK NODES ON HOLE BOUNDARY
C
      LXN(1, NNNOLD + 1) = LLLOLD + 1
      LXN(2, NNNOLD + 1) = -(LLLOLD + 2*NH)
      LXN(3, NNNOLD + 1) = LLLOLD + NH + 1
      LXN(4, NNNOLD + 1) = 0
      NNN = NNNOLD + 1
      DO 160 I = 2, NH
         NNN = NNN + 1
         LXN(1, NNN) = LLLOLD + I
         LXN(2, NNN) = -(LLLOLD + NH + I - 1)
         LXN(3, NNN) = LLLOLD + NH + I
         LXN(4, NNN) = 0
  160 CONTINUE
C
C  GENERATE ELEMENTS
C
      DO 180 I = LLLOLD + 1, LLL
         DO 170 J = 1, 2
            KXL(J, I) = 0
  170    CONTINUE
  180 CONTINUE
C
      I1 = ISTART
      I2 = I1 + 1
      IF (I2 .GT. NH) I2 = 1
      DO 190 I = 1, NH - 1
         CALL FNDLIN (MXND, LXN, NXH(I1), NXH(I2), LINE, ERR)
         IF (ERR) GO TO 200
         KKK = KKK + 1
         LXK(1, KKK) = LINE
         LXK(2, KKK) = LLLOLD + I
         LXK(3, KKK) = LLLOLD + I + 1
         LXK(4, KKK) = LLLOLD + I + NH
         I1 = I1 + 1
         IF (I1 .GT. NH) I1 = 1
         I2 = I1 + 1
         IF (I2 .GT. NH) I2 = 1
C
         IF (KXL(1, LINE) .EQ. 0) THEN
            KXL(1, LINE) = KKK
         ELSE IF (KXL(2, LINE) .EQ. 0) THEN
            KXL(2, LINE) = KKK
         ELSE
            CALL MESAGE ('KXL TABLE FULL')
         END IF
         KXL(2, LLLOLD + I) = KKK
         KXL(1, LLLOLD + I + 1) = KKK
         KXL(1, LLLOLD + I + NH) = KKK
  190 CONTINUE
      CALL FNDLIN (MXND, LXN, NXH(I1), NXH(I2), LINE, ERR)
      IF (ERR) GO TO 200
      KKK = KKK + 1
      LXK(1, KKK) = LINE
      LXK(2, KKK) = LLLOLD + NH
      LXK(3, KKK) = LLLOLD + 1
      LXK(4, KKK) = LLLOLD + 2*NH
C
      IF (KXL(1, LINE) .EQ. 0) THEN
         KXL(1, LINE) = KKK
      ELSE IF (KXL(2, LINE) .EQ. 0) THEN
         KXL(2, LINE) = KKK
      ELSE
         CALL MESAGE ('KXL TABLE FULL')
      END IF
      KXL(2, LLLOLD + NH) = KKK
      KXL(1, LLLOLD + 1) = KKK
      KXL(1, LLLOLD + 2*NH) = KKK
C
  200 CONTINUE
      RETURN
      END
