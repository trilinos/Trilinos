C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE FXNUID (NREGN, IGROUP, MR, MS, ML, NSPR, ILINE, ISIDE,
     &   NLPS, IFLINE, ILLIST, LCON, ISLIST, IFSIDE, LINKR, LINKS,
     &   LINKL, NNN, MAXNL, MXND, LISTL, NUID, NXL, LXN, INDX, NOROOM,
     &   ERR)
C***********************************************************************

C     FXNUID - FIX NUID'S:  RESETS NUID'S OF INTERIOR LINES IN GROUPS
C                           TO ZERO

C***********************************************************************

      DIMENSION IGROUP(NREGN), NSPR(MR), ILINE(ML), ISIDE(MS), NLPS(MS)
      DIMENSION IFLINE(MS), ILLIST(MS*3), LCON(3, ML), ISLIST(MR*4)
      DIMENSION IFSIDE(MR), LINKR(2, MR), LINKS(2, MS), LINKL(2, ML)
      DIMENSION LISTL(MAXNL), NUID(MXND), NXL(2, 3*MXND), LXN(4, MXND)
      DIMENSION INDX(MXND), LINES(20)

      LOGICAL ADDLNK, ERR, LDUP, NOROOM

C  GET LIST OF LINES

      NOROOM = .FALSE.
      ERR = .FALSE.

      ADDLNK = .FALSE.
      N1 = 1
      DO 100 I = 1, NREGN
         CALL LTSORT (MR, LINKR, IGROUP(I), IPTR, ADDLNK)
         IF (IPTR .GT. 0) THEN
            CALL LLIST (MS, ML, MAXNL, NSPR(IPTR), NL, IGROUP(I),
     &         LISTL(N1), ILINE, ISIDE, NLPS, IFLINE, ILLIST, LCON,
     &         ISLIST(IFSIDE(IPTR)), LINKS, LINKL, ERR)
            N1 = N1 + NL
            IF (N1 .GT. MAXNL) THEN
               CALL MESAGE ('IN FXNUID, LINE LIST OVERFLOW')
               NOROOM = .TRUE.
               ERR = .TRUE.
               RETURN
            END IF
         END IF
  100 CONTINUE
      NUML = N1 - 1

C  SORT THE LINE LIST

      IF (NUML .GT. 1) THEN
         DO 110 I = 1, NUML
            INDX(I) = I
  110    CONTINUE
         CALL INDEXI_FQ (NUML, LISTL, NUML, INDX)
      ELSE
         RETURN
      END IF

C  IDENTIFY INTERIOR LINES

      I1 = 1
  120 CONTINUE
      LDUP = .FALSE.
      IF (I1 .LT. NUML) THEN
         I2 = I1 + 1
  130    CONTINUE
         IF (I2 .LE. NUML) THEN
            IF (LISTL(INDX(I1)) .EQ. LISTL(INDX(I2))) THEN
               INDX(I2) = 0
               I2 = I2 + 1
               LDUP = .TRUE.
               GO TO 130
            ELSE
               IF (.NOT.LDUP) INDX(I1) = 0
               I1 = I2
               GO TO 120
            END IF
         END IF
      END IF

C  FORM SORTED LINE LIST IN INDX THEN COPY IT BACK TO LISTL

      N1 = 0
      DO 140 I = 1, NUML
         IF (INDX(I) .GT. 0) THEN
            N1 = N1 + 1
            INDX(N1) = LISTL(INDX(I))
         END IF
  140 CONTINUE
      NUML = N1
      DO 150 I = 1, NUML
         LISTL(I) = INDX(I)
  150 CONTINUE

C  SORT NUID'S ON LINE'S FOR SPEEDY LOOKUP

      N1 = 0
      DO 160 I = 1, NNN
         IF (NUID(I) .GT. 1000000000) THEN
            N1 = N1 + 1
            INDX(N1) = I
         END IF
  160 CONTINUE
      NUMN = N1
      IF (NUMN .GT. 1) CALL INDEXI_FQ (NNN, NUID, NUMN, INDX)

C  LOOP FOR INTERIOR LINES

      DO 220 I = 1, NUML
         KEY = 1000000000 + LISTL(I)*100000

C  FIND LOW POINT

         IBOT = 0
         DO 170 J = 1, NUMN
            IBOT = J
            IF (NUID(INDX(J)) .GE. KEY) GO TO 180
  170    CONTINUE
  180    CONTINUE

C  CHECK INDIVIDUAL POINTS BETWEEN LOW + 1 AND HIGH - 1

         KEY1 = KEY/100000
         DO 190 J = IBOT, NUMN
            IF (NUID(INDX(J))/100000 .EQ. KEY1) THEN
               LXN(2, INDX(J)) = ABS(LXN(2, INDX(J)))
               INDX(J) = 0
            ELSE IF (NUID(INDX(J)) .GT. KEY) THEN
               GO TO 200
            END IF
  190    CONTINUE
  200    CONTINUE

C  COMPACT NUID'S INDEX LIST

         N1 = 0
         DO 210 J = 1, NUMN
            IF (INDX(J) .GT. 0) THEN
               N1 = N1 + 1
               INDX(N1) = INDX(J)
            END IF
  210    CONTINUE
         NUMN = N1

  220 CONTINUE

C  CHECK ALL POINT NUID'S TO MAKE SURE THEY ARE ON BOUNDARY

      DO 240 I = 1, NNN
         IF (NUID(I) .GT. 0 .AND. NUID(I) .LT. 100000) THEN
            NODE = I
            CALL GETLXN (MXND, LXN, NODE, LINES, NL, ERR)
            IF (NL .GT. 20) THEN
               CALL MESAGE ('IN FXNUID, TOO MANY LINES/NODE')
               NOROOM = .TRUE.
               ERR = .TRUE.
               RETURN
            END IF

            KOUNT = 0
            DO 230 J = 1, NL
               I1 = NXL(1, LINES(J)) + NXL(2, LINES(J)) - I
               IF (LXN(2, I1) .LT. 0) KOUNT = KOUNT + 1
  230       CONTINUE
            IF (KOUNT .LT. 2) LXN(2, I) = ABS(LXN(2, I))
         END IF
  240 CONTINUE

      RETURN

      END
