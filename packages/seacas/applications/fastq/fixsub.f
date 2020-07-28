C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE FIXSUB (MXND, NNNOLD, NNN, LLLOLD, LLL, KKKOLD, KKK,
     &   XN, YN, NUID, LXK, KXL, NXL, LXN, INDX, IAVAIL, NAVAIL, FINAL)
C***********************************************************************

C  SUBROUTINE FIXSUB = FIXES THE KXL, LXK, NXL, AND LXN ARRAYS FOR
C                      SUBREGIONS - TAKES OUT DUPLICATE LINES AND NODES

C***********************************************************************

      DIMENSION XN(MXND), YN(MXND), NUID(MXND), INDX(MXND)
      DIMENSION LXK(4, MXND), KXL(2, 3*MXND), NXL(2, 3*MXND)
      DIMENSION LXN(4, MXND)
      DIMENSION LINES(20)

      LOGICAL ERR, FINAL, FOUND, NOROOM

C  GENERATE A LIST OF NODES ON THE PERIMETER IN SORTED ORDER

      NPER = 0
      DO 100 I = 1, NNNOLD
         IF (NUID(I) .NE. 0) THEN
            NPER = NPER + 1
            INDX(NPER) = I
         END IF
  100 CONTINUE
      CALL INDEXI_FQ (NNNOLD, NUID, NPER, INDX)

C  GO THROUGH ALL THE BOUNDARY NODES IN THE LIST CHECKING FOR DUPLICATES

      I = NNNOLD + 1
  110 CONTINUE
      IF (NUID(I) .NE. 0) THEN

C  SEE IF ANOTHER NODE EXISTS WITH THE SAME NUID

         CALL LOWFND (MXND, NUID, NPER, INDX, I, IOLD)

C  IF ANOTHER NODE EXISTS, THEN START CHECKING LINES

         IF (IOLD .GT. 0) THEN
            CALL GETLXN (MXND, LXN, IOLD, LINES, KEND, ERR)

C  CHECK ALL THE LINES ATTACHED TO THE OLD NODE, TO SEE IF THEY ARE
C  THE SAME LINE ATTACHED TO THE NODE BEING CHECKED.

C  IF THE SAME LINE EXISTS, DELETE THE LINE, AND MARK THE NODE FOR
C  LATER DELETION

            DO 200 J = 1, 3
               IF (LXN(J, I) .NE. 0) THEN
                  L = ABS(LXN(J, I))
                  N1 = NXL(1, L) + NXL(2, L) - I
                  FOUND = .FALSE.
                  DO 180 K = 1, KEND
                     LOLD = LINES(K)
                     N2 = NXL(1, LOLD) + NXL(2, LOLD) - IOLD
                     IF ((NUID(N2) .EQ. NUID(N1)) .AND.
     &                  (NUID(N2) .NE. 0)) THEN

C  THE SAME LINE HAS BEEN FOUND - CHANGE REFERENCES TO THE LATEST
C  NODE TO REFERENCES TO THE OLD NODE

                        KXL(2, LOLD) = KXL(1, L)
                        KELEM = KXL(1, L)
                        DO 120 II = 1, 4
                           IF (LXK(II, KELEM) .EQ. L)
     &                        LXK(II, KELEM) = LOLD
  120                   CONTINUE
                        CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL,
     &                     N1, L, NNN, ERR, NOROOM)

C  NOW RENUMBER THE REMAINING LINES IN THE KXL AND NXL ARRAYS

                        LLL = LLL - 1
                        DO 130 II = L, LLL
                           KXL(1, II) = KXL(1, II + 1)
                           KXL(2, II) = KXL(2, II + 1)
                           NXL(1, II) = NXL(1, II + 1)
                           NXL(2, II) = NXL(2, II + 1)
  130                   CONTINUE
                        KXL(1, LLL + 1) = 0
                        KXL(2, LLL + 1) = 0
                        NXL(1, LLL + 1) = 0
                        NXL(2, LLL + 1) = 0

C  NOW RENUMBER ANY REFERENCES TO LINES ABOVE L IN THE LXK AND
C  THE LXN ARRAYS

                        DO 150 II = 1, NNN
                           DO 140 JJ = 1, 3
                              IF (ABS(LXN(JJ, II)) .EQ. L) THEN
                                 LXN(JJ, II) = LOLD
                              ELSE IF (ABS(LXN(JJ, II)) .GT. L) THEN
                                 LXN(JJ, II) = ABS(LXN(JJ, II)) - 1
                              END IF
  140                      CONTINUE
                           IF (LXN(4, II) .EQ. L) THEN
                              LXN(4, II) = LOLD
                           ELSE IF (LXN(4, II) .GT. L) THEN
                              LXN(4, II) = LXN(4, II) - 1
                           END IF
  150                   CONTINUE
                        DO 170 II = KKKOLD + 1, KKK
                           DO 160 JJ = 1, 4
                              IF (LXK(JJ, II) .EQ. L) THEN
                                 LXK(JJ, II) = LOLD
                              ELSE IF (LXK(JJ, II) .GT. L) THEN
                                 LXK(JJ, II) = LXK(JJ, II) - 1
                              END IF
  160                      CONTINUE
  170                   CONTINUE
                        FOUND = .TRUE.
                        GOTO 190
                     END IF
  180             CONTINUE

C  END OF CHECK FOR THE SAME LINE - JUST ADD THE LINE TO THE IOLD NODE
C  IF THERE IS A PLACE FOR THE LINE (I.E. THE MAXIMUM IS FOUR/NODE WITH
C  THIS SCHEME).

  190             CONTINUE
                  IF (.NOT.FOUND) THEN
                     CALL ADDLXN (MXND, LXN, NUID, NAVAIL, IAVAIL,
     &                  IOLD, L, NNN, ERR, NOROOM)
                     IF (NXL(1, L) .EQ. I) THEN
                        NXL(1, L) = IOLD
                     ELSE
                        NXL(2, L) = IOLD
                     END IF
                  END IF
               END IF
  200       CONTINUE

C  ALL THE OLD LINES HAVE BEEN GONE THROUGH - NOW DELETE THE NODE I

            DO 210 J = 1, I - 1
               IF (LXN(4, J) .LT. -I) THEN
                  LXN(4, J) = LXN(4, J) + 1
               END IF
  210       CONTINUE

            NNN = NNN - 1
            DO 230 J = I, NNN
               DO 220 K = 1, 3
                  LXN(K, J) = LXN(K, J + 1)
  220          CONTINUE
               IF (LXN(4, J + 1) .LT. 0) THEN
                  LXN(4, J) = LXN(4, J + 1) + 1
               ELSE
                  LXN(4, J) = LXN(4, J + 1)
               END IF
               XN(J) = XN(J + 1)
               YN(J) = YN(J + 1)
               NUID(J) = NUID(J + 1)
  230       CONTINUE
            LXN(1, NNN + 1) = 0
            LXN(2, NNN + 1) = 0
            LXN(3, NNN + 1) = 0
            LXN(4, NNN + 1) = 0
            NUID(NNN + 1) = 0
            DO 240 J = LLLOLD + 1, LLL
               IF (NXL(1, J) .GE. I)NXL(1, J) = NXL(1, J) - 1
               IF (NXL(2, J) .GE. I)NXL(2, J) = NXL(2, J) - 1
  240       CONTINUE
         ELSE
            I = I + 1
         END IF
      ELSE
         I = I + 1
      END IF
      IF (I .LE. NNN) GO TO 110

C  IF THIS IS THE FINAL SUBREGION TO BE ADDED, THEN FLAG
C  THE LXN ARRAY FOR TRULY EXTERIOR NODES, AND CLEAR THE TEMPORARY
C  NUID'S OF THE SUBREGION ONLY BOUNDARY NODES

      IF (FINAL) THEN
         DO 250 I = 1, NNN
            IF ((ABS(NUID(I)) .GT. 1000000000) .OR.
     &         ((ABS(NUID(I)) .LT. 100000) .AND.
     &         (NUID(I) .NE. 0))) THEN
               LXN(2, I) = -ABS(LXN(2, I))
            ELSE
               NUID(I) = 0
               LXN(2, I) = ABS(LXN(2, I))
            END IF
  250    CONTINUE

C  LINK-UP AVAILABLE LXN SPACE

         IAVAIL = NNN + 1
         NAVAIL = MXND - NNN
         DO 260 I = IAVAIL, MXND
            LXN(1, I) = 0
            LXN(2, I) = 0
            LXN(3, I) = 0
            LXN(4, I) = I + 1
  260    CONTINUE
      END IF
      RETURN

      END
