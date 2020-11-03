C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE GETSBC (MXND, MXNPER, NPER, NL, ML, MAXSBC, MAXPRM,
     &   NPRM, NID, LISTL, XN, YN, NUID, LXK, KXL, NXL, LSTSBC, NPERIM,
     &   KSBC, LCON, ISBOUN, LINKL, NSPF, IFSB, LISTSB, LINKSB, LLL,
     &   BAR, ERR)
C***********************************************************************

C  SUBROUTINE GETSBC = GETS THE SIDE BOUNDARY LIST

C***********************************************************************

      DIMENSION NID (MXNPER, MAXPRM), NPERIM (MAXPRM)
      DIMENSION LISTL (NL), XN (MXND), YN (MXND), NUID (MXND)
      DIMENSION LXK (4, MXND), KXL (2, MXND*3), NXL (2, MXND*3)
      DIMENSION LCON (3, ML), ISBOUN (ML), LINKL (2, ML)
      DIMENSION NSPF (ML), IFSB (ML), LISTSB (2, ML), LINKSB (2, ML)
      DIMENSION LSTSBC (MAXSBC), NODES (4)

      LOGICAL EXTER, ERR, CCW, BAR, ADDLNK

      ERR = .TRUE.
      CCW = .TRUE.
      ADDLNK = .FALSE.
      NPERIM (1) = NPER

      DO 110 II = 1, NPRM
         DO 100 I = 1, NPERIM(II)
            IF (BAR) THEN
               NID (I, II) = IABS (NUID (I))
            ELSE
               IF (NID (I, II) .LT. 0) NID (I, II) = - NID (I, II)
            ENDIF
  100    CONTINUE
  110 CONTINUE

C  SORT THROUGH AND PICK OFF ELEMENTS WITH SIDE BOUNDARY CONDITIONS

      DO 240 I = 1, LLL
         IF (BAR) THEN
            I1 = LXK (1, I)
            I2 = LXK (2, I)
         ELSE
            I1 = NXL (1, I)
            I2 = NXL (2, I)
         ENDIF

C  SEE IF THE LINE IS CLEARLY INTERIOR

         IF (I1 .GT. 0 .AND. I2 .GT. 0) THEN
           if ((NUID (I1) .NE. 0) .AND. (NUID (I2) .NE. 0)) THEN
            LTEST = 0
            EXTER = .FALSE.

C  CHECK AGAINST THE PERIMETER LIST TO SEE IF IT IS TRULY EXTERIOR

            DO 130 JJ  =  1, NPRM
               DO 120 J = 1, NPERIM (JJ)
                  IF (ABS (NUID (I1)) .EQ. NID (J, JJ)) THEN
                     IF (J .EQ. 1) THEN
                        J1 = J + 1
                        J2 = NPERIM(JJ)
                     ELSEIF (J .EQ. NPERIM(JJ)) THEN
                        J1 = J - 1
                        J2 = 1
                     ELSE
                        J1 = J - 1
                        J2 = J + 1
                     ENDIF
                     IF ( (ABS (NUID (I2)) .EQ. NID (J1, JJ)) .OR.
     &                  (ABS (NUID (I2)) .EQ. NID (J2, JJ)) )
     &                  EXTER = .TRUE.
                     GOTO 140
                  ENDIF
  120          CONTINUE
  130       CONTINUE
  140       CONTINUE
            IF (EXTER) THEN

C  FIND THE LINE NUMBER IT BELONGS TO

               IF (ABS (NUID (I1)) .GT. 1000000000) THEN
                  LTEST =  (ABS (NUID (I1)) - 1000000000) / 100000
               ELSEIF (ABS (NUID (I2)) .GT. 1000000000) THEN
                  LTEST =  (ABS (NUID (I2)) - 1000000000) / 100000
               ELSE
                  NSUM = ABS (NUID (I1)) + ABS (NUID (I2))
                  DO 150 J = 1, NL
                     CALL LTSORT (ML, LINKL, LISTL (J), K, ADDLNK)
                     IF ((LCON (1, K) + LCON (2, K)) .EQ. NSUM) THEN
                        IF (( (LCON (1, K) .EQ. ABS (NUID (I1))) .AND.
     +                     (LCON (2, K) .EQ. ABS (NUID (I2)))) .OR.
     +                     ((LCON (1, K) .EQ. ABS (NUID (I2))) .AND.
     +                     (LCON (2, K) .EQ. ABS (NUID (I1))))) THEN
                           LTEST = LISTL (J)
                           GOTO 160
                        ENDIF
                     ENDIF
  150             CONTINUE
  160             CONTINUE
               ENDIF

C  FIND THE ELEMENT BOUNDARY FLAG IF THERE IS ONE

               IF (LTEST.LE.0) THEN
                  CALL MESAGE (' ERROR IN SEARCHING NXL FOR '//
     &               'ELEMENT BCC')
                  RETURN
               ELSE
                  CALL LTSORT (ML, LINKL, LTEST, J, ADDLNK)
                  IF (ISBOUN (J) .GT. 0) THEN
                     IFLAG = ISBOUN (J)

C  CHECK TO MAKE SURE LINE IS LINKED TO FLAG
C  AND GET THE NEXT LINK  (NFLAG)

                     CALL LTSORT (ML, LINKSB, IFLAG, L, ADDLNK)
                     DO 170 JJ = IFSB (L), IFSB (L) + NSPF (L) - 1
                        IF (LISTSB (1, JJ) .LT. 0) THEN
                           CALL MESAGE ('PROBLEMS WITH SIDES IN '//
     &                        'FLAG LIST IN GETSBC')
                        ELSE
                           IF (LISTSB (1, JJ) .EQ. LTEST) THEN
                              NFLAG = LISTSB (2, JJ)
                              GOTO 180
                           ENDIF
                        ENDIF
  170                CONTINUE
                     WRITE (*, 10000)IFLAG
                     RETURN
  180                CONTINUE
                     IF (BAR) THEN
                        NELEM = I
                     ELSE
                        NELEM = KXL (1, I)
                        IF (NELEM .EQ. 0)NELEM = KXL (2, I)
                     ENDIF
                     KSBC = KSBC + 1
                     LSTSBC (KSBC) = - IFLAG
                     KSBC = KSBC + 1
                     if (ksbc .gt. maxsbc) stop 'maxsbc error'
                     LSTSBC (KSBC) = NELEM

C  GET THE CORRECT ELEMENT SIDE

                     IF (BAR) THEN
                        JSIDE = 1
                     ELSE
                        CALL GNXKA (MXND, XN, YN, NELEM, NODES, AREA,
     &                     LXK, NXL, CCW)
                        DO 190 J = 1, 4
                           IF (I1 .EQ. NODES (J)) THEN
                              JP1 = J + 1
                              JM1 = J - 1
                              IF (JP1 .EQ. 5)JP1 = 1
                              IF (JM1 .EQ. 0)JM1 = 4
                              IF (I2 .EQ. NODES (JP1)) THEN
                                 JSIDE = J
                                 GOTO 200
                              ELSEIF (I2 .EQ. NODES (JM1)) THEN
                                 JSIDE = JM1
                                 GOTO 200
                              ENDIF
                           ENDIF
  190                   CONTINUE
                        WRITE (*, 10010)NELEM
                        RETURN
  200                   CONTINUE
                     ENDIF
                     KSBC = KSBC + 1
                     LSTSBC (KSBC) = JSIDE

C  SEE IF ANY MORE FLAGS ARE ATTACHED TO THIS SIDE

  210                CONTINUE
                     IF (NFLAG .GT. 0) THEN

C  CHECK TO MAKE SURE LINE IS LINKED TO FLAG
C  AND GET THE NEXT LINK  (NFLAG)

                        IFLAG = NFLAG
                        CALL LTSORT (ML, LINKSB, IFLAG, L, ADDLNK)
                        DO 220 JJ = IFSB (L), IFSB (L) + NSPF (L)
                           IF (LISTSB (1, JJ) .EQ. LTEST) THEN
                              NFLAG = LISTSB (2, JJ)
                              GOTO 230
                           ENDIF
  220                   CONTINUE
                        WRITE (*, 10000)IFLAG
                        RETURN
  230                   CONTINUE
                        KSBC = KSBC + 1
                        LSTSBC (KSBC) = - IFLAG
                        KSBC = KSBC + 1
                        LSTSBC (KSBC) = NELEM
                        KSBC = KSBC + 1
                        LSTSBC (KSBC) = JSIDE
                        GOTO 210
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
       END IF
  240 CONTINUE

      ERR = .FALSE.
      RETURN

10000 FORMAT (' SIDE BOUNDARY FLAG', I5, ' IS NOT PROPERLY LINKED')
10010 FORMAT (' ERROR FINDING CORRECT BOUNDARY SIDE ON ELEMENT', I5)
      END
