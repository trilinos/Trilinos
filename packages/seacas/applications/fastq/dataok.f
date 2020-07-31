C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE DATAOK (MP, ML, MS, MR, L, KNUM, COOR, ILINE, LTYPE,
     &   KNINT, LCON, NLPS, IFLINE, ILLIST, NSPR, IFSIDE, ISLIST, LINKP,
     &   LINKL, LINKS, SIZE, ERRCHK, ERR)
C***********************************************************************

C  SUBROUTINE FILLOK = CHECKS TO MAKE SURE NONEXISTENT DATA IS NOT
C                      BEING REFERENCED IN THE REGION DEFINITIONS

C***********************************************************************

      DIMENSION COOR (2, MP), LINKP (2, MP)
      DIMENSION ILINE (ML), LTYPE (ML), KNINT (ML), LCON (3, ML)
      DIMENSION LINKL (2, ML)
      DIMENSION NLPS (MS), IFLINE (MS), ILLIST (MS*3), LINKS (2, MS)
      DIMENSION NSPR (MR), IFSIDE (MR), ISLIST (MR*4)

      LOGICAL ERR, ADDLNK, ERRCHK

      ERR = .TRUE.
      ADDLNK = .FALSE.

      DO 130 I = IFSIDE (L), IFSIDE (L) + NSPR (L)-1

C  CHECK TO MAKE SURE REGION'S SIDE DEFINITIONS ARE ALL THERE

         IF (ISLIST (I).GT.0)THEN
            II = ISLIST (I)
            CALL LTSORT (MS, LINKS, II, IPNTR, ADDLNK)
            IF (IPNTR.LE.0) THEN
               IF (ERRCHK) THEN
                  WRITE (*, 10000)KNUM, II
                  RETURN
               ELSE
                  GOTO 120
               ENDIF
            END IF

C  CHECK TO MAKE SURE SIDE'S LINE DEFINITIONS ARE ALL THERE

            CALL LTSORT (MS, LINKS, II, JJ, ADDLNK)
            DO 110 J = IFLINE (JJ), IFLINE (JJ) + NLPS (JJ)-1
               KK = ILLIST (J)
               CALL LTSORT (ML, LINKL, KK, LL, ADDLNK)
               IF ((KK.LE.0) .OR. (LL.LE.0)) THEN
                  IF (ERRCHK) THEN
                     WRITE (*, 10010)II, KK
                     RETURN
                  ELSE
                     GOTO 100
                  ENDIF
               END IF

C  CHECK TO MAKE SURE LINE'S POINT DEFINITIONS ARE ALL THERE

               I1 = LCON (1, LL)
               I2 = LCON (2, LL)
               I3 = LCON (3, LL)
               CALL LTSORT (MP, LINKP, I1, J1, ADDLNK)
               CALL LTSORT (MP, LINKP, I2, J2, ADDLNK)
               IF (I3 .NE. 0) THEN
                  CALL LTSORT (MP, LINKP, IABS (I3), J3, ADDLNK)
               ELSE
                  J3 = 0
               END IF

               IF ((I1.LE.0) .OR. (J1.LE.0)) THEN
                  IF (ERRCHK) THEN
                     WRITE (*, 10030)KK, I1
                     RETURN
                  ELSE
                     GOTO 100
                  ENDIF
               ELSEIF ((I2.LE.0) .OR. (J2.LE.0)) THEN
                  IF (ERRCHK) THEN
                     WRITE (*, 10030)KK, I2
                     RETURN
                  ELSE
                     GOTO 100
                  ENDIF
               ELSEIF ((LTYPE (LL) .NE. 1) .AND. ((I3 .EQ. 0) .OR.
     &            (J3.LE.0))) THEN
                  IF (ERRCHK) THEN
                     WRITE (*, 10030)KK, I3
                     RETURN
                  ELSE
                     GOTO 100
                  ENDIF
               END IF

C  CHECK TO INSURE AN INTEGRAL ASSIGNMENT

               IF (IABS (KNINT (LL)) .EQ. 0) THEN
                  IF (I3 .LT. 0)J3 = -J3
                  CALL LINLEN (MP, COOR, LINKP, KNUM, ILINE(LL),
     &               LTYPE(LL), I3, J1, J2, J3, DIST, ERR)
                  IF (ERR) THEN
                     IF (ERRCHK) THEN
                        WRITE (*, 10020)KK, IABS (KNINT (LL))
                        RETURN
                     ELSE
                        GOTO 100
                     ENDIF
                  ELSE
                     IF (SIZE .LE. 0.) THEN
                        KNINT (LL) = 1
                     ELSE
                        KNINT (LL) = MAX0 (NINT (DIST/SIZE), 1)
                     END IF
                  END IF
               END IF
  100          CONTINUE
  110       CONTINUE

C  CHECK TO MAKE SURE REGION'S LINE DEFINITIONS ARE ALL THERE

         ELSEIF (ISLIST (I) .LT. 0) THEN
            KK = IABS (ISLIST (I))
            CALL LTSORT (ML, LINKL, KK, LL, ADDLNK)
            IF ( (KK .LE. 0)  .OR. (LL .LE. 0) ) THEN
               IF (ERRCHK) THEN
                  WRITE (*, 10010)KNUM, KK
                  RETURN
               ELSE
                  GOTO 120
               ENDIF
            END IF

C  CHECK TO MAKE SURE LINE'S POINT DEFINITIONS ARE ALL THERE

            I1 = LCON (1, LL)
            I2 = LCON (2, LL)
            I3 = LCON (3, LL)
            CALL LTSORT (MP, LINKP, I1, J1, ADDLNK)
            CALL LTSORT (MP, LINKP, I2, J2, ADDLNK)
            IF (I3 .NE. 0) THEN
               CALL LTSORT (MP, LINKP, IABS (I3), J3, ADDLNK)
            ELSE
               J3 = 0
            END IF

            IF ((I1.LE.0) .OR. (J1.LE.0)) THEN
               IF (ERRCHK) THEN
                  WRITE (*, 10030)KK, I1
                  RETURN
               ELSE
                  GOTO 120
               ENDIF
            ELSEIF ((I2.LE.0) .OR. (J2.LE.0)) THEN
               IF (ERRCHK) THEN
                  WRITE (*, 10030)KK, I2
                  RETURN
               ELSE
                  GOTO 120
               ENDIF
            ELSEIF ((LTYPE (LL) .NE. 1) .AND. ((I3 .EQ. 0) .OR.
     &         (J3.LE.0))) THEN
               IF (ERRCHK) THEN
                  WRITE (*, 10030)KK, I3
                  RETURN
               ELSE
                  GOTO 120
               ENDIF
            END IF

C  CHECK TO MAKE SURE INTERVAL ASSIGNMENT IS HANDLED

            IF (IABS (KNINT (LL)) .EQ. 0) THEN

C**MBS/29-JUN-1989/ DO NOT NEGATE POINTER TO CENTER OF CLOCKWISE ARC
C              IF (I3 .LT. 0)J3 = -J3
               CALL LINLEN (MP, COOR, LINKP, KNUM, ILINE(LL),
     &            LTYPE(LL), I3, J1, J2, J3, DIST, ERR)
               IF (ERR) THEN
                  IF (ERRCHK) THEN
                     WRITE (*, 10020)KK, IABS (KNINT (LL))
                     RETURN
                  ELSE
                     GOTO 120
                  ENDIF
               ELSE
                  IF (SIZE .LE. 0.) THEN
                     KNINT (LL) = 1
                  ELSE
                     KNINT (LL) = MAX0 (NINT (DIST/SIZE), 1)
                  END IF
               END IF
            END IF

C  A ZERO SIDE NUMBER HAS BEEN FOUND

         ELSE
            IF (ERRCHK) THEN
               WRITE (*, 10000)KNUM, ISLIST (I)
            ELSE
               GOTO 120
            ENDIF
         END IF
  120    CONTINUE
  130 CONTINUE

C  ALL DEFINITIONS ARE IN ORDER

      ERR = .FALSE.
      RETURN

10000 FORMAT (' FOR REGION:', I5, ' SIDE:', I5, ' DOES NOT EXIST')
10010 FORMAT (' FOR SIDE:', I5, ' LINE:', I5, ' DOES NOT EXIST')
10020 FORMAT (' FOR LINE:', I5, ' INTERVAL OF:', I5, ' IS NOT WORKING')
10030 FORMAT (' FOR LINE:', I5, ' POINT:', I5, ' DOES NOT EXIST')

      END
