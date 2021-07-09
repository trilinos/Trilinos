C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE BISECT (MP, ML, MS, MR, MSNAP, N, IPOINT, COOR, IPBOUN,
     &   ILINE, LTYPE, NINT, FACTOR, LCON, ILBOUN, ISBOUN, ISIDE, NLPS,
     &   IFLINE, ILLIST, IREGN, IMAT, NSPR, IFSIDE, ISLIST, IRPB, IPBF,
     &   NPPF, IFPB, LISTPB, ILBF, NLPF, IFLB, LISTLB, ISBF, NSPF, IFSB,
     &   LISTSB, LINKP, LINKL, LINKS, LINKR, LINKM, LINKPB, LINKLB,
     &   LINKSB, IHOLDP, IHOLDL, IHOLDS, IHOLDR, IHOLDM, IHOLD2, IHOLD3,
     &   IWTPBF, IWTLBF, IWTSBF, IRGFLG, X, Y, MERGE, LASTP, ALPHA,
     &   NOROOM, SNAP, SNAPDX, NSNAP)
C***********************************************************************

C  SUBROUTINE BISECT = FINDS CLOSEST PERPENDICULAR BISECTED LINE

C***********************************************************************

C  SUBROUTINE CALLED BY:
C     INPUT  = INPUTS MESH DEFINITIONS FROM THE LIGHT TABLE

C***********************************************************************

C  SUBROUTINES CALLED:
C     DLPARA = DETERMINES LINE PARAMETERS FROM TWO POINTS

C***********************************************************************

      DIMENSION IPOINT(MP), COOR(2, MP), IPBOUN(MP)
      DIMENSION ILINE(ML), LTYPE(ML), NINT(ML), FACTOR(ML), LCON(3, ML)
      DIMENSION ILBOUN(ML), ISBOUN(ML)
      DIMENSION ISIDE(MS), NLPS(MS), IFLINE(MS), ILLIST(MS*3)
      DIMENSION IREGN(MR), IMAT(MR), NSPR(MR), IFSIDE(MR), ISLIST(MR*4)
      DIMENSION IRPB(MR)
      DIMENSION IPBF(MP), NPPF(MP), IFPB(MP), LISTPB(2, MP)
      DIMENSION ILBF(ML), NLPF(ML), IFLB(ML), LISTLB(2, ML)
      DIMENSION ISBF(ML), NSPF(ML), IFSB(ML), LISTSB(2, ML)
      DIMENSION IWTPBF(3, MP), IWTLBF(3, ML), IWTSBF(3, ML)
      DIMENSION LINKP(2, MP), LINKL(2, ML), LINKS(2, MS), LINKR(2, MR)
      DIMENSION LINKM(2, (MS + MR)), LINKPB(2, MP), LINKLB(2, ML)
      DIMENSION LINKSB(2, ML)
      DIMENSION IHOLDP(2*MP), IHOLDL(ML*2), IHOLDR(2, MR), IHOLD2(2, ML)
      DIMENSION IHOLDM(2, (MS + MR)), IHOLD3(2, ML), IHOLDS(2, MS)
      DIMENSION IRGFLG(MR)
      DIMENSION N(29), SNAPDX(2, MSNAP), NSNAP(2)

      LOGICAL BIFIND, ADDLNK, NOROOM, ALPHA, TEST, MERGE, NUMPLT, NEWNUM
      LOGICAL ADDCEN, CHANGE, ADJUST, ADDOLD, SIDEOK, GETMAX, SNAP

      NUMPLT = .FALSE.

C  FIND THE CLOSEST LINE ABOVE THE POINT INPUT

      ADDLNK = .FALSE.
      TEST = .FALSE.
      GETMAX = .FALSE.
      N3OLD = N(3)
      CALL CLOSEL (MP, ML, N, COOR, ILINE, LTYPE, LCON, LINKP, LINKL,
     &   X, Y, BIFIND, IFIND, ADDCEN, XCHOLD, YCHOLD)

C  INPUT THE NEW POINT

      IF (BIFIND) THEN
         CALL LTSORT (ML, LINKL, IFIND, L1, ADDLNK)
         LT = LTYPE(L1)
         KNUM = ILINE(L1)
         CALL LINEPR (ML, MP, LINKP, LCON, L1, I1, I2, I3, J1, J2, J3)
         IF (.NOT.ALPHA) THEN
            CALL PLTSTD (1, 0.)
            CALL DLINE (MP, ML, COOR, LINKP, KNUM, LT, I1, I2, I3,
     &         NUMPLT, DUM1, DUM2, TEST, GETMAX, DUM3, DUM4, DUM5, DUM6)
            CALL PLTSTD (1, 7.)
            CALL PLTFLU
         END IF
CC
CC  ADJUST THE NEW POINT TO A SNAP GRID IF APPROPRIATE
CC
C         IF (SNAP) THEN
C            XTEST = X
C            YTEST = Y
C            CALL SNAPPT (MSNAP, SNAPDX, NSNAP, XTEST, YTEST)
CC
CC  CHECK TO MAKE SURE IT DOESN'T SNAP TO THE LINE END POINTS OR THE
CC  CENTER POINT OF THE ARC
CC
C            DIVX = MAX (ABS(XTEST), ABS(COOR(1, J1)), ABS(COOR(1, J2)),
C     &         ABS(COOR(1, J3)))
C            IF (DIVX .EQ. 0.0) DIVX = 1.0
C            DIVY = MAX (ABS(YTEST), ABS(COOR(2, J1)), ABS(COOR(2, J2)),
C     &         ABS(COOR(2, J3)))
C            IF (DIVY .EQ. 0.0) DIVY = 1.0
C            IF (.NOT.(((ABS(XTEST - COOR(1, J1))/DIVX .LT. 1.0E-6) .AND.
C     &         (ABS(YTEST - COOR(2, J1))/DIVY .LT. 1.0E-6))
C     &         .OR.
C     &         ((ABS(XTEST - COOR(1, J2))/DIVX .LT. 1.0E-6) .AND.
C     &         (ABS(YTEST - COOR(2, J2))/DIVY .LT. 1.0E-6))
C     &         .OR.
C     &         ((ABS(XTEST - COOR(1, J3))/DIVX .LT. 1.0E-6) .AND.
C     &         (ABS(YTEST - COOR(2, J3))/DIVY .LT. 1.0E-6)) ) ) THEN
C               X = XTEST
C               Y = YTEST
C            END IF
C         END IF
         JJ = N(18) + 1
         CALL INPOIN (MP, N(1), N(18), JJ, X, Y, NHOLDP, IHOLDP, IPOINT,
     &      COOR, IPBOUN, LINKP, MERGE, NOROOM)
         IF (NOROOM)RETURN
         LASTP = JJ
         CALL LTSORT (MP, LINKP, LASTP, I, ADDLNK)
         IF (.NOT.ALPHA) THEN
            CALL PLTSTD (1, 3.)
#if NeedsDoubleEscape
            CALL MPD2SY (1, COOR(1, I), COOR(2, I), '\\CX')
#else
            CALL MPD2SY (1, COOR(1, I), COOR(2, I), '\CX')
#endif
            CALL PLTSTD (1, 7.)
            CALL PLTFLU
         END IF

C  INPUT THE NEW CENTER IF NEEDED

         IF (ADDCEN) THEN
            JJ = N(18) + 1
            CALL INPOIN (MP, N(1), N(18), JJ, XCHOLD, YCHOLD, NHOLDP,
     &         IHOLDP, IPOINT, COOR, IPBOUN, LINKP, MERGE, NOROOM)
            IF (NOROOM)RETURN
            CALL LTSORT (MP, LINKP, JJ, I, ADDLNK)
            IF (.NOT.ALPHA) THEN
               CALL PLTSTD (1, 3.)
#if NeedsDoubleEscape
               CALL MPD2SY (1, COOR(1, I), COOR(2, I), '\\CX')
#else
               CALL MPD2SY (1, COOR(1, I), COOR(2, I), '\CX')
#endif
               CALL PLTSTD (1, 7.)
               CALL PLTFLU
            END IF
            I3 = JJ
         END IF

C  CALCULATE HOW TO DIVIDE UP THE ORIGINAL INERVALS

         IF (NINT(L1) .EQ. 0) THEN
            NN1 = 0
            NN2 = 0
         ELSE
            D1 = SQRT((COOR(1, J1) - COOR(1, J3))**2 +
     &         (COOR(2, J1) - COOR(2, J3))**2)
            D2 = SQRT((COOR(1, J2) - COOR(1, J3))**2 +
     &         (COOR(2, J2) - COOR(2, J3))**2)
            NN1 = (DBLE(NINT(L1))*D1)/(D1 + D2)
            NN1 = MAX0(1, NN1)
            NN2 = MAX0(1, NINT(L1) - NN1)
         END IF

C  INPUT THE ORIGINAL LINE

         FACT = FACTOR(L1)
         ILB = ILBOUN(L1)
         ISB = ISBOUN(L1)
         KNUM = 0
         IF (LT .GT. 3)LT = 3
         CALL INLINE (ML, N(2), N(19), IFIND, LT, I1, LASTP, I3, NN1,
     &      FACT, NHOLDL, IHOLDL, ILINE, LTYPE, NINT, FACTOR, LCON,
     &      ILBOUN, ISBOUN, LINKL, MERGE, NOROOM)
         IF (NOROOM)RETURN

C  PLOT THE LINE

         IF (.NOT.ALPHA) THEN
            CALL DLINE (MP, ML, COOR, LINKP, KNUM, LT, I1, LASTP, I3,
     &         NUMPLT, DUM1, DUM2, TEST, GETMAX, DUM3, DUM4, DUM5, DUM6)
            CALL PLTFLU
         END IF

C  INPUT THE NEW LINE

         JJ = N(19) + 1
         CALL INLINE (ML, N(2), N(19), JJ, LT, LASTP, I2, I3, NN2, FACT,
     &      NHOLDL, IHOLDL, ILINE, LTYPE, NINT, FACTOR, LCON, ILBOUN,
     &      ISBOUN, LINKL, MERGE, NOROOM)
         IF (NOROOM)RETURN

C  INPUT THE NEW LINE'S FLAGS AND THEN LINK UP FLAGS FOR BOTH LINES

         SIDEOK = .TRUE.
         NHOLD2 = 0
         IF (ILB .NE. 0) THEN
            ADDOLD = .TRUE.
            CALL INBOUN (ML, ILB, 1, N(19), N(26), N(13), N(14), IDUM,
     &         MERGE, NOROOM, NEWNUM, NHOLD2, IHOLD2, ILBF, NLPF, IFLB,
     &         LISTLB, LINKLB, IWTLBF, JHOLD, ADDOLD)
            IF (NOROOM)RETURN
            SIDEOK = .TRUE.
            CALL LINKBC (ML, MS, 1, N(13), N(2), N(26), N(13), N(14),
     &         N(20), ILBF, IFLB, NLPF, LISTLB, NLPS, IFLINE, ILLIST,
     &         ILBOUN, LINKLB, IWTLBF, LINKL, LINKS, SIDEOK, NOROOM)
            IF (NOROOM)RETURN
         END IF
         NHOLD3 = 0
         IF (ISB .NE. 0) THEN
            ADDOLD = .TRUE.
            CALL INBOUN (ML, ISB, 1, N(19), N(27), N(15), N(16), IDUM,
     &         MERGE, NOROOM, NEWNUM, NHOLD3, IHOLD3, ISBF, NSPF, IFSB,
     &         LISTSB, LINKSB, IWTSBF, JHOLD, ADDOLD)
            SIDEOK = .TRUE.
            CALL LINKBC (ML, MS, 1, N(15), N(2), N(27), N(15), N(16),
     &         N(20), ISBF, IFSB, NSPF, LISTSB, NLPS, IFLINE, ILLIST,
     &         ISBOUN, LINKSB, IWTSBF, LINKL, LINKS, SIDEOK, NOROOM)
            IF (NOROOM)RETURN
         END IF

C  PLOT THE NEW LINE

         IF (.NOT.ALPHA) THEN
            CALL DLINE (MP, ML, COOR, LINKP, KNUM, LT, LASTP, I2, I3,
     &         NUMPLT, DUM1, DUM2, TEST, GETMAX, DUM3, DUM4, DUM5, DUM6)
            CALL PLTFLU
         END IF

C  UPDATE ANY REGION OR SIDE DEFINITIONS TO INCLUDE THE NEW LINES

         IEND = N(22)
         DO 120 I = 1, IEND
            CALL LTSORT (MR, LINKR, I, II, ADDLNK)
            IF (II .GT. 0) THEN
               ADJUST = .FALSE.
               KKOUNT = 0

C  GET THE MATCH POINTS FROM THE LAST LINE IN THE REGION

               JBEGIN = IFSIDE(II)
               JEND = IFSIDE(II) + NSPR(II) - 1
               CALL LTSORT (MS, LINKS, ISLIST(JEND), IPNTR, ADDLNK)
               IF ((ISLIST(JEND) .GT. 0) .AND. (IPNTR .GT. 0)) THEN
                  KFIRST = IFLINE(IPNTR) + NLPS(IPNTR) - 1
                  CALL LTSORT (ML, LINKL, ILLIST(KFIRST), KK, ADDLNK)
                  IF (KK .GT. 0.) THEN
                     J1END = LCON(1, KK)
                     J2END = LCON(2, KK)
                  ELSE
                     J1END = 0
                     J2END = 0
                  END IF
                  KEND = IFLINE(IPNTR) + NLPS(IPNTR) - 1
                  CALL LTSORT (ML, LINKL, ILLIST(KEND), KK, ADDLNK)
                  IF (KK .GT. 0.) THEN
                     I1END = LCON(1, KK)
                     I2END = LCON(2, KK)
                  ELSE
                     I1END = 0
                     I2END = 0
                  END IF
               ELSE IF (ISLIST(JEND) .LT. 0) THEN
                  JJ = IABS(ISLIST(JEND))
                  CALL LTSORT (ML, LINKL, JJ, KK, ADDLNK)
                  IF (KK .GT. 0.) THEN
                     I1END = LCON(1, KK)
                     I2END = LCON(2, KK)
                  ELSE
                     I1END = 0
                     I2END = 0
                  END IF
                  J1END = 0
                  J2END = 0
               ELSE
                  I1END = 0
                  I2END = 0
                  J1END = 0
                  J2END = 0
               END IF

C  NOW LOOP THROUGH THE SIDES/LINES AND CHECK FOR REPLACEMENT

               DO 110 J = JBEGIN, JEND
                  KKOUNT = KKOUNT + 1
                  IHOLDL(KKOUNT) = ISLIST(J)
                  CALL LTSORT (MS, LINKS, ISLIST(J), IPNTR, ADDLNK)
                  IF ((ISLIST(J) .GT. 0) .AND. (IPNTR .GT. 0) .AND.
     &               (ISIDE(IPNTR) .GT. 0)) THEN

C  CHECK THE SIDE FOR LINE INCLUSION IF THE SIDE EXISTS, AND IT
C  HAS NOT ALREADY BEEN FIXED (ISIDE(IPNTR)  .LT.  0)

                     JJ = ISLIST(J)
                     KOUNT = 0
                     CHANGE = .FALSE.
                     KBEGIN = IFLINE(IPNTR)
                     KEND = IFLINE(IPNTR) + NLPS(IPNTR) - 1

C  GET THE BEGINNING OF THE SIDE FOR SWITCHING

                     CALL LTSORT (ML, LINKL, ILLIST(KBEGIN), KK, ADDLNK)
                     IF (KK .GT. 0.) THEN
                        J1END = LCON(1, KK)
                        J2END = LCON(2, KK)
                     ELSE
                        J1END = 0
                        J2END = 0
                     END IF
                     DO 100 K = KBEGIN, KEND
                        KOUNT = KOUNT + 1
                        IHOLDP(KOUNT) = ILLIST(K)
                        CALL LTSORT (ML, LINKL, ILLIST(K), KK, ADDLNK)
                        IF ((KK .GT. 0.) .AND.
     &                     (ILLIST(K) .EQ. IFIND)) THEN
                           CHANGE = .TRUE.
                           KOUNT = KOUNT + 1
                           IHOLDP(KOUNT) = N(19)

C  CHECK TO MAKE SURE THE TWO LINES DON'T NEED SWITCHING

                           IF ((I2 .EQ. I1END) .OR.
     &                        (I2 .EQ. I2END)) THEN
                              I1END = I1
                              I2END = LASTP
                              IHOLDP(KOUNT) = IHOLDP(KOUNT - 1)
                              IHOLDP(KOUNT - 1) = N(19)
                           ELSE
                              I1END = LASTP
                              I2END = I2
                           END IF
                        ELSE IF (KK .GT. 0) THEN
                           I1END = LCON(1, KK)
                           I2END = LCON(2, KK)
                        ELSE
                           I1END = 0
                           I2END = 0
                        END IF
  100                CONTINUE

C  INCLUDE THE LINE IN THE SIDE LIST IF NEEDED

                     IF (CHANGE) THEN
                        ID = JJ
                        CALL INSIDE (MS, N(3), N(4), N(20), JJ, IHOLDP,
     &                     KOUNT, ISIDE, NLPS, IFLINE, ILLIST, LINKS,
     &                     NHOLDS, IHOLDS, MERGE, NOROOM)
                        IF (NOROOM)RETURN

C  FLAG THE SIDE DATA AS ALREADY FIXED

                        CALL LTSORT (MS, LINKS, ID, KK, ADDLNK)
                        IF (KK .GT. 0)ISIDE(KK) = -IABS(ISIDE(KK))
                     END IF

C  CHECK THE LINE FOR INCLUSION

                  ELSE IF (ISLIST(J) .LT. 0) THEN
                     JJ = IABS(ISLIST(J))
                     CALL LTSORT (ML, LINKL, JJ, KK, ADDLNK)
                     IF ((KK .GT. 0.) .AND. (JJ .EQ. IFIND)) THEN
                        ADJUST = .TRUE.
                        KKOUNT = KKOUNT + 1
                        IHOLDL(KKOUNT) = -N(19)

C  CHECK TO MAKE SURE THE TWO LINES DON'T NEED SWITCHING

                        IF ((I2 .EQ. I1END) .OR. (I2 .EQ. I2END) .OR.
     &                     (I2 .EQ. J1END) .OR. (I2 .EQ. J2END)) THEN
                           I1END = I1
                           I2END = LASTP
                           IHOLDL(KKOUNT) = IHOLDL(KKOUNT - 1)
                           IHOLDL(KKOUNT - 1) = -N(19)
                        ELSE
                           I1END = LASTP
                           I2END = I2
                        END IF
                     ELSE IF (KK .GT. 0) THEN
                        I1END = LCON(1, KK)
                        I2END = LCON(2, KK)
                     ELSE
                        I1END = 0
                        I2END = 0
                     END IF
                     J1END = 0
                     J2END = 0
                  END IF
  110          CONTINUE
               IF (ADJUST) THEN
                  IMTRL = IMAT(II)
                  ID = I
                  CALL INREGN (MS, MR, N(7), N(8), N(22), N(23), ID,
     &               IMTRL, IHOLDL, KKOUNT, IREGN, IMAT, NSPR, IFSIDE,
     &               ISLIST, LINKR, LINKM, NHOLDR, IHOLDR, NHOLDM,
     &               IHOLDM, IRGFLG, MERGE, NOROOM)
                  IF (NOROOM)RETURN
               END IF
            END IF
  120    CONTINUE

C  RESTORE ALL THE FLAGGED SIDE DATA

         DO 130 I = N3OLD + 1, N(3)
            ISIDE(I) = IABS(ISIDE(I))
  130    CONTINUE
      ELSE
         LASTP = 0
      END IF
      RETURN
      END
