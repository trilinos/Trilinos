C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE WRFSQ (IUNIT, MP, ML, MS, MR, MSNAP, MSC, MCOM, ICOM,
     &   JCOM, CIN, RIN, IIN, KIN, N, IPOINT, COOR, IPBOUN, ILINE,
     &   LTYPE, NINT, FACTOR, LCON, ILBOUN, ISBOUN, ISIDE, NLPS,
     &   IFLINE, ILLIST, IBARST, JMAT, JCENT, NLPB, JFLINE,
     &   JLLIST, IREGN, IMAT, NSPR, IFSIDE, ISLIST, IRPB, IPBF, NPPF,
     &   IFPB, LISTPB, ILBF, NLPF, IFLB, LISTLB, ISBF, NSPF, IFSB,
     &   LISTSB, LINKP, LINKL, LINKS, LINKB, LINKR, LINKSC, LINKPB,
     &   LINKLB, LINKSB, IWTPBF, IWTLBF, IWTSBF, RSIZE, IFHOLE, NHPR,
     &   IHLIST, IRGFLG, ISCHM, SCHEME, NUMBER, DEFSCH, DEFSIZ, TITLE,
     &   OPTIM, THREE, EIGHT, NINE, SNAP, SNAPDX, NSNAP, REGWRT, BARWRT)
C***********************************************************************

C  SUBROUTINE WRFSQ  =  WRITES FASTQ CARD FILE

C***********************************************************************

C  SUBROUTINE CALLED BY:
C     FASTQ  =  A PROGRAM TO QUICKLY PREPARE FASTQ INPUT

C***********************************************************************

      DIMENSION IPOINT(MP), COOR(2, MP), IPBOUN(MP)
      DIMENSION ILINE(ML), LTYPE(ML), NINT(ML), FACTOR(ML), LCON(3, ML)
      DIMENSION ILBOUN(ML), ISBOUN(ML)
      DIMENSION ISIDE(MS), NLPS(MS), IFLINE(MS), ILLIST(MS*3)
      DIMENSION IBARST(MS), JMAT(MS), JCENT(MS), NLPB(MS), JFLINE(MS)
      DIMENSION JLLIST(MS*3)
      DIMENSION IREGN(MR), IMAT(MR), NSPR(MR), IFSIDE(MR), ISLIST(MR*4)
      DIMENSION IRPB(MR), RSIZE(MR), IFHOLE(MR), NHPR(MR), IHLIST(MR*2)
      DIMENSION ISCHM(MSC), SCHEME(MSC), IRGFLG(MR)
      DIMENSION IPBF(MP), NPPF(MP), IFPB(MP), LISTPB(2, MP)
      DIMENSION ILBF(ML), NLPF(ML), IFLB(ML), LISTLB(2, ML)
      DIMENSION ISBF(ML), NSPF(ML), IFSB(ML), LISTSB(2, ML)
      DIMENSION IWTPBF(3, MP), IWTLBF(3, MP), IWTSBF(3, MP)
      DIMENSION LINKP(2, MP), LINKL(2, ML), LINKS(2, MS), LINKB(2, MS)
      DIMENSION LINKR(2, MR), LINKSC(2, MR)
      DIMENSION LINKPB(2, MP), LINKLB(2, ML), LINKSB(2, ML)
      DIMENSION NUMBER(MSC), SNAPDX(2, MSNAP), NSNAP(2)
      DIMENSION KIN(MCOM), IIN(MCOM), RIN(MCOM)
      DIMENSION N(29), ID (13)

      CHARACTER*72 SCHEME, DEFSCH, TITLE, DUMMY, CIN(MCOM)
      CHARACTER NUMBER*80, TYPE(7)*5

      LOGICAL IANS, OPTIM, THREE, EIGHT, NINE, ADDLNK, SNAP
      LOGICAL REGWRT, BARWRT, FLAG, GOWRIT, STAR

      DATA TYPE/'  STR', ' CORN', ' CIRC', ' CIRM', ' PARA', ' CIRR',
     &   ' ELIP'/

      ADDLNK = .FALSE.
      GOWRIT = .FALSE.
      XADD = 0.
      YADD = 0.

C  WRITE OUT ONLY THE REGIONS OF INTEREST IF THE REGWRT FLAG HAS BEEN SET

      IF (REGWRT) THEN

C  SEE IF A SHIFT OF THE REGION IS NEEDED

         CALL INTRUP ('SHIFT REGION', IANS, MCOM, ICOM, JCOM, CIN, IIN,
     &      RIN, KIN)
         IF (IANS) THEN
            IZ = 0
            CALL FREFLD (IZ, IZ, 'X SHIFT: ', MCOM, IOSTAT, IFOUND, KIN,
     &         CIN, IIN, RIN)
            IF (IFOUND .GT. 0) THEN
               XADD = RIN (1)
            ELSE
               XADD = 0.
            ENDIF
            CALL FREFLD (IZ, IZ, 'Y SHIFT: ', MCOM, IOSTAT, IFOUND, KIN,
     &         CIN, IIN, RIN)
            IF (IFOUND .GT. 0) THEN
               YADD = RIN (1)
            ELSE
               YADD = 0.
            ENDIF
         ENDIF
         FLAG = .FALSE.
         CALL FLAGD (MP, N (18), LINKP, IPOINT, FLAG)
         CALL FLAGD (ML, N (19), LINKL, ILINE, FLAG)
         CALL FLAGD (MS, N (20), LINKS, ISIDE, FLAG)
         CALL FLAGD (MS, N (21), LINKB, IBARST, FLAG)
         CALL FLAGD (MR, N (22), LINKR, IREGN, FLAG)
         CALL MESAGE ('WRITE REGIONS FROM <I1> TO <I2>')
         CALL MESAGE ('HIT RETURN TO END INPUT')
  100    CONTINUE
         IF (ICOM.GT.JCOM)THEN
            CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN, CIN,
     &         IIN, RIN)
            ICOM = 1
         ENDIF
         CALL GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2, IFOUND)

         IF (IFOUND .GT. 0) THEN
            IF (I1 .GT. 0) THEN
               CALL CHECK (I1, I2, N (22))

C  FLAG ALL DATA ASSOCIATED WITH THE REGIONS

               DO 150 I = I1, I2
                  CALL LTSORT (MR, LINKR, I, II, ADDLNK)
                  IF (II.GT.0) THEN
                     GOWRIT = .TRUE.
                     IREGN (II) = -IABS (IREGN (II))
                     DO 140 J = IFSIDE (II), IFSIDE (II)+NSPR (II)-1

C  FLAG SIDE DATA

                        CALL LTSORT (MS, LINKS, ISLIST (J), JJ, ADDLNK)
                        IF ((ISLIST (J) .GT. 0) .AND. (JJ .GT. 0)) THEN
                           ISIDE (JJ) = -IABS (ISIDE (JJ))
                           DO 120 K = IFLINE (JJ), IFLINE (JJ) +
     &                        NLPS (JJ)-1
                              CALL LTSORT (ML, LINKL, ILLIST (K), KK,
     &                           ADDLNK)
                              IF (KK .GT. 0 )THEN
                                 ILINE (KK) = -IABS (ILINE (KK))
                                 DO 110 L = 1, 3
                                    IF (LCON (L, KK) .GT. 0) THEN
                                       CALL LTSORT (MP, LINKP,
     &                                    LCON (L, KK), LL, ADDLNK)
                                       IF (LL. GT. 0) THEN
                                          IPOINT (LL) =
     &                                       -IABS (IPOINT (LL))
                                       ENDIF
                                    ENDIF
  110                            CONTINUE
                              ENDIF
  120                      CONTINUE

C  FLAG LINE DATA

                        ELSE
                           JJ = IABS (ISLIST (J))
                           CALL LTSORT (ML, LINKL, JJ, KK, ADDLNK)
                           IF (KK.GT.0)THEN
                              ILINE (KK) = -IABS (ILINE (KK))
                              DO 130 L = 1, 3
                                 IF (LCON (L, KK) .GT. 0) THEN
                                    CALL LTSORT (MP, LINKP,
     &                                 LCON (L, KK), LL, ADDLNK)
                                    IF (LL .GT. 0) THEN
                                       IPOINT (LL) =
     &                                    - IABS (IPOINT (LL))
                                    ENDIF
                                 ENDIF
  130                         CONTINUE
                           ENDIF
                        ENDIF
  140                CONTINUE
                  ENDIF
  150          CONTINUE
               GOTO 100
            ENDIF
         ENDIF

C  WRITE OUT THE BARSET DATA

      ELSEIF (BARWRT) THEN
         FLAG = .FALSE.
         CALL FLAGD (MP, N (18), LINKP, IPOINT, FLAG)
         CALL FLAGD (ML, N (19), LINKL, ILINE, FLAG)
         CALL FLAGD (MS, N (20), LINKS, ISIDE, FLAG)
         CALL FLAGD (MS, N (21), LINKB, IBARST, FLAG)
         CALL FLAGD (MR, N (22), LINKR, IREGN, FLAG)
         GOWRIT = .FALSE.
         CALL MESAGE ('WRITE BARSETS FROM <I1> TO <I2>')
         CALL MESAGE ('HIT RETURN TO END INPUT')
  160    CONTINUE
         IF (ICOM .GT. JCOM) THEN
            CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN, CIN,
     &         IIN, RIN)
            ICOM = 1
         ENDIF
         CALL GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2, IFOUND)
         IF (IFOUND .GT. 0) THEN
            IF (I1 .GT. 0) THEN
               CALL CHECK (I1, I2, N (21))

C  FLAG ALL LINES ASSOCIATED WITH THE BARSETS

               DO 190 I = I1, I2
                  CALL LTSORT (MS, LINKB, I, II, ADDLNK)
                  IF (II .GT. 0) THEN
                     GOWRIT = .TRUE.
                     IBARST (II) = -IABS (IBARST (II))
                     DO 180 J = JFLINE (II), JFLINE (II) + NLPB (II) - 1
                        JJ = IABS (JLLIST (J))
                        CALL LTSORT (ML, LINKL, JJ, KK, ADDLNK)
                        IF (KK .GT. 0) THEN
                           ILINE (KK) = -IABS (ILINE (KK))
                           DO 170 L = 1, 3
                              IF (LCON (L, KK) .GT. 0) THEN
                                 CALL LTSORT (MP, LINKP, LCON (L, KK),
     &                              LL, ADDLNK)
                                 IF (LL .GT. 0) THEN
                                    IPOINT (LL) = -IABS (IPOINT (LL))
                                 ENDIF
                              ENDIF
  170                      CONTINUE
                        ENDIF
  180                CONTINUE
                  ENDIF
  190          CONTINUE
               GOTO 160
            ENDIF
         ENDIF

C  OTHERWISE FLAG ALL THE DATA TO BE WRITTEN

      ELSE
         FLAG = .TRUE.
         CALL FLAGD (MP, N (18), LINKP, IPOINT, FLAG)
         CALL FLAGD (ML, N (19), LINKL, ILINE, FLAG)
         CALL FLAGD (MS, N (20), LINKS, ISIDE, FLAG)
         CALL FLAGD (MS, N (21), LINKB, IBARST, FLAG)
         CALL FLAGD (MR, N (22), LINKR, IREGN, FLAG)
         GOWRIT = .TRUE.
      ENDIF

      IF (.NOT. GOWRIT) THEN
         CALL MESAGE ('** NO DATA HAS BEEN WRITTEN **')
         GOTO 510
      ENDIF

C  WRITE OUT THE TITLE

      CALL STRLNG (TITLE, LEN)
      WRITE (IUNIT, 10010) TITLE(1:LEN)

C  WRITE OUT THE POINTS IN ORDER

      DO 200 I = 1, N(18)
         CALL LTSORT (MP, LINKP, I, J, ADDLNK)
         IF ((J .GT. 0) .AND. (IPOINT (J) .LT. 0)) THEN
            WRITE (IUNIT, 10020) IABS (IPOINT(J)),
     &         COOR(1, J) + XADD, COOR(2, J) + YADD
         END IF
  200 CONTINUE

C  WRITE OUT THE LINES IN ORDER

      DO 210 I = 1, N(19)
         CALL LTSORT (ML, LINKL, I, J, ADDLNK)
         IF ((J .GT. 0) .AND. (ILINE (J) .LT. 0)) THEN
            WRITE (IUNIT, 10030) IABS (ILINE(J)), TYPE(LTYPE(J)),
     &         LCON(1, J), LCON(2, J), LCON(3, J), NINT(J), FACTOR(J)
         END IF
  210 CONTINUE

C  WRITE OUT THE SIDES IN ORDER

      DO 230 I = 1, N(20)
         CALL LTSORT (MS, LINKS, I, J, ADDLNK)
         IF ((J .GT. 0) .AND. (ISIDE (J) .LT. 0)) THEN
            N2 = IFLINE(J) - 1
  220       CONTINUE
            N1 = N2 + 1
            IF (N1 .EQ. IFLINE(J)) THEN
               N2 = N1 + 9
            ELSE
               N2 = N1 + 10
            END IF
            N2 = MIN0(N2, IFLINE(J) + NLPS(J) - 1)
            IF (N2 .LT. IFLINE(J) + NLPS(J) - 1) THEN
               IF (N1 .EQ. IFLINE(J)) THEN
                  WRITE (IUNIT, 10040) ' SIDE  ', IABS (ISIDE(J)),
     &               (ILLIST(K), K = N1, N2)
               ELSE
                  WRITE (IUNIT, 10050) (ILLIST(K), K = N1, N2)
               END IF
               GO TO 220
            ELSE
               IF (N1 .EQ. IFLINE(J)) THEN
                  WRITE (IUNIT, 10060) 'SIDE  ', IABS (ISIDE(J)),
     &               (ILLIST(K), K = N1, N2)
               ELSE
                  WRITE (IUNIT, 10070) (ILLIST(K), K = N1, N2)
               END IF
            END IF
         END IF
  230 CONTINUE

C  WRITE OUT THE BAR SETS IN ORDER

      DO 250 I = 1, N(21)
         CALL LTSORT (MS, LINKB, I, J, ADDLNK)
         IF ((J .GT. 0) .AND. (IBARST (J) .LT. 0)) THEN
            N2 = JFLINE(J) - 1
  240       CONTINUE
            N1 = N2 + 1
            IF (N1 .EQ. JFLINE(J)) THEN
               N2 = N1 + 7
            ELSE
               N2 = N1 + 10
            END IF
            N2 = MIN0(N2, JFLINE(J) + NLPB(J) - 1)
            IF (N2 .LT. JFLINE(J) + NLPB(J) - 1) THEN
               IF (N1 .EQ. JFLINE(J)) THEN
                  WRITE (IUNIT, 10040) ' BARSET', IABS (IBARST(J)),
     &               JMAT(J), JCENT(J), (JLLIST(K), K = N1, N2)
               ELSE
                  WRITE (IUNIT, 10050) (JLLIST(K), K = N1, N2)
               END IF
               GO TO 240
            ELSE
               IF (N1 .EQ. JFLINE(J)) THEN
                  WRITE (IUNIT, 10060) 'BARSET', IABS (IBARST(J)),
     &               JMAT(J), JCENT(J), (JLLIST(K), K = N1, N2)
               ELSE
                  WRITE (IUNIT, 10070) (JLLIST(K), K = N1, N2)
               END IF
            END IF
         END IF
  250 CONTINUE

C  WRITE OUT THE REGIONS IN ORDER

      DO 270 I = 1, N(22)
         CALL LTSORT (MR, LINKR, I, J, ADDLNK)
         IF ((J .GT. 0) .AND. (IRGFLG(J) .LE. -1) .AND.
     &      (IREGN (J) .LE. 0)) THEN
            N2 = IFSIDE(J) - 1
  260       CONTINUE
            N1 = N2 + 1
            IF (N1 .EQ. IFSIDE(J)) THEN
               N2 = N1 + 8
            ELSE
               N2 = N1 + 10
            END IF
            N2 = MIN0(N2, IFSIDE(J) + NSPR(J) - 1)
            IF (N2 .LT. IFSIDE(J) + NSPR(J) - 1) THEN
               IF (N1 .EQ. IFSIDE(J)) THEN
                  WRITE (IUNIT, 10040) ' REGION', IABS (IREGN(J)),
     &               IMAT(J), (ISLIST(K), K = N1, N2)
               ELSE
                  WRITE (IUNIT, 10050) (ISLIST(K), K = N1, N2)
               END IF
               GO TO 260
            ELSE
               IF (N1 .EQ. IFSIDE(J)) THEN
                  WRITE (IUNIT, 10060) 'REGION', IABS (IREGN(J)),
     &               IMAT(J), (ISLIST(K), K = N1, N2)
               ELSE
                  WRITE (IUNIT, 10070) (ISLIST(K), K = N1, N2)
               END IF
            END IF

C  WRITE OUT THE REGION INTERVAL SIZE DATA

            IF (RSIZE(J) .GT. 0.)
     &         WRITE (IUNIT, 10080) 'SIZE  ', RSIZE(J),
     &         IABS (IREGN(J))
         END IF
  270 CONTINUE
      IF (DEFSIZ .GT. 0.) WRITE (IUNIT, 10080) 'SIZE  ', DEFSIZ

C  WRITE OUT THE HOLES IN ORDER

      DO 290 I = 1, N(22)
         CALL LTSORT (MR, LINKR, I, J, ADDLNK)
         IF ((NHPR(J) .GT. 0) .AND. (IREGN (J) .LT. 0)) THEN
            N2 = IFHOLE(J) - 1
  280       CONTINUE
            N1 = N2 + 1
            IF (N1 .EQ. IFHOLE(J)) THEN
               N2 = N1 + 9
            ELSE
               N2 = N1 + 10
            END IF
            N2 = MIN0(N2, IFHOLE(J) + NHPR(J) - 1)
            IF (N2 .LT. IFHOLE(J) + NHPR(J) - 1) THEN
               IF (N1 .EQ. IFHOLE(J)) THEN
                  WRITE (IUNIT, 10040) ' HOLE  ', IABS (IREGN(J)),
     &               (IHLIST(K), K = N1, N2)
               ELSE
                  WRITE (IUNIT, 10050) (IHLIST(K), K = N1, N2)
               END IF
               GO TO 280
            ELSE
               IF (N1 .EQ. IFHOLE(J)) THEN
                  WRITE (IUNIT, 10060) ' HOLE  ', IABS (IREGN(J)),
     &               (IHLIST(K), K = N1, N2)
               ELSE
                  WRITE (IUNIT, 10070) (IHLIST(K), K = N1, N2)
               END IF
            END IF
         END IF
  290 CONTINUE

C  WRITE OUT THE GROUPS IN ORDER

      DO 310 I = 1, N(22)
         CALL LTSORT (MR, LINKR, I, J, ADDLNK)
         IF ((J .GT. 0) .AND. (IRGFLG(J) .GE. 1)) THEN
            N2 = IFSIDE(J) - 1
  300       CONTINUE
            N1 = N2 + 1
            IF (N1 .EQ. IFSIDE(J)) THEN
               N2 = N1 + 9
            ELSE
               N2 = N1 + 10
            END IF
            N2 = MIN0(N2, IFSIDE(J) + NSPR(J) - 1)
            IF (N2 .LT. IFSIDE(J) + NSPR(J) - 1) THEN
               IF (N1 .EQ. IFSIDE(J)) THEN
                  WRITE (IUNIT, 10040) ' GROUP ', IABS (IREGN(J)),
     &               (ISLIST(K), K = N1, N2)
               ELSE
                  WRITE (IUNIT, 10050) (ISLIST(K), K = N1, N2)
               END IF
               GO TO 300
            ELSE
               IF (N1 .EQ. IFSIDE(J)) THEN
                  WRITE (IUNIT, 10060) 'GROUP ', IABS (IREGN(J)),
     &               (ISLIST(K), K = N1, N2)
               ELSE
                  WRITE (IUNIT, 10070) (ISLIST(K), K = N1, N2)
               END IF
            END IF
         END IF
  310 CONTINUE

C  WRITE OUT THE SCHEMES IN ORDER

      DO 320 I = 1, N(24)
         CALL LTSORT (MR, LINKSC, I, J, ADDLNK)
         IF (J .GT. 0) THEN
            CALL LTSORT (MR, LINKR, ISCHM (J), JJ, ADDLNK)
            IF ((JJ .GT. 0) .AND. (IREGN (JJ) .LT. 0)) THEN
               CALL STRLNG (SCHEME(J), LEN)
               DUMMY = SCHEME(J)(1:LEN)
               WRITE (IUNIT, 10110) ISCHM(J), DUMMY(1:LEN)
            END IF
         END IF
  320 CONTINUE
      CALL STRLNG (DEFSCH, LEN)
      IZERO = 0
      WRITE (IUNIT, 10110) IZERO, DEFSCH(1:LEN)

C  WRITE OUT THE BODY LIST

      N2 = 0
      KMAX = 11
  330 CONTINUE
      N1 = N2 + 1
      IF (N(9) .GT. 0) THEN
         KOUNT = 0
         DO 340 I = N1, N(9)
            IF (IRPB(I) .GT. 0) THEN
               CALL LTSORT (MR, LINKR, IRPB(I), J, ADDLNK)
            ELSE IF (IRPB(I) .LT. 0) THEN
               CALL LTSORT (MS, LINKB, IABS(IRPB(I)), J, ADDLNK)
            ELSE
               J = 0
            ENDIF
            IF ( ((IRPB(I) .GT. 0) .AND. (J .GT. 0) .AND.
     &         (IREGN (J) .LT. 0)) .OR.
     &         ((IRPB(I) .LT. 0) .AND. (J .GT. 0) .AND.
     &         (IBARST (J) .LT. 0)) ) THEN
               KOUNT = KOUNT + 1
               ID (KOUNT) = IRPB (I)
            ENDIF
            N2 = I
            IF (KOUNT .EQ. KMAX) GOTO 350
  340    CONTINUE
  350    CONTINUE
         IF (N2 .LT. N(9)) THEN
            IF (N1 .EQ. 1) THEN
               WRITE (IUNIT, 10040) ' BODY  ', (ID (I), I = 1, KOUNT)
            ELSE
               WRITE (IUNIT, 10050) (ID(I), I = 1, KOUNT)
            END IF
            GO TO 330
         ELSEIF (KOUNT .GT. 0) THEN
            IF (N1 .EQ. 1) THEN
               WRITE (IUNIT, 10060) 'BODY  ', (ID (I), I = 1, KOUNT)
            ELSE
               WRITE (IUNIT, 10070) (ID (I), I = 1, KOUNT)
            END IF
         ELSE
            WRITE (IUNIT, 10000)
         END IF
      END IF

C  WRITE OUT THE POINT BOUNDARY FLAGS IN ORDER

      DO 390 I = 1, N(25)
         CALL LTSORT (MP, LINKPB, I, J, ADDLNK)
         IF (J .GT. 0) THEN
            STAR = .FALSE.
            N2 = IFPB(J) - 1
  360       CONTINUE
            N1 = N2 + 1
            IF (N1 .EQ. IFPB(J)) THEN
               KMAX = 10
            ELSE
               KMAX = 11
            END IF

            KOUNT = 0
            DO 370 K = N1, IFPB(J) + NPPF(J) - 1
               CALL LTSORT (MP, LINKP, LISTPB (1, K), JJ, ADDLNK)
               IF ((JJ .GT. 0) .AND. (IPOINT (JJ) .LT. 0)) THEN
                  KOUNT = KOUNT + 1
                  ID (KOUNT) = LISTPB (1, K)
               ENDIF
               N2 = K
               IF (KOUNT .EQ. KMAX) GOTO 380
  370       CONTINUE
  380       CONTINUE

            IF (N2 .LT. IFPB(J) + NPPF(J) - 1) THEN
               STAR = .TRUE.
               IF (N1 .EQ. IFPB(J)) THEN
                  WRITE (IUNIT, 10040) ' POINBC', IPBF(J),
     &               (ID (K), K = 1, KOUNT)
               ELSE
                  WRITE (IUNIT, 10050) (ID (K), K = 1, KOUNT)
               END IF
               GO TO 360
            ELSEIF (KOUNT .GT. 0) THEN
               IF (N1 .EQ. IFPB(J)) THEN
                  WRITE (IUNIT, 10060) ' POINBC', IPBF(J),
     &               (ID (K), K = 1, KOUNT)
               ELSE
                  WRITE (IUNIT, 10070) (ID (K), K = 1, KOUNT)
               END IF
            ELSEIF (STAR) THEN
               WRITE (IUNIT, 10000)
            END IF
            IF (IWTPBF(1, J) .GT. 0) THEN
               WRITE (IUNIT, 10100) 'POINT', IPBF(J), IWTPBF(1, J),
     &            IWTPBF(2, J)
            END IF
         END IF
  390 CONTINUE

C  WRITE OUT THE LINE BOUNDARY FLAGS IN ORDER

      DO 430 I = 1, N(26)
         CALL LTSORT (ML, LINKLB, I, J, ADDLNK)
         IF (J .GT. 0) THEN
            STAR = .FALSE.
            N2 = IFLB(J) - 1
  400       CONTINUE
            N1 = N2 + 1
            IF (N1 .EQ. IFLB(J)) THEN
               KMAX = 10
            ELSE
               KMAX = 11
            END IF

            KOUNT = 0
            DO 410 K = N1, IFLB(J) + NLPF(J) - 1
               CALL LTSORT (ML, LINKL, LISTLB (1, K), JJ, ADDLNK)
               IF ((JJ .GT. 0) .AND. (ILINE (JJ) .LT. 0)) THEN
                  KOUNT = KOUNT + 1
                  ID (KOUNT) = LISTLB (1, K)
               ENDIF
               N2 = K
               IF (KOUNT .EQ. KMAX) GOTO 420
  410       CONTINUE
  420       CONTINUE

            IF (N2 .LT. IFLB(J) + NLPF(J) - 1) THEN
               STAR = .TRUE.
               IF (N1 .EQ. IFLB(J)) THEN
                  WRITE (IUNIT, 10040) ' NODEBC', ILBF(J),
     &               (ID (K), K = 1, KOUNT)
               ELSE
                  WRITE (IUNIT, 10050) (ID (K), K = 1, KOUNT)
               END IF
               GO TO 400
            ELSEIF (KOUNT .GT. 0) THEN
               IF (N1 .EQ. IFLB(J)) THEN
                  WRITE (IUNIT, 10060) ' NODEBC', ILBF(J),
     &               (ID (K), K = 1, KOUNT)
               ELSE
                  WRITE (IUNIT, 10070) (ID (K), K = 1, KOUNT)
               END IF
            ELSEIF (STAR) THEN
               WRITE (IUNIT, 10000)
            END IF
            IF (IWTLBF(1, J) .NE. 0) THEN
               IF (IWTLBF(3, J) .GT. 0) THEN
                  WRITE (IUNIT, 10090) '  LINE', ILBF(J), IWTLBF(1, J),
     &               IWTLBF(2, J), IWTLBF(3, J)
               ELSE
                  WRITE (IUNIT, 10100) '  LINE', ILBF(J), IWTLBF(1, J),
     &               IWTLBF(2, J)
               END IF
            END IF
         END IF
  430 CONTINUE

C  WRITE OUT THE SIDE BOUNDARY FLAGS IN ORDER

      DO 470 I = 1, N(27)
         CALL LTSORT (ML, LINKSB, I, J, ADDLNK)
         IF (J .GT. 0) THEN
            STAR = .FALSE.
            N2 = IFSB(J) - 1
  440       CONTINUE
            N1 = N2 + 1
            IF (N1 .EQ. IFSB(J)) THEN
               KMAX = 10
            ELSE
               KMAX = 11
            END IF

            KOUNT = 0
            DO 450 K = N1, IFSB(J) + NSPF(J) - 1
               CALL LTSORT (ML, LINKL, LISTSB (1, K), JJ, ADDLNK)
               IF ((JJ .GT. 0) .AND. (ILINE (JJ) .LT. 0)) THEN
                  KOUNT = KOUNT + 1
                  ID (KOUNT) = LISTSB (1, K)
               ENDIF
               N2 = K
               IF (KOUNT .EQ. KMAX) GOTO 460
  450       CONTINUE
  460       CONTINUE

            IF (N2 .LT. IFSB(J) + NSPF(J) - 1) THEN
               STAR = .TRUE.
               IF (N1 .EQ. IFSB(J)) THEN
                  WRITE (IUNIT, 10040) ' ELEMBC', ISBF(J),
     &               (ID (K), K = 1, KOUNT)
               ELSE
                  WRITE (IUNIT, 10050) (ID (K), K = 1, KOUNT)
               END IF
               GO TO 440
            ELSEIF (KOUNT .GT. 0) THEN
               IF (N1 .EQ. IFSB(J)) THEN
                  WRITE (IUNIT, 10060) ' ELEMBC', ISBF(J),
     &               (ID (K), K = 1, KOUNT)
               ELSE
                  WRITE (IUNIT, 10070) (ID (K), K = 1, KOUNT)
               END IF
            ELSEIF (STAR) THEN
               WRITE (IUNIT, 10000)
            END IF
            IF (IWTSBF(1, J) .NE. 0) THEN
               IF (IWTSBF(3, J) .GT. 0) THEN
                  WRITE (IUNIT, 10090) ' SIDE', ISBF(J), IWTSBF(1, J),
     &               IWTSBF(2, J), IWTSBF(3, J)
               ELSE
                  WRITE (IUNIT, 10100) ' SIDE', ISBF(J), IWTSBF(1, J),
     &               IWTSBF(2, J)
               END IF
            END IF
         END IF
  470 CONTINUE

C  WRITE OUT THE RENUMBERING CARDS

      IF (OPTIM) THEN
         IF (N(28) .GT. 0) THEN
            DO 480 I = 1, N(28)
               DUMMY = NUMBER(I)(6:77)
               CALL STRLNG (DUMMY, LEN)
               WRITE (IUNIT, 10120) NUMBER(I)(1:5), NUMBER(I)(6:6 + LEN)
  480       CONTINUE
         ELSE
            WRITE (IUNIT, 10130)
         END IF
      END IF

C  WRITE OUT THREE NODE, EIGHT NODE, OR NINE NODE FLAG

      IF (THREE) THEN
         WRITE (IUNIT, 10140)
      ENDIF
      IF (EIGHT) THEN
         WRITE (IUNIT, 10150)
      ELSE IF (NINE) THEN
         WRITE (IUNIT, 10160)
      END IF

C  WRITE DIGITIZER SNAP-TO-GRID FLAG

      IF ((NSNAP(1) .GT. 0) .OR. (NSNAP(2) .GT. 0)) THEN
         IF (SNAP) THEN
            WRITE (IUNIT, 10180)
         ELSE
            WRITE (IUNIT, 10190)
         END IF

C  WRITE X-GRID LINES

         IF (NSNAP(1) .GT. 0) THEN
            N2 = 0
  490       CONTINUE
            N1 = N2 + 1
            N2 = N1 + 4
            N2 = MIN0(N2, NSNAP(1))
            IF (N2 .LT. NSNAP(1)) THEN
               IF (N1 .EQ. 1) THEN
                  WRITE (IUNIT, 10200) 'XGRID ',
     &               (SNAPDX(1, K), K = N1, N2)
               ELSE
                  WRITE (IUNIT, 10210) (SNAPDX(1, K), K = N1, N2)
               END IF
               GO TO 490
            ELSE
               IF (N1 .EQ. 1) THEN
                  WRITE (IUNIT, 10220) ' XGRID ',
     &               (SNAPDX(1, K), K = N1, N2)
               ELSE
                  WRITE (IUNIT, 10230) (SNAPDX(1, K), K = N1, N2)
               END IF
            END IF

C  WRITE Y-GRID LINES

            IF (NSNAP(2) .GT. 0) THEN
               N2 = 0
  500          CONTINUE
               N1 = N2 + 1
               N2 = N1 + 4
               N2 = MIN0(N2, NSNAP(2))
               IF (N2 .LT. NSNAP(2)) THEN
                  IF (N1 .EQ. 1) THEN
                     WRITE (IUNIT, 10200) ' YGRID ',
     &                  (SNAPDX(2, K), K = N1, N2)
                  ELSE
                     WRITE (IUNIT, 10210) (SNAPDX(2, K), K = N1, N2)
                  END IF
                  GO TO 500
               ELSE
                  IF (N1 .EQ. 1) THEN
                     WRITE (IUNIT, 10220) 'YGRID ',
     &                  (SNAPDX(2, K), K = N1, N2)
                  ELSE
                     WRITE (IUNIT, 10230) (SNAPDX(2, K), K = N1, N2)
                  END IF
               END IF
            END IF
         END IF
      END IF

C  WRITE EXIT

      WRITE (IUNIT, 10170)

      CALL MESAGE ('FASTQ DATA FILE SUCCESSFULLY WRITTEN')
  510 CONTINUE
      FLAG=.FALSE.
      CALL FLAGD (MP, N(18), LINKP, IPOINT, FLAG)
      CALL FLAGD (ML, N(19), LINKL, ILINE, FLAG)
      CALL FLAGD (MS, N(20), LINKS, ISIDE, FLAG)
      CALL FLAGD (MS, N(21), LINKB, IBARST, FLAG)
      CALL FLAGD (MR, N(22), LINKR, IREGN, FLAG)
      RETURN

10000 FORMAT ('      ')
10010 FORMAT (' TITLE ', /, ' ', A)
10020 FORMAT (' POINT ', I5, 2(5X, 1PE14.7))
10030 FORMAT (' LINE  ', I5, 1X, A5, 4(1X, I5), 2X, F6.4)
10040 FORMAT (A7, 11(1X, I5), ' *')
10050 FORMAT (7X, 11(1X, I5), ' *')
10060 FORMAT (A7, 11(1X, I5))
10070 FORMAT (7X, 11(1X, I5))
10080 FORMAT (A7, 1X, 1PE14.7, 2X, I5)
10090 FORMAT (' WEIGHT', 1X, A5, 4(1X, I5))
10100 FORMAT (' WEIGHT', 1X, A5, 3(1X, I5))
10110 FORMAT (' SCHEME', I5, 1X, A)
10120 FORMAT (' RENUM ', A5, 1X, A)
10130 FORMAT (' RENUM ')
10140 FORMAT (' THREE ')
10150 FORMAT (' EIGHT ')
10160 FORMAT (' NINE  ')
10170 FORMAT (' EXIT  ')
10180 FORMAT (' SNAP ON')
10190 FORMAT (' SNAP OFF')
10200 FORMAT (A7, 5(1X, 1PE13.6), ' *')
10210 FORMAT (7X, 5(1X, 1PE13.6), ' *')
10220 FORMAT (A7, 5(1X, 1PE13.6))
10230 FORMAT (7X, 5(1X, 1PE13.6))

      END
