C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE DELFSQ (MP, ML, MS, MR, MSC, MCOM, ICOM, JCOM, CIN,
     &   RIN, IIN, KIN, N, IPBOUN, ILBOUN, ISBOUN, NLPS, IFLINE, ILLIST,
     &   NSPR, IFSIDE, ISLIST, IRPB, IPBF, NPPF, IFPB, LISTPB, ILBF,
     &   NLPF, IFLB, LISTLB, ISBF, NSPF, IFSB, LISTSB, LINKP, LINKL,
     &   LINKS, LINKB, LINKR, LINKSC, LINKPB, LINKLB, LINKSB, IWTPBF,
     &   IWTLBF, IWTSBF, IFHOLE, NHPR, IHLIST, IRGFLG, NUMBER, DEFSCH,
     &   OPTIM, VAXVMS, WROTE, TIME1, BATCH, VERSN)
C***********************************************************************

C  SUBROUTINE DELFSQ = DELETES POINTS, LINES, REGIONS, SCHEMES, AND
C                      BOUNDARY DEFINITIONS

C***********************************************************************

C  SUBROUTINE CALLED BY:
C     FASTQ = A PROGRAM TO QUICKLY PREPARE QMESH INPUT

C***********************************************************************

C  VARIABLES USED:
C     IANS   = LOGICAL RESPONSE FROM YES-NO QUESTION
C     ANS    = CHARACTER RESPONSE FOR MENU CHOICE

C***********************************************************************

      DIMENSION IPBOUN(MP), ILBOUN(ML), ISBOUN(ML)
      DIMENSION NLPS(MS), IFLINE(MS), ILLIST(MS*3)
      DIMENSION NSPR(MR), IFSIDE(MR), ISLIST(4*MR), IRPB(MR)
      DIMENSION IPBF(MP), NPPF(MP), IFPB(MP), LISTPB(2, MP)
      DIMENSION ILBF(ML), NLPF(ML), IFLB(ML), LISTLB(2, ML)
      DIMENSION ISBF(ML), NSPF(ML), IFSB(ML), LISTSB(2, ML)
      DIMENSION IWTPBF(3, MP), IWTLBF(3, ML), IWTSBF(3, ML)
      DIMENSION LINKP(2, MP), LINKL(2, ML), LINKS(2, MS), LINKB(2, MS)
      DIMENSION LINKR(2, MR), LINKSC(2, MR), LINKPB(2, MP)
      DIMENSION LINKLB(2, ML), LINKSB(2, ML), NUMBER(MSC)
      DIMENSION IFHOLE(MR), NHPR(MR), IHLIST(MR*2), IRGFLG(MR)
      DIMENSION N(29)
      DIMENSION KIN(MCOM), CIN(MCOM), IIN(MCOM), RIN(MCOM)

      CHARACTER*72 DEFSCH, CIN
      CHARACTER VERSN*9, NUMBER*80

      LOGICAL OPTIM, ADDLNK, VAXVMS, WROTE, SIDEOK, BATCH, NOROOM
      LOGICAL LGROUP

      ADDLNK = .TRUE.
      IZ = 0

  100 CONTINUE
      IF (ICOM .GT. JCOM) THEN
         CALL MESAGE (' ')
         CALL FREFLD (IZ, IZ, 'ENTER DELETE OPTION: ', MCOM, IOSTAT,
     &      JCOM, KIN, CIN, IIN, RIN)
         ICOM = 1
      END IF

      IF ((CIN(ICOM)(1:1) .EQ. 'P') .OR. (CIN(ICOM)(1:1) .EQ. 'p')) THEN
         ICOM = ICOM + 1
         IF (N(1) .GT. 0) THEN
            CALL MESAGE (' ')
            CALL MESAGE ('DELETE POINTS <I1> THROUGH <I2>:')
            CALL MESAGE ('HIT A RETURN TO END')
  110       CONTINUE
            IF (ICOM .GT. JCOM) THEN
               CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN, CIN,
     &            IIN, RIN)
               ICOM = 1
            END IF
            CALL GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2,
     &         IFOUND)
            IF (IFOUND .GT. 0) THEN
               CALL CHECK (I1, I2, N(18))
               DO 120 I = I1, I2
                  CALL LTSORT (MP, LINKP, I, IZ, ADDLNK)
  120          CONTINUE
               GO TO 110
            END IF
         ELSE
            CALL MESAGE (' ')
            CALL MESAGE ('*-----------------------------------*')
            CALL MESAGE ('* NO POINTS IN THE CURRENT DATABASE *')
            CALL MESAGE ('*-----------------------------------*')
         END IF

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'L') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'l')) THEN
         ICOM = ICOM + 1
         IF (N(2) .GT. 0) THEN
            CALL MESAGE (' ')
            CALL MESAGE ('DELETE LINES <I1> THROUGH <I2>:')
            CALL MESAGE ('HIT A RETURN TO END')
  130       CONTINUE
            IF (ICOM .GT. JCOM) THEN
               CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN, CIN,
     &            IIN, RIN)
               ICOM = 1
            END IF
            CALL GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2,
     &         IFOUND)
            IF (IFOUND .GT. 0) THEN
               CALL CHECK (I1, I2, N(19))
               DO 140 I = I1, I2
                  CALL LTSORT (ML, LINKL, I, IZ, ADDLNK)
  140          CONTINUE
               GO TO 130
            END IF
         ELSE
            CALL MESAGE (' ')
            CALL MESAGE ('*----------------------------------*')
            CALL MESAGE ('* NO LINES IN THE CURRENT DATABASE *')
            CALL MESAGE ('*----------------------------------*')
         END IF

C  DELETE BAR SET DEFINITIONS

      ELSE IF ((CIN(ICOM)(1:2) .EQ. 'BA') .OR.
     &   (CIN(ICOM)(1:2) .EQ. 'ba')) THEN
         ICOM = ICOM + 1
         IF (N(5) .GT. 0) THEN
            CALL MESAGE (' ')
            CALL MESAGE ('DELETE BAR SETS <I1> THROUGH <I2>:')
            CALL MESAGE ('HIT A RETURN TO END')
  150       CONTINUE
            IF (ICOM .GT. JCOM) THEN
               CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN, CIN,
     &            IIN, RIN)
               ICOM = 1
            END IF
            CALL GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2,
     &         IFOUND)
            IF (IFOUND .GT. 0) THEN
               CALL CHECK (I1, I2, N(21))
               DO 160 I = I1, I2
                  CALL LTSORT (MS, LINKB, I, IZ, ADDLNK)
  160          CONTINUE
               GO TO 150
            END IF
         ELSE
            CALL MESAGE (' ')
            CALL MESAGE ('*-------------------------------------*')
            CALL MESAGE ('* NO BAR SETS IN THE CURRENT DATABASE *')
            CALL MESAGE ('*-------------------------------------*')
         END IF

C  DELETE THE RENUMBERING CARDS

      ELSE IF ((CIN(ICOM)(1:3) .EQ. 'REN') .OR.
     &   (CIN(ICOM)(1:3) .EQ. 'ren')) THEN
         ICOM = ICOM + 1
         IF (N(28) .GT. 0) THEN
            CALL MESAGE (' ')
            I1 = 1
            I2 = N(28)
            ITEST = 20
            DO 170 I = I1, I2
               WRITE(*, 10040) I, NUMBER(I)(1:72)
               IF (NUMBER(I)(73:80) .NE.'        ') THEN
                  WRITE(*, 10050) NUMBER(I)(73:80)
               END IF
  170       CONTINUE
            CALL MESAGE (' ')
            CALL MESAGE ('DELETE RENUMBER CARDS <I1> THROUGH <I2>:')
            CALL MESAGE ('HIT A RETURN TO END')
  180       CONTINUE
            IF (ICOM .GT. JCOM) THEN
               CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN, CIN,
     &            IIN, RIN)
               ICOM = 1
            END IF
            CALL GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2,
     &         IFOUND)
            IF (IFOUND .GT. 0) THEN
               CALL CHECK (I1, I2, N(28))
               DO 200 I = I1, I2
                  DO 190 J = I, N(28)
                     NUMBER(J) = NUMBER(J + 1)
  190             CONTINUE
                  N(28) = N(28) - 1
  200          CONTINUE
               IF (N(28) .EQ. 0) OPTIM = .FALSE.
               GO TO 180
            END IF
         ELSE
            CALL MESAGE (' ')
            CALL MESAGE
     &         ('*-------------------------------------------*')
            IF (.NOT. OPTIM) CALL MESAGE
     &         ('*        NO RENUMBER CONTROL CARDS          *')
            CALL MESAGE
     &         ('*        OPTIMIZATION HAS BEEN DISABLED     *')
            CALL MESAGE
     &         ('*-------------------------------------------*')
            CALL MESAGE (' ')
            OPTIM = .FALSE.
         END IF

C  DELETE THE REGIONS

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'R') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'r')) THEN
         ICOM = ICOM + 1
         LGROUP = .TRUE.
         DO 210 J = 1, N(7)
            IF (IRGFLG(J) .LE. -1) THEN
               LGROUP = .FALSE.
               GO TO 220
            END IF
  210    CONTINUE
  220    CONTINUE
         IF (.NOT.LGROUP) THEN
            CALL MESAGE (' ')
            CALL MESAGE ('DELETE REGIONS <I1> THROUGH <I2>:')
            CALL MESAGE ('HIT A RETURN TO END')
  230       CONTINUE
            IF (ICOM .GT. JCOM) THEN
               CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN, CIN,
     &            IIN, RIN)
               ICOM = 1
            END IF
            CALL GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2,
     &         IFOUND)
            IF (IFOUND .GT. 0) THEN
               CALL CHECK (I1, I2, N(22))
               DO 320 I = I1, I2
                  ADDLNK = .FALSE.
                  CALL LTSORT (MR, LINKR, I, IPTR, ADDLNK)
                  IF ((IPTR .GT. 0) .AND. (IRGFLG(IPTR) .LE. 0)) THEN

C  DELETE REGION FROM BODY LIST

                     DO 250 J = 1, N(9)
                        IF (IRPB(J) .EQ. I) THEN
                           DO 240 K = J + 1, N(9)
                              IRPB(K - 1) = IRPB(K)
  240                      CONTINUE
                           N(9) = N(9) - 1
                        END IF
  250                CONTINUE

C  DELETE REGION FROM GROUPS

                     DO 280 J = 1, N(7)
                        IF (IRGFLG(J) .GE. 1) THEN
                           K1 = IFSIDE(IPTR)
                           K2 = K1 + NSPR(IPTR) - 1
                           DO 270 K = K1, K2
                              IF (ISLIST(K) .EQ. I) THEN
                                 DO 260 L = K + 1, K2
                                    ISLIST(K - 1) = ISLIST(K)
  260                            CONTINUE
                                 NSPR(IPTR) = NSPR(IPTR) - 1
                              END IF
  270                      CONTINUE
                        END IF
  280                CONTINUE

C  DELETE REGION FROM HOLES

                     DO 310 J = 1, N(7)
                        IF (NHPR(J) .GE. 1) THEN
                           K1 = IFHOLE(IPTR)
                           K2 = K1 + NHPR(IPTR) - 1
                           DO 300 K = K1, K2
                              IF (IHLIST(K) .EQ. I) THEN
                                 DO 290 L = K + 1, K2
                                    IHLIST(K - 1) = IHLIST(K)
  290                            CONTINUE
                                 NHPR(IPTR) = NHPR(IPTR) - 1
                              END IF
  300                      CONTINUE
                        END IF
  310                CONTINUE

C  DELETE LINK TO REGION

                     ADDLNK = .TRUE.
                     CALL LTSORT (MR, LINKR, I, IZ, ADDLNK)
                  END IF
  320          CONTINUE
               ADDLNK = .TRUE.
               GO TO 230
            END IF
         ELSE
            CALL MESAGE (' ')
            CALL MESAGE ('*------------------------------------*')
            CALL MESAGE ('* NO REGIONS IN THE CURRENT DATABASE *')
            CALL MESAGE ('*------------------------------------*')
         END IF

C  DELETE THE GROUPS

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'G') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'g')) THEN
         ICOM = ICOM + 1
         LGROUP = .FALSE.
         DO 330 J = 1, N(7)
            IF (IRGFLG(J) .GE. 1) THEN
               LGROUP = .TRUE.
               GO TO 340
            END IF
  330    CONTINUE
  340    CONTINUE
         IF (LGROUP) THEN
            CALL MESAGE (' ')
            CALL MESAGE ('DELETE GROUPS<I1> THROUGH <I2>:')
            CALL MESAGE ('HIT A RETURN TO END')
  350       CONTINUE
            IF (ICOM .GT. JCOM) THEN
               CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN, CIN,
     &            IIN, RIN)
               ICOM = 1
            END IF
            CALL GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2,
     &         IFOUND)
            IF (IFOUND .GT. 0) THEN
               CALL CHECK (I1, I2, N(22))
               DO 380 I = I1, I2
                  ADDLNK = .FALSE.
                  CALL LTSORT (MR, LINKR, I, IPTR, ADDLNK)
                  IF ((IPTR .GT. 0) .AND. (IRGFLG(IPTR) .GE. 1)) THEN

C  DELETE GROUP FROM BODY LIST

                     DO 370 J = 1, N(9)
                        IF (IRPB(J) .EQ. I) THEN
                           DO 360 K = J + 1, N(9)
                              IRPB(K - 1) = IRPB(K)
  360                      CONTINUE
                           N(9) = N(9) - 1
                        END IF
  370                CONTINUE

C  DELETE LINK TO GROUP

                     ADDLNK = .TRUE.
                     CALL LTSORT (MR, LINKR, I, IZ, ADDLNK)
                  END IF
  380          CONTINUE
               ADDLNK = .TRUE.
               GO TO 350
            END IF
         ELSE
            CALL MESAGE (' ')
            CALL MESAGE ('*-----------------------------------*')
            CALL MESAGE ('* NO GROUPS IN THE CURRENT DATABASE *')
            CALL MESAGE ('*-----------------------------------*')
         END IF

C  DELETE THE HOLES

      ELSE IF ((CIN(ICOM)(1:2) .EQ. 'HO') .OR.
     &   (CIN(ICOM)(1:2) .EQ. 'ho')) THEN
         ICOM = ICOM + 1
         IF (N(29) .GT. 0) THEN
            CALL MESAGE (' ')
            CALL MESAGE ('DELETE HOLES FOR REGIONS <I1> THROUGH <I2>:')
            CALL MESAGE ('HIT A RETURN TO END')
  390       CONTINUE
            IF (ICOM .GT. JCOM) THEN
               CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN, CIN,
     &            IIN, RIN)
               ICOM = 1
            END IF
            CALL GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2,
     &         IFOUND)
            IF (IFOUND .GT. 0) THEN
               CALL CHECK (I1, I2, N(22))
               ADDLNK = .FALSE.
               DO 400 I = I1, I2
                  CALL LTSORT (MR, LINKR, I, L, ADDLNK)
                  IF (L .GT. 0) THEN
                     IF (NHPR(L) .GT. 0) THEN
                        NHPR(L) = 0
                     ELSE
                        WRITE (*, 10000) I
                     END IF
                  ELSE
                     WRITE (*, 10010) I
                  END IF
  400          CONTINUE
               ADDLNK = .TRUE.
               GO TO 390
            END IF
         ELSE
            CALL MESAGE (' ')
            CALL MESAGE ('*----------------------------------*')
            CALL MESAGE ('* NO HOLES IN THE CURRENT DATABASE *')
            CALL MESAGE ('*----------------------------------*')
         END IF

C  DELETE SCHEMES

      ELSE IF ((CIN(ICOM)(1:2) .EQ. 'SC') .OR.
     &   (CIN(ICOM)(1:2) .EQ. 'sc')) THEN
         ICOM = ICOM + 1
         IF (N(10) .GT. 0) THEN
            CALL MESAGE (' ')
            CALL MESAGE ('DELETE SCHEMES <I1> THROUGH <I2>:')
            CALL MESAGE ('HIT A RETURN TO END')
  410       CONTINUE
            IF (ICOM .GT. JCOM) THEN
               CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN, CIN,
     &            IIN, RIN)
               ICOM = 1
            END IF
            CALL GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2,
     &         IFOUND)
            IF (IFOUND .GT. 0) THEN
               CALL CHECK (I1, I2, N(24))
               DO 420 I = I1, I2
                  CALL LTSORT (MR, LINKSC, I, IZ, ADDLNK)
  420          CONTINUE
               GO TO 410
            END IF
         ELSE
            CALL MESAGE (' ')
            CALL MESAGE
     &         ('*---------------------------------------------*')
            CALL MESAGE
     &         ('* ONLY DEFAULT SCHEME IN THE CURRENT DATABASE *')
            CALL MESAGE
     &         ('*---------------------------------------------*')
            WRITE(*, 10020) DEFSCH
            CALL MESAGE (' ')
         END IF

C  SPAWN A PROCESS

      ELSE IF ((CIN(ICOM)(1:2) .EQ. 'SP') .OR.
     &   (CIN(ICOM)(1:2) .EQ. 'sp')) THEN
         ICOM = ICOM + 1
         CALL SPAWN (VAXVMS)

C  DELETE SIDE DEFINITIONS

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'S') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 's')) THEN
         ICOM = ICOM + 1
         IF (N(3) .GT. 0) THEN
            CALL MESAGE (' ')
            CALL MESAGE ('DELETE SIDES <I1> THROUGH <I2>:')
            CALL MESAGE ('HIT A RETURN TO END')
  430       CONTINUE
            IF (ICOM .GT. JCOM) THEN
               CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN, CIN,
     &            IIN, RIN)
               ICOM = 1
            END IF
            CALL GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2,
     &         IFOUND)
            IF (IFOUND .GT. 0) THEN
               CALL CHECK (I1, I2, N(20))
               DO 440 I = I1, I2
                  CALL LTSORT (MS, LINKS, I, IZ, ADDLNK)
  440          CONTINUE
               GO TO 430
            END IF
         ELSE
            CALL MESAGE (' ')
            CALL MESAGE ('*----------------------------------*')
            CALL MESAGE ('* NO SIDES IN THE CURRENT DATABASE *')
            CALL MESAGE ('*----------------------------------*')
         END IF

C  DELETE BOUNDARY CONDITIONS

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'B') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'b')) THEN
         ICOM = ICOM + 1
         IF (ICOM .GT. JCOM) THEN
            CALL MESAGE
     &         ('THE FOLLOWING BOUNDARY FLAG TYPES ARE AVAILABLE')
            CALL MESAGE
     &         ('        P*OINT FLAGS    - FOR NODES AT POINTS')
            CALL MESAGE
     &         ('        N*ODE FLAGS     - FOR NODES ON A BOUNDARY')
            CALL MESAGE ('        E*LEMENT FLAGS  - FOR ELEMENT SIDES '
     &         //'ON A BOUNDARY')
            CALL FREFLD (IZ, IZ, 'TYPE OF BOUNDARY FLAG TO BE '//
     &         'DELETED FROM: ', MCOM, IOSTAT, JCOM, KIN, CIN, IIN, RIN)
            ICOM = 1
         END IF
         IF ((CIN(ICOM)(1:1) .EQ. 'P') .OR.
     &      (CIN(ICOM)(1:1) .EQ. 'p')) THEN
            ICOM = ICOM + 1
            IF (N(11) .GT. 0) THEN
               CALL MESAGE (' ')
               CALL MESAGE ('DELETE POINT BOUNDARY FLAGS <I1> '//
     &            'THROUGH <I2>:')
               CALL MESAGE ('HIT A RETURN TO END')
  450          CONTINUE
               IF (ICOM .GT. JCOM) THEN
                  CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN,
     &               CIN, IIN, RIN)
                  ICOM = 1
               END IF
               CALL GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2,
     1            IFOUND)
               IF (IFOUND .GT. 0) THEN
                  CALL CHECK (I1, I2, N(25))
                  DO 460 I = I1, I2
                     CALL LTSORT (MP, LINKPB, I, IZ, ADDLNK)
  460             CONTINUE
                  SIDEOK = .FALSE.
                  CALL LINKBC (MP, MS, 1, N(11), N(1), N(25), N(11),
     &               N(12), N(20), IPBF, IFPB, NPPF, LISTPB, NLPS,
     &               IFLINE, ILLIST, IPBOUN, LINKPB, IWTPBF, LINKP,
     &               LINKS, SIDEOK, NOROOM)
                  GO TO 450
               END IF
            ELSE
               CALL MESAGE (' ')
               CALL MESAGE
     &            ('*-----------------------------------------*')
               CALL MESAGE
     &            ('* NO POINT BOUNDARY FLAGS IN THE DATABASE *')
               CALL MESAGE
     &            ('*-----------------------------------------*')
               CALL MESAGE (' ')
            END IF
         ELSE IF ((CIN(ICOM)(1:1) .EQ. 'N') .OR.
     &      (CIN(ICOM)(1:1) .EQ. 'n')) THEN
            ICOM = ICOM + 1
            IF (N(13) .GT. 0) THEN
               CALL MESAGE (' ')
               CALL MESAGE
     &            ('DELETE NODE BOUNDARY FLAGS <I1> THROUGH <I2>:')
               CALL MESAGE ('HIT A RETURN TO END')
  470          CONTINUE
               IF (ICOM .GT. JCOM) THEN
                  CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN,
     &               CIN, IIN, RIN)
                  ICOM = 1
               END IF
               CALL GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2,
     &            IFOUND)
               IF (IFOUND .GT. 0) THEN
                  CALL CHECK (I1, I2, N(26))
                  DO 480 I = I1, I2
                     CALL LTSORT (ML, LINKLB, I, IZ, ADDLNK)
  480             CONTINUE

C  RELINK UP THE LINES TO THEIR ASSOCIATED FLAGS

                  SIDEOK = .TRUE.
                  CALL LINKBC (ML, MS, 1, N(13), N(2), N(26), N(13),
     &               N(14), N(20), ILBF, IFLB, NLPF, LISTLB, NLPS,
     &               IFLINE, ILLIST, ILBOUN, LINKLB, IWTLBF, LINKL,
     &               LINKS, SIDEOK, NOROOM)
                  GO TO 470
               END IF
            ELSE
               CALL MESAGE (' ')
               CALL MESAGE
     &            ('*----------------------------------------*')
               CALL MESAGE
     &            ('* NO NODE BOUNDARY FLAGS IN THE DATABASE *')
               CALL MESAGE
     &            ('*----------------------------------------*')
               CALL MESAGE (' ')
            END IF
         ELSE IF ((CIN(ICOM)(1:1) .EQ. 'E') .OR.
     &      (CIN(ICOM)(1:1) .EQ. 'e')) THEN
            ICOM = ICOM + 1
            IF (N(15) .GT. 0) THEN
               CALL MESAGE (' ')
               CALL MESAGE
     &            ('DELETE ELEMENT BOUNDARY FLAGS <I1> THROUGH <I2>:')
               CALL MESAGE ('HIT A RETURN TO END')
  490          CONTINUE
               IF (ICOM .GT. JCOM) THEN
                  CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN,
     &               CIN, IIN, RIN)
                  ICOM = 1
               END IF
               CALL GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2,
     &            IFOUND)
               IF (IFOUND .GT. 0) THEN
                  CALL CHECK (I1, I2, N(27))
                  DO 500 I = I1, I2
                     CALL LTSORT (ML, LINKSB, I, IZ, ADDLNK)
  500             CONTINUE

C  RELINK UP THE LINES TO THEIR ASSOCIATED FLAGS

                  SIDEOK = .TRUE.
                  CALL LINKBC (ML, MS, 1, N(15), N(2), N(27), N(15),
     &               N(16), N(20), ISBF, IFSB, NSPF, LISTSB, NLPS,
     &               IFLINE, ILLIST, ISBOUN, LINKSB, IWTSBF, LINKL,
     &               LINKS, SIDEOK, NOROOM)
                  GO TO 490
               END IF
            ELSE
               CALL MESAGE (' ')
               CALL MESAGE
     &            ('*-------------------------------------------*')
               CALL MESAGE
     &            ('* NO ELEMENT BOUNDARY FLAGS IN THE DATABASE *')
               CALL MESAGE
     &            ('*-------------------------------------------*')
               CALL MESAGE (' ')
            END IF
         END IF

C  EXIT OPTION - EXITS FASTQ

      ELSE IF ((CIN(ICOM)(1:2) .EQ. 'EX') .OR.
     &   (CIN(ICOM)(1:2) .EQ. 'ex')) THEN
         ICOM = ICOM + 1
         CALL STRLNG (CIN(ICOM), LEN)
         IF (((LEN .GT. 1) .AND. (CIN(ICOM)(2:2) .NE. 'X')) .OR.
     &      ((LEN .GT. 1) .AND. (CIN(ICOM)(2:2) .NE. 'x'))) THEN
            CALL HELP_FQ(6)
         ELSE
            CALL FEXIT (WROTE, MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN,
     &         TIME1, BATCH, VERSN)
         ENDIF
         GO TO 100
      ELSE IF (CIN(ICOM)(1:1) .EQ. ' ') THEN
         ICOM = ICOM + 1
         RETURN

C  PRINT HELP MESAGE

      ELSE
         ICOM = ICOM + 1
         CALL HELP_FQ(6)

      END IF
      GO TO 100

10000 FORMAT (' NO HOLES DEFINED IN REGION:', I6)
10010 FORMAT (' UNDEFINED REGION:', I6)
10020 FORMAT (' DEFLT: ', A72)
10040 FORMAT (1X, I5, 2X, A72)
10050 FORMAT (8X, A8)

      END
