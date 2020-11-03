C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE LIST (MP, ML, MS, MR, MSC, MCOM, ICOM, JCOM, CIN, RIN,
     &   IIN, KIN, N, IPOINT, COOR, IPBOUN, ILINE, LTYPE, NINT, FACTOR,
     &   LCON, ILBOUN, ISBOUN, ISIDE, NLPS, IFLINE, ILLIST, IBARST,
     &   JMAT, JCENT, NLPB, JFLINE, JLLIST, IREGN, IMAT, NSPR, IFSIDE,
     &   ISLIST, IRPB, IPBF, NPPF, IFPB, LISTPB, ILBF, NLPF, IFLB,
     &   LISTLB, ISBF, NSPF, IFSB, LISTSB, LINKP, LINKL, LINKS, LINKB,
     &   LINKR, LINKSC, LINKPB, LINKLB, LINKSB, IWTPBF, IWTLBF, IWTSBF,
     &   RSIZE, IFHOLE, NHPR, IHLIST, IRGFLG, ISCHM, SCHEME, NUMBER,
     &   DEFSCH, DEFSIZ, TITLE, OPTIM, THREE, EIGHT, NINE, VAXVMS,
     &   WROTE, TIME1, VERSN, BATCH)
C***********************************************************************

C  SUBROUTINE LIST = LISTS POINTS, LINES, REGIONS, SCHEMES, AND BOUNDARY
C                    DEFINITIONS

C***********************************************************************

C  SUBROUTINE CALLED BY:
C     FASTQ = A PROGRAM TO QUICKLY PREPARE QMESH INPUT

C***********************************************************************

C  SUBROUTINES CALLED:
C     CHECK  = CHECKS 2 VALUES FOR BEING OUT OF PRESCRIBED BOUNDS

C***********************************************************************

C  VARIABLES USED:
C     IANS   = LOGICAL RESPONSE FROM YES-NO QUESTION

C***********************************************************************

      DIMENSION IPOINT(MP), COOR(2, MP), IPBOUN(MP)
      DIMENSION ILINE(ML), LTYPE(ML), NINT(ML), FACTOR(ML), LCON(3, ML)
      DIMENSION ILBOUN(ML), ISBOUN(ML)
      DIMENSION ISIDE(MS), NLPS(MS), IFLINE(MS), ILLIST(MS*3)
      DIMENSION IBARST(MS), JMAT(MS), JCENT(MS), NLPB(MS), JFLINE(MS)
      DIMENSION JLLIST(MS*3)
      DIMENSION IREGN(MR), IMAT(MR), NSPR(MR), IFSIDE(MR), ISLIST(MR*4)
      DIMENSION IRPB(MR), ISCHM(MSC), SCHEME(MSC), RSIZE(MR)
      DIMENSION IPBF(MP), NPPF(MP), IFPB(MP), LISTPB(2, MP)
      DIMENSION ILBF(ML), NLPF(ML), IFLB(ML), LISTLB(2, ML)
      DIMENSION ISBF(ML), NSPF(ML), IFSB(ML), LISTSB(2, ML)
      DIMENSION IWTPBF(3, MP), IWTLBF(3, ML), IWTSBF(3, ML)
      DIMENSION LINKP(2, MP), LINKL(2, ML), LINKS(2, MS)
      DIMENSION LINKB(2, MS), LINKR(2, MR), LINKSC(2, MR)
      DIMENSION LINKPB(2, MP), LINKLB(2, ML), LINKSB(2, ML)
      DIMENSION IFHOLE(MR), NHPR(MR), IHLIST(MR*2), IRGFLG(MR)
      DIMENSION NUMBER(MSC), N(29)
      DIMENSION KIN(MCOM), CIN(MCOM), IIN(MCOM), RIN(MCOM)

      CHARACTER*72 SCHEME, DEFSCH, CIN, DUMMY*10, VERSN*9
      CHARACTER*72 TITLE, NUMBER*80, CHOICE*7

      LOGICAL IANS, OPTIM, ADDLNK, EIGHT, NINE, VAXVMS, WROTE, BATCH
      LOGICAL LGROUP, THREE

      IZ = 0
      ADDLNK = .FALSE.
      BATCH = .FALSE.

  100 CONTINUE
      IF (ICOM .GT. JCOM) THEN
         CALL MESAGE (' ')
         CALL FREFLD (IZ, IZ, 'ENTER LIST OPTION: ', MCOM, IOSTAT, JCOM,
     &      KIN, CIN, IIN, RIN)
         ICOM = 1
      END IF

C  LIST OUT THE POINTS

      IF ((CIN(ICOM)(1:1) .EQ. 'P') .OR. (CIN(ICOM)(1:1) .EQ. 'p')) THEN
         ICOM = ICOM+1
         IF (N(1) .GT. 0) THEN
            CALL MESAGE (' ')
            CALL INTRUP ('LIST ALL POINTS', IANS, MCOM, ICOM, JCOM, CIN,
     &         IIN, RIN, KIN)
            IF (IANS) THEN
               I1 = 1
               I2 = N(18)
            ELSE
               CALL MESAGE ('LIST POINTS <I1> THROUGH <I2>:')
               IF (ICOM .GT. JCOM) THEN
                  CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN,
     &               CIN, IIN, RIN)
                  ICOM = 1
               END IF
               CALL GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2,
     &            IFOUND)
               IF (IFOUND .GT. 0) THEN
                  CALL CHECK(I1, I2, N(18))
               ELSE
                  GO TO 120
               END IF
            END IF
            WRITE(*, 10000)
            DO 110 I = I1, I2
               CALL LTSORT (MP, LINKP, I, K, ADDLNK)
               IF (K .GT. 0) THEN
                  WRITE(*, 10010)IPOINT(K), COOR(1, K), COOR(2, K),
     &               IPBOUN(K)
               END IF
  110       CONTINUE
         ELSE
            CALL MESAGE (' ')
            CALL MESAGE ('*-----------------------------------*')
            CALL MESAGE ('* NO POINTS IN THE CURRENT DATABASE *')
            CALL MESAGE ('*-----------------------------------*')
         END IF
  120    CONTINUE

C  LIST OUT THE LINES

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'L') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'l')) THEN
         ICOM = ICOM+1
         IF (N(2) .GT. 0) THEN
            CALL MESAGE (' ')
            CALL INTRUP ('LIST ALL LINES', IANS, MCOM, ICOM, JCOM, CIN,
     &         IIN, RIN, KIN)
            IF (IANS) THEN
               I1 = 1
               I2 = N(19)
            ELSE
               CALL MESAGE ('LIST LINES <I1> THROUGH <I2>:')
               IF (ICOM .GT. JCOM) THEN
                  CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN,
     &               CIN, IIN, RIN)
                  ICOM = 1
               END IF
               CALL GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2,
     &            IFOUND)
               IF (IFOUND .GT. 0) THEN
                  CALL CHECK(I1, I2, N(19))
               ELSE
                  GO TO 140
               END IF
            END IF
            WRITE(*, 10020)
            DO 130 I = I1, I2
               CALL LTSORT (ML, LINKL, I, K, ADDLNK)
               IF (K .GT. 0) THEN
                  IF (LTYPE(K) .EQ. 1) THEN
                     WRITE(*, 10040)ILINE(K), (LCON(L, K), L = 1, 2),
     &                  NINT(K), ILBOUN(K), ISBOUN(K), FACTOR(K)
                  ELSE
                     WRITE(*, 10030)ILINE(K), (LCON(L, K), L = 1, 3),
     &                  NINT(K), ILBOUN(K), ISBOUN(K), FACTOR(K)
                  END IF
               END IF
  130       CONTINUE
         ELSE
            CALL MESAGE (' ')
            CALL MESAGE ('*----------------------------------*')
            CALL MESAGE ('* NO LINES IN THE CURRENT DATABASE *')
            CALL MESAGE ('*----------------------------------*')
         END IF
  140    CONTINUE

C  LIST OUT THE SIDES

      ELSE IF ((CIN(ICOM)(1:2) .EQ. 'SI') .OR.
     &   (CIN(ICOM)(1:2) .EQ. 'si')) THEN
         ICOM = ICOM+1
         IF (N(3) .GT. 0) THEN
            CALL MESAGE (' ')
            CALL INTRUP ('LIST ALL SIDES', IANS, MCOM, ICOM, JCOM, CIN,
     &         IIN, RIN, KIN)
            IF (IANS) THEN
               I1 = 1
               I2 = N(20)
            ELSE
               CALL MESAGE ('LIST SIDES <I1> THROUGH <I2>:')
               IF (ICOM .GT. JCOM) THEN
                  CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN,
     &               CIN, IIN, RIN)
                  ICOM = 1
               END IF
               CALL GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2,
     &            IFOUND)
               IF (IFOUND .GT. 0) THEN
                  CALL CHECK(I1, I2, N(20))
               ELSE
                  GO TO 170
               END IF
            END IF
            WRITE(*, 10090)
            DO 160 I = I1, I2
               CALL LTSORT (MS, LINKS, I, K, ADDLNK)
               IF (K .GT. 0) THEN
                  K1 = IFLINE(K)
  150             CONTINUE
                  K2 = MIN0(K1+10, IFLINE(K)+NLPS(K)-1)
                  IF (K1 .EQ. IFLINE(K)) THEN
                     WRITE(*, 10110)ISIDE(K), (ILLIST(L), L = K1, K2)
                  ELSE
                     WRITE(*, 10120) (ILLIST(L), L = K1, K2)
                  END IF
                  IF (K2 .LT. IFLINE(K)+NLPS(K)-1) THEN
                     K1 = K2+1
                     GO TO 150
                  END IF
               END IF
  160       CONTINUE
         ELSE
            CALL MESAGE (' ')
            CALL MESAGE ('*----------------------------------*')
            CALL MESAGE ('* NO SIDES IN THE CURRENT DATABASE *')
            CALL MESAGE ('*----------------------------------*')
         END IF
  170    CONTINUE

C  SPAWN A PROCESS

      ELSE IF ((CIN(ICOM)(1:2) .EQ. 'SP') .OR.
     &   (CIN(ICOM)(1:2) .EQ. 'sp')) THEN
         ICOM = ICOM+1
         CALL SPAWN(VAXVMS)

C  LIST OUT SCHEMES

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'S') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 's')) THEN
         ICOM = ICOM+1
         IF (N(10) .GT. 0) THEN
            CALL MESAGE (' ')
            CALL INTRUP ('LIST ALL SCHEMES', IANS, MCOM, ICOM, JCOM,
     &         CIN, IIN, RIN, KIN)
            IF (IANS) THEN
               I1 = 1
               I2 = N(24)
            ELSE
               CALL MESAGE ('LIST SCHEMES <I1> THROUGH <I2>:')
               IF (ICOM .GT. JCOM) THEN
                  CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN,
     &               CIN, IIN, RIN)
                  ICOM = 1
               END IF
               CALL GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2,
     &            IFOUND)
               IF (IFOUND .GT. 0) THEN
                  CALL CHECK(I1, I2, N(24))
               ELSE
                  GO TO 190
               END IF
            END IF
            CALL MESAGE (' ')
            WRITE(*, 10170)
            DO 180 I = I1, I2
               CALL LTSORT (MR, LINKSC, I, K, ADDLNK)
               IF (K .GT. 0) THEN
                  WRITE(*, 10190)ISCHM(K), SCHEME(K)
               END IF
  180       CONTINUE
         ELSE
            CALL MESAGE (' ')
            CALL MESAGE
     &         ('*---------------------------------------------*')
            CALL MESAGE
     &         ('* ONLY DEFAULT SCHEME IN THE CURRENT DATABASE *')
            CALL MESAGE
     &         ('*---------------------------------------------*')
            WRITE(*, 10180)DEFSCH
            CALL MESAGE (' ')
         END IF
  190    CONTINUE

C  LIST OUT THE BAR SETS

      ELSE IF ((CIN(ICOM)(1:2) .EQ. 'BA') .OR.
     &   (CIN(ICOM)(1:2) .EQ. 'ba')) THEN
         ICOM = ICOM+1
         IF (N(5) .GT. 0) THEN
            CALL MESAGE (' ')
            CALL INTRUP ('LIST ALL BAR SETS', IANS, MCOM, ICOM, JCOM,
     &         CIN, IIN, RIN, KIN)
            IF (IANS) THEN
               I1 = 1
               I2 = N(21)
            ELSE
               CALL MESAGE ('LIST BAR SETS <I1> THROUGH <I2>:')
               IF (ICOM .GT. JCOM) THEN
                  CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN,
     &               CIN, IIN, RIN)
                  ICOM = 1
               END IF
               CALL GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2,
     &            IFOUND)
               IF (IFOUND .GT. 0) THEN
                  CALL CHECK(I1, I2, N(21))
               ELSE
                  GO TO 220
               END IF
            END IF
            CALL MESAGE (' ')
            WRITE(*, 10130)
            DO 210 I = I1, I2
               CALL LTSORT (MS, LINKB, I, K, ADDLNK)
               IF (K .GT. 0) THEN
                  K1 = JFLINE(K)
  200             CONTINUE
                  K2 = MIN0(K1+10, JFLINE(K)+NLPB(K)-1)
                  IF (K1 .EQ. JFLINE(K)) THEN
                     WRITE(*, 10140)IBARST(K), JMAT(K), JCENT(K),
     &                  (JLLIST(L), L = K1, K2)
                  ELSE
                     WRITE(*, 10150) (JLLIST(L), L = K1, K2)
                  END IF
                  IF (K2 .LT. JFLINE(I)+NLPB(K)-1) THEN
                     K1 = K2+1
                     GO TO 200
                  END IF
               END IF
  210       CONTINUE
         ELSE
            CALL MESAGE (' ')
            CALL MESAGE ('*-------------------------------------*')
            CALL MESAGE ('* NO BAR SETS IN THE CURRENT DATABASE *')
            CALL MESAGE ('*-------------------------------------*')
         END IF
  220    CONTINUE

C  LIST OUT THE RENUMBERING CARDS

      ELSE IF ((CIN(ICOM)(1:3) .EQ. 'REN') .OR.
     &   (CIN(ICOM)(1:3) .EQ. 'ren')) THEN
         ICOM = ICOM+1
         IF (N(28) .GT. 0) THEN
            CALL MESAGE (' ')
            I1 = 1
            I2 = N(28)
            CALL MESAGE (' ')
            WRITE(*, 10220)
            DO 230 I = I1, I2
               WRITE(*, 10230)I, NUMBER(I)(1:72)
               IF (NUMBER(I) (73:80) .NE. '        ') THEN
                  WRITE(*, 10240)NUMBER(I) (73:80)
               END IF
  230       CONTINUE
         ELSE IF (OPTIM) THEN
            CALL MESAGE (' ')
            CALL MESAGE ('*------------------------------------------*')
            CALL MESAGE ('* NO RENUMBER CARDS - OPTIMIZATION ENABLED *')
            CALL MESAGE ('*------------------------------------------*')
            CALL MESAGE (' ')
         ELSE
            CALL MESAGE (' ')
            CALL MESAGE
     &         ('*-------------------------------------------*')
            CALL MESAGE
     &         ('* NO RENUMBER CARDS - OPTIMIZATION DISABLED *')
            CALL MESAGE
     &         ('*-------------------------------------------*')
            CALL MESAGE (' ')
         END IF

C  LIST OUT THE REGIONS

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'R') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'r')) THEN
         ICOM = ICOM+1
         LGROUP = .TRUE.
         DO 240 I = 1, N(7)
            IF (IRGFLG(I) .LE. -1) THEN
               LGROUP = .FALSE.
               GO TO 250
            END IF
  240    CONTINUE
  250    CONTINUE
         IF (.NOT.LGROUP) THEN
            CALL MESAGE (' ')
            CALL INTRUP ('LIST ALL REGIONS', IANS, MCOM, ICOM, JCOM,
     &         CIN, IIN, RIN, KIN)
            IF (IANS) THEN
               I1 = 1
               I2 = N(22)
            ELSE
               CALL MESAGE ('LIST REGIONS <I1> THROUGH <I2>:')
               IF (ICOM .GT. JCOM) THEN
                  CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN,
     &               CIN, IIN, RIN)
                  ICOM = 1
               END IF
               CALL GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2,
     &            IFOUND)
               IF (IFOUND .GT. 0) THEN
                  CALL CHECK(I1, I2, N(22))
               ELSE
                  GO TO 280
               END IF
            END IF
            CALL MESAGE (' ')
            WRITE(*, 10050)
            DO 270 I = I1, I2
               CALL LTSORT (MR, LINKR, I, K, ADDLNK)
               IF ((K .GT. 0) .AND. (IRGFLG(K) .LE. 0)) THEN
                  IF (RSIZE(K) .GT. 0) THEN
                     RRSIZE = RSIZE(K)
                  ELSE
                     RRSIZE = DEFSIZ
                  END IF
                  CALL LTSORT (MR, LINKSC, K, IPNTR, ADDLNK)
                  IF ((N(24) .GE. I) .AND. (IPNTR .GT. 0)) THEN
                     DUMMY = SCHEME(IPNTR)(1:10)
                  ELSE
                     DUMMY = DEFSCH(1:10)
                  END IF
                  K1 = IFSIDE(K)
  260             CONTINUE
                  K2 = MIN0(K1+4, IFSIDE(K)+NSPR(K)-1)
                  IF (K1 .EQ. IFSIDE(K)) THEN
                     WRITE(*, 10060)IREGN(K), IMAT(K), RRSIZE, DUMMY,
     &                  (ISLIST(L), L = K1, K2)
                  ELSE
                     WRITE(*, 10070) (ISLIST(L), L = K1, K2)
                  END IF
                  IF (K2 .LT. IFSIDE(K)+NSPR(K)-1) THEN
                     K1 = K2+1
                     GO TO 260
                  END IF
               END IF
  270       CONTINUE
         ELSE
            CALL MESAGE (' ')
            CALL MESAGE ('*------------------------------------*')
            CALL MESAGE ('* NO REGIONS IN THE CURRENT DATABASE *')
            CALL MESAGE ('*------------------------------------*')
         END IF
  280    CONTINUE

C  LIST OUT THE GROUPS

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'G') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'g')) THEN
         ICOM = ICOM+1
         LGROUP = .FALSE.
         DO 290 I = 1, N(7)
            IF (IRGFLG(I) .GE. 1) THEN
               LGROUP = .TRUE.
               GO TO 300
            END IF
  290    CONTINUE
  300    CONTINUE
         IF (LGROUP) THEN
            CALL MESAGE (' ')
            CALL INTRUP ('LIST ALL GROUPS', IANS, MCOM, ICOM, JCOM,
     &         CIN, IIN, RIN, KIN)
            IF (IANS) THEN
               I1 = 1
               I2 = N(22)
            ELSE
               CALL MESAGE ('LIST GROUPS <I1> THROUGH <I2>:')
               IF (ICOM .GT. JCOM) THEN
                  CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN,
     &               CIN, IIN, RIN)
                  ICOM = 1
               END IF
               CALL GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2,
     &            IFOUND)
               IF (IFOUND .GT. 0) THEN
                  CALL CHECK(I1, I2, N(22))
               ELSE
                  GO TO 330
               END IF
            END IF
            CALL MESAGE (' ')
            WRITE(*, 10080)
            DO 320 I = I1, I2
               CALL LTSORT (MR, LINKR, I, K, ADDLNK)
               IF ((K .GT. 0) .AND. (IRGFLG(K) .GE. 1)) THEN
                  K1 = IFSIDE(K)
  310             CONTINUE
                  K2 = MIN0(K1+10, IFSIDE(K)+NSPR(K)-1)
                  IF (K1 .EQ. IFSIDE(K)) THEN
                     WRITE(*, 10110) IREGN(K), (ISLIST(L), L = K1, K2)
                  ELSE
                     WRITE(*, 10120) (ISLIST(L), L = K1, K2)
                  END IF
                  IF (K2 .LT. IFSIDE(K)+NSPR(K)-1) THEN
                     K1 = K2+1
                     GO TO 310
                  END IF
               END IF
  320       CONTINUE
         ELSE
            CALL MESAGE (' ')
            CALL MESAGE ('*-----------------------------------*')
            CALL MESAGE ('* NO GROUPS IN THE CURRENT DATABASE *')
            CALL MESAGE ('*-----------------------------------*')
         END IF
  330    CONTINUE

C  LIST OUT THE HOLES

      ELSE IF ((CIN(ICOM)(1:2) .EQ. 'HO') .OR.
     &   (CIN(ICOM)(1:2) .EQ. 'ho')) THEN
         ICOM = ICOM+1
         IF (N(29) .GT. 0) THEN
            CALL MESAGE (' ')
            CALL INTRUP ('LIST ALL REGIONS WITH HOLES', IANS, MCOM,
     &         ICOM, JCOM, CIN, IIN, RIN, KIN)
            IF (IANS) THEN
               I1 = 1
               I2 = N(7)
            ELSE
               CALL MESAGE ('LIST HOLES IN REGIONS <I1> THROUGH <I2>:')
               IF (ICOM .GT. JCOM) THEN
                  CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN,
     &               CIN, IIN, RIN)
                  ICOM = 1
               END IF
               CALL GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2,
     &            IFOUND)
               IF (IFOUND .GT. 0) THEN
                  CALL CHECK(I1, I2, N(7))
               ELSE
                  GO TO 360
               END IF
            END IF
            CALL MESAGE (' ')
            WRITE(*, 10100)
            DO 350 I = I1, I2
               CALL LTSORT (MR, LINKR, I, K, ADDLNK)
               IF (K .GT. 0) THEN
                  K1 = IFHOLE(K)
  340             CONTINUE
                  K2 = MIN0(K1+10, IFHOLE(K)+NHPR(K)-1)
                  IF (K2 .GE. K1) THEN
                     IF (K1 .EQ. IFHOLE(K)) THEN
                        WRITE(*, 10110)IREGN(K), (IHLIST(L), L = K1, K2)
                     ELSE
                        WRITE(*, 10120) (IHLIST(L), L = K1, K2)
                     END IF
                     IF (K2 .LT. IFHOLE(K)+NHPR(K)-1) THEN
                        K1 = K2+1
                        GO TO 340
                     END IF
                  END IF
               END IF
  350       CONTINUE
         ELSE
            CALL MESAGE (' ')
            CALL MESAGE ('*----------------------------------*')
            CALL MESAGE ('* NO HOLES IN THE CURRENT DATABASE *')
            CALL MESAGE ('*----------------------------------*')
         END IF
  360    CONTINUE

C  LIST OUT THE REGIONS IN THE BODY

      ELSE IF ((CIN(ICOM)(1:3) .EQ. 'BOD') .OR.
     &   (CIN(ICOM)(1:3) .EQ. 'bod')) THEN
         ICOM = ICOM+1
         J1 = 1
         WRITE(*, 10200)
  370    CONTINUE
         IF ((N(9)-J1+1) .GT. 13) THEN
            J2 = J1+12
            WRITE(*, 10210) (IRPB(J), J = J1, J2)
            J1 = J2+1
            GO TO 370
         END IF
         WRITE(*, 10210) (IRPB(J), J = J1, N(9))
      ELSE IF (CIN(ICOM)(1:1) .EQ. ' ') THEN
         ICOM = ICOM+1
         CALL MESAGE (' ')
         RETURN

C  LIST OUT BOUNDARY CONDITIONS

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'B') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'b')) THEN
         ICOM = ICOM+1
         CALL MESAGE (' ')
         CHOICE = 'POINT  '
         CALL LISTBF(MP, N(25), CHOICE, LINKPB, IPBF, NPPF, IFPB,
     &      LISTPB, IWTPBF)
         CHOICE = 'NODE   '
         CALL LISTBF(ML, N(26), CHOICE, LINKLB, ILBF, NLPF, IFLB,
     &      LISTLB, IWTLBF)
         CHOICE = 'ELEMENT'
         CALL LISTBF(ML, N(27), CHOICE, LINKSB, ISBF, NSPF, IFSB,
     &      LISTSB, IWTSBF)

C  LIST OUT THE THREE NODE QUAD GENERATION FLAG

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'T') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 't')) THEN
         ICOM = ICOM+1
         IF (THREE) THEN
            CALL MESAGE (' ')
            CALL MESAGE ('*--------------------------------------*')
            CALL MESAGE ('* THREE NODE BAR GENERATION - ENABLED  *')
            CALL MESAGE ('*--------------------------------------*')
            CALL MESAGE (' ')
         ELSE
            CALL MESAGE (' ')
            CALL MESAGE ('*---------------------------------------*')
            CALL MESAGE ('* THREE NODE BAR GENERATION - DISABLED  *')
            CALL MESAGE ('*---------------------------------------*')
            CALL MESAGE (' ')
         END IF

C  LIST OUT THE EIGHT NODE QUAD GENERATION FLAG

      ELSE IF ((CIN(ICOM)(1:2) .EQ. 'EI') .OR.
     &   (CIN(ICOM)(1:2) .EQ. 'ei')) THEN
         ICOM = ICOM+1
         IF (EIGHT) THEN
            CALL MESAGE (' ')
            CALL MESAGE ('*--------------------------------------*')
            CALL MESAGE ('* EIGHT NODE QUAD GENERATION - ENABLED *')
            CALL MESAGE ('*--------------------------------------*')
            CALL MESAGE (' ')
         ELSE
            CALL MESAGE (' ')
            CALL MESAGE ('*---------------------------------------*')
            CALL MESAGE ('* EIGHT NODE QUAD GENERATION - DISABLED *')
            CALL MESAGE ('*---------------------------------------*')
            CALL MESAGE (' ')
         END IF

C  LIST OUT THE NINE NODE QUAD GENERATION FLAG

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'N') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'n')) THEN
         ICOM = ICOM+1
         IF (NINE) THEN
            CALL MESAGE (' ')
            CALL MESAGE ('*-------------------------------------*')
            CALL MESAGE ('* NINE NODE QUAD GENERATION - ENABLED *')
            CALL MESAGE ('*-------------------------------------*')
            CALL MESAGE (' ')
         ELSE
            CALL MESAGE (' ')
            CALL MESAGE ('*--------------------------------------*')
            CALL MESAGE ('* NINE NODE QUAD GENERATION - DISABLED *')
            CALL MESAGE ('*--------------------------------------*')
            CALL MESAGE (' ')
         END IF

C  EXIT OPTION - EXITS FASTQ

      ELSE IF ((CIN(ICOM)(1:2) .EQ. 'EX') .OR.
     &   (CIN(ICOM)(1:2) .EQ. 'ex')) THEN
         ICOM = ICOM+1
         CALL STRLNG (CIN(ICOM), LEN)
         IF (((LEN .GT. 1) .AND. (CIN(ICOM)(2:2) .NE. 'X')) .OR.
     &      ((LEN .GT. 1) .AND. (CIN(ICOM)(2:2) .NE. 'x'))) THEN
            CALL HELP_FQ(4)
         ELSE
            CALL FEXIT(WROTE, MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN,
     &         TIME1, BATCH, VERSN)
         ENDIF
         GO TO 100

C  PRINT HELP MESAGE

      ELSE
         ICOM = ICOM+1
         CALL HELP_FQ(4)

      END IF
      GO TO 100

10000 FORMAT(
     &   '  POINT            X(R)               Y(Z)        BOUNDARY',/,
     &   '   NO.         COORDINATE         COORDINATE        FLAG',/,
     &   '  -----       ------------      -------------      --------')
10010 FORMAT(1X, I5, 3X, F15.4, 4X, F15.4, 4X, I5)
10020 FORMAT('   LINE  BEGINNING   ENDING   CENTER  NO. OF   NODE BC  ',
     &   'ELEM BC',/,
     &   '    NO.    NODE       NODE     NODE  INTERVALS   FLAG   ',
     &   '  FLAG    FACTOR',/,
     &   '  -----  ---------   ------   ------ --------- -------- ',
     &   '-------- --------')
10030 FORMAT(1X, 7(I5, 4X), F10.5)
10040 FORMAT(1X, 3(I5, 4X), '  -----', 2X, 3(I5, 4X), F10.5)
10050 FORMAT(' REGION  MAT.  INTERVAL    REGION',/,
     &   '   NO.    NO.    SIZE      SCHEME       REGION SIDE/',
     &   'LINE LISTING',/,
     &   ' ------  ---- ---------- -----------  ',
     &   '------------------------------------')
10060 FORMAT(I5, 2X, I5, 2X, F10.5, 2X, A10, 2X, 7I6)
10070 FORMAT(38X, 9I6)
10080 FORMAT(' GROUP'/,
     &   '   NO.                 REGION LISTING',/,
     &   ' ------  ----------------------------------------------')
10090 FORMAT('  SIDE',/,
     &   '   NO.               SIDE/LINE LISTING', /,
     &   ' ------  ----------------------------------------------')
10100 FORMAT(' REGION',/,
     &   '   NO.              HOLE REGION LISTING', /,
     &   ' ------  ----------------------------------------------')
10110 FORMAT(I5, 5X, 11I6)
10120 FORMAT(10X, 11I6)
10130 FORMAT('  BAR SET  MAT.  REFR     BAR SET LINE', /,
     &   '    NO.    NO.   NODE        LISTING',/,
     &   ' -------- ------ ------ ',
     &   '-----------------------------------')
10140 FORMAT(I6, 2X, I5, 2X, I5, 2X, 11I5)
10150 FORMAT(22X, 11I5)
10170 FORMAT('  FOR', /, ' REGION                  SCHEME', /,
     &   ' ------  ----------------------------------------------')
10180 FORMAT(' DEFLT: ', A72)
10190 FORMAT(1X, I5, 2X, A72)
10200 FORMAT('     REGIONS IN THE BODY',/,
     &   ' ----------------------------------------------------')
10210 FORMAT(1X, 13I6)
10220 FORMAT(' CARD           RENUMBERING CARD', /, '  NO.', /,
     &   ' -----  -------------------------------------------')
10230 FORMAT(1X, I5, 2X, A72)
10240 FORMAT(8X, A8)

      END
