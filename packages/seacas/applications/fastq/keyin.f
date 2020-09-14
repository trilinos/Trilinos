C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE KEYIN (MP, ML, MS, MR, MSC, MA, MCOM, ICOM, JCOM, CIN,
     &   RIN, IIN, KIN, IDUMP, N, IPOINT, COOR, IPBOUN, ILINE, LTYPE,
     &   NINT, FACTOR, LCON, ILBOUN, ISBOUN, ISIDE, NLPS, IFLINE,
     &   ILLIST, IBARST, JMAT, JCENT, NLPB, JFLINE, JLLIST, IREGN, IMAT,
     &   NSPR, IFSIDE, ISLIST, IRPB, IPBF, NPPF, IFPB, LISTPB, ILBF,
     &   NLPF, IFLB, LISTLB, ISBF, NSPF, IFSB, LISTSB, LINKP, LINKL,
     &   LINKS, LINKB, LINKR, LINKM, LINKSC, LINKPB, LINKLB, LINKSB,
     &   IHOLDP, IHOLDL, IHOLDS, IHOLDB, IHOLDR, IHOLDM, IHOLD1, IHOLD2,
     &   IHOLD3, IWTPBF, IWTLBF, IWTSBF, RSIZE, IFHOLE, NHPR, IHLIST,
     &   IRGFLG, ISCHM, SCHEME, NUMBER, DEFSCH, DEFSIZ, TITLE, OPTIM,
     &   THREE, EIGHT, NINE, NOROOM, VAXVMS, WROTE, TIME1, VERSN, BATCH)
C***********************************************************************

C  SUBROUTINE KEYIN = INPUTS MESH DEFINITIONS FROM THE KEYBOARD

C***********************************************************************

C  SUBROUTINE CALLED BY:
C     FASTQ = A PROGRAM TO QUICKLY GENERATE QUADRILATERAL MESHES

C***********************************************************************

      PARAMETER (NIN = 80)

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
      DIMENSION LINKP(2, MP), LINKL(2, ML), LINKS(2, MS), LINKB(2, MS)
      DIMENSION LINKR(2, MR), LINKM(2, (MS + MR)), LINKSC(2, MR)
      DIMENSION LINKPB(2, MP), LINKLB(2, ML), LINKSB(2, ML)
      DIMENSION IHOLDP(2, MP), IHOLDL(2, ML), IHOLDS(2, MS)
      DIMENSION IHOLDR(2, MR), IHOLDB(2, MS), IHOLDM(2, (MS + MR))
      DIMENSION IHOLD1(2, MP), IHOLD2(2, ML), IHOLD3(2, ML)
      DIMENSION IFHOLE(MR), NHPR(MR), IHLIST(MR*2), IRGFLG(MR)
      DIMENSION NUMBER(MSC)
      DIMENSION N(29), NOLD(29), III(1)
      DIMENSION KIN(MCOM), IIN(MCOM), RIN(MCOM), JIN(NIN)

      CHARACTER*72 SCHEME, DEFSCH, CIN(MCOM), VERSN*9
      CHARACTER*72 TITLE, HOLD, NUMBER*80

      LOGICAL IANS, OPTIM, NOROOM, ADDOLD, MERGE, NEWNUM, DOLINK, ADDLNK
      LOGICAL THREE, EIGHT, NINE, VAXVMS, WROTE, SIDEOK, BATCH

      IZ = 0
      MERGE = .FALSE.
      DOLINK = .TRUE.
      NOROOM = .FALSE.
      ADDLNK = .FALSE.

      DO 100 I = 1, 29
         NOLD(I) = N(I)
  100 CONTINUE

  110 CONTINUE
      IF (ICOM .GT. JCOM) THEN
         CALL MESAGE (' ')
         CALL FREFLD (IZ, IZ, 'ENTER KEYIN OPTION: ', MCOM, IOSTAT,
     &      JCOM, KIN, CIN, IIN, RIN)
         ICOM = 1
      END IF

C  INPUT A POINT INTO THE DATABASE

      IF ((CIN(ICOM)(1:1) .EQ. 'P') .OR. (CIN(ICOM)(1:1) .EQ. 'p')) THEN
         ICOM = ICOM + 1
         CALL MESAGE ('ENTER POINT DATA IN THE FOLLOWING FORMAT:')
         CALL MESAGE ('[ POINT NO., X, Y ]')
         CALL MESAGE ('HIT RETURN TO END INPUT')
         ICOM = JCOM + 1
  120    CONTINUE
         CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, IFOUND, KIN, CIN, IIN,
     &      RIN)
         IF (IFOUND .GT. 0) THEN
            JJ = IIN(1)
            IF ((JJ .LE. 0) .OR. (JJ .GE. 10000)) THEN
               WRITE(*, 10000)JJ
               GO TO 120
            END IF
            ADDLNK = .FALSE.
            CALL LTSORT(MP, LINKP, JJ, IPNTR, ADDLNK)
            CALL INPOIN(MP, N(1), N(18), JJ, RIN(2), RIN(3), NHOLDP,
     &         IHOLDP, IPOINT, COOR, IPBOUN, LINKP, MERGE, NOROOM)
            IF (NOROOM) GO TO 400

C  REPLACE THE FLAGS OF A REDEFINED POINT

            IF (IPNTR .GT. 0) THEN
               CALL LTSORT(MP, LINKP, JJ, JPNTR, ADDLNK)
               IPBOUN(JPNTR) = IPBOUN(IPNTR)
            END IF
            GO TO 120
         END IF

C  ENTER A LINE INTO THE DATABASE

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'L') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'l')) THEN
         ICOM = ICOM + 1
         IF (ICOM .GT. JCOM) THEN
            CALL MESAGE (' ')
            CALL MESAGE ('THE FOLLOWING LINE TYPES ARE AVAILABLE:')
            CALL MESAGE ('     S*TRAIGHT = STRAIGHT LINE')
            CALL MESAGE ('     CI*RCULAR = CIRCULAR CCW ARC ABOUT A '//
     &         'CENTER')
            CALL MESAGE ('     3*CIRCULAR = CIRCULAR ARC WITH 3RD '//
     &         'ARC POINT')
            CALL MESAGE ('     R*CIRCULAR = CIRCULAR ARC WITH RADIUS')
            CALL MESAGE ('     E*LIPSE    = CCW ELIPSE ABOUT A CENTER')
            CALL MESAGE ('     CO*RNER    = 2 LINE SEGMENTS JOINED')
            CALL MESAGE ('     P*ARABOLA  = PARABOLIC SHAPED LINE')
            CALL FREFLD (IZ, IZ, 'WHICH LINE TYPE WOULD YOU LIKE TO '//
     &         'ENTER:', MCOM, IOSTAT, JCOM, KIN, CIN, IIN, RIN)
            ICOM = 1
         END IF
         IF ((CIN(ICOM)(1:1) .EQ. 'S') .OR.
     &      (CIN(ICOM)(1:1) .EQ. 's')) THEN
            ICOM = ICOM + 1
            CALL MESAGE ('ENTER STRAIGHT LINE DATA IN THE FOLLOWING '//
     &         'FORMAT:')
            CALL MESAGE ('[ LINE NO., POINT 1, POINT 2, NO. '//
     &         'INTERVALS, FACTOR ]')
            IT = 1
         ELSE IF ((CIN(ICOM)(1:2) .EQ. 'CI') .OR.
     &      (CIN(ICOM)(1:2) .EQ. 'ci')) THEN
            ICOM = ICOM + 1
            CALL MESAGE ('NOTE: IF A CW ARC IS DESIRED, ENTER')
            CALL MESAGE ('      THE CENTER POINT AS NEGATIVE')
            CALL MESAGE ('ENTER CIRCULAR ARC LINE DATA IN THE '//
     &         'FOLLOWING FORMAT:')
            CALL MESAGE ('[ LINE NO., POINT 1, POINT 2, CENTER, NO. '//
     &         'INTERVALS, FACTOR ]')
            IT = 3
         ELSE IF (CIN(ICOM)(1:1) .EQ. '3') THEN
            ICOM = ICOM + 1
            CALL MESAGE ('ENTER THIRD POINT ARC DATA IN THE '//
     &         'FOLLOWING FORMAT:')
            CALL MESAGE ('[ LINE NO., POINT 1, POINT 2, POINT 3, NO. '//
     &         'INTERVALS, FACTOR ]')
            IT = 4
         ELSE IF ((CIN(ICOM)(1:1) .EQ. 'R') .OR.
     &      (CIN(ICOM)(1:1) .EQ. 'r')) THEN
            ICOM = ICOM + 1
            CALL MESAGE ('NOTE: THE RADIUS IS ASSUMED TO BE CONTAINED')
            CALL MESAGE ('      IN THE X COORDINATE OF POINT 3. THE')
            CALL MESAGE ('      ARC CENTER IS ASSUMED TO THE LEFT OF A')
            CALL MESAGE ('      LINE FROM POINT 1 TO POINT 2 (OPPOSITE')
            CALL MESAGE ('      IF POINT 3 IS ENTERED NEGATIVE).')
            CALL MESAGE ('ENTER CIRCULAR ARC W/RADIUS DATA IN THE '//
     &         'FOLLOWING FORMAT:')
            CALL MESAGE ('[ LINE NO., POINT 1, POINT 2, CENTER, NO. '//
     &         'INTERVALS, FACTOR ]')
            IT = 7
         ELSE IF ((CIN(ICOM)(1:1) .EQ. 'E') .OR.
     &      (CIN(ICOM)(1:1) .EQ. 'e')) THEN
            ICOM = ICOM + 1
            CALL MESAGE ('NOTE: THE TWO POINTS ON THE ELIPSE CANNOT'//
     &         ' BE COLINEAR WITH THE')
            CALL MESAGE ('      CENTER POINT IN THIS DEFINITION.')
            CALL MESAGE ('NOTE: IF A CW ARC IS DESIRED, ENTER')
            CALL MESAGE ('      THE CENTER POINT AS NEGATIVE.')
            CALL MESAGE ('ENTER ELIPSE ABOUT A CENTER DATA IN THE '//
     &         'FOLLOWING FORMAT:')
            CALL MESAGE ('[ LINE NO., POINT 1, POINT 2, POINT 3, NO. '//
     &         'INTERVALS, FACTOR ]')
            IT = 7
         ELSE IF ((CIN(ICOM)(1:2) .EQ. 'CO') .OR.
     &      (CIN(ICOM)(1:2) .EQ. 'co')) THEN
            ICOM = ICOM + 1
            CALL MESAGE ('NOTE: A CORNER LINE CONTAINS TO STRAIGHT '//
     &         'LINE')
            CALL MESAGE ('      SEGMENTS JOINED AT POINT 3')
            CALL MESAGE ('ENTER CORNER LINE DATA IN THE FOLLOWING '//
     &         'FORMAT:')
            CALL MESAGE ('[ LINE NO., POINT 1, POINT 2, POINT 3, NO. '//
     &         'INTERVALS, FACTOR ]')
            IT = 2
         ELSE IF ((CIN(ICOM)(1:1) .EQ. 'P') .OR.
     &      (CIN(ICOM)(1:1) .EQ. 'p')) THEN
            ICOM = ICOM + 1
            CALL MESAGE ('NOTE: POINT 3 IS THE TIP OF THE PARABOLA, '//
     &         'AND')
            CALL MESAGE ('      POINT 1 AND POINT 2 MUST BE EQUAL ARC')
            CALL MESAGE ('      LENGTHS AWAY (ISOCELES TRIANGLE)')
            CALL MESAGE ('ENTER PARABOLIC LINE DATA IN THE FOLLOWING '//
     &         'FORMAT:')
            CALL MESAGE ('[ LINE NO., POINT 1, POINT 2, POINT 3, NO. '//
     &         'INTERVALS, FACTOR ]')
            IT = 5
         ELSE
            ICOM = ICOM + 1
            GO TO 110
         END IF
         CALL MESAGE ('HIT RETURN TO END INPUT')
         ICOM = JCOM + 1
  130    CONTINUE
         CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, IFOUND, KIN, CIN, IIN,
     &      RIN)
         IF (IFOUND .GT. 0) THEN
            JJ = IIN(1)
            IF ((JJ .LE. 0) .OR. (JJ .GE. 10000)) THEN
               WRITE(*, 10010)JJ
               GO TO 130
            ELSE IF (((IT .EQ. 1).AND.(IFOUND .LT. 3)) .OR.
     &         ((IT.NE.1).AND.(IFOUND .LT. 4))) THEN
               WRITE(*, 10020)JJ
               GO TO 130
            END IF
            IZERO = 0
            ADDLNK = .FALSE.
            CALL LTSORT(ML, LINKL, JJ, IPNTR, ADDLNK)
            IF (IT .EQ. 1) THEN
               CALL INLINE(ML, N(2), N(19), JJ, IT, IIN(2), IIN(3),
     &            IZERO, IIN(4), RIN(5), NHOLDL, IHOLDL, ILINE, LTYPE,
     &            NINT, FACTOR, LCON, ILBOUN, ISBOUN, LINKL, MERGE,
     &            NOROOM)
            ELSE
               CALL INLINE(ML, N(2), N(19), JJ, IT, IIN(2), IIN(3),
     &            IIN(4), IIN(5), RIN(6), NHOLDL, IHOLDL, ILINE, LTYPE,
     &            NINT, FACTOR, LCON, ILBOUN, ISBOUN, LINKL, MERGE,
     &            NOROOM)
            END IF
            IF (NOROOM) GO TO 400

C  LINK UP THE OLD FLAGS TO THE NEW LINE

            IF (IPNTR .GT. 0) THEN
               ADDLNK = .FALSE.
               CALL LTSORT(ML, LINKL, JJ, JPNTR, ADDLNK)
               ILBOUN(JPNTR) = ILBOUN(IPNTR)
               ISBOUN(JPNTR) = ISBOUN(IPNTR)
            END IF
            GO TO 130
         END IF

C  ENTER A REGION INTERVAL SIZE INTO THE DATABASE

      ELSE IF ((CIN(ICOM)(1:3) .EQ. 'SIZ') .OR.
     &   (CIN(ICOM)(1:3) .EQ. 'siz')) THEN
         ICOM = ICOM + 1
         CALL MESAGE ('NOTE: ENTER A DEFAULT SIZE BY SPECIFYING')
         CALL MESAGE ('      A SIZE WITH NO REGIONS.')
         CALL MESAGE ('ENTER REGION SIZE DATA IN THE FOLLOWING FORMAT:')
         CALL MESAGE ('[ SIZE, REGION 1, REGION 2, ..., REGION N ]')
         CALL MESAGE ('HIT RETURN TO END INPUT')
         ICOM = JCOM + 1
  140    CONTINUE
         CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, IFOUND, KIN, CIN, IIN,
     &      RIN)
         IF (IFOUND .GT. 0) THEN
            IF (IFOUND .LT. 2) THEN
               DEFSIZ = RIN(1)
            ELSE
               ADDLNK = .FALSE.
               DO 150 IRSZ = 2, IFOUND
                  CALL LTSORT(MR, LINKR, IIN(IRSZ), JJ, ADDLNK)
                  IF (JJ .GE. 0) THEN
                     RSIZE(JJ) = AMAX1(RIN(1), 0.)
                  ELSE
                     WRITE(*, 10030)IIN(IRSZ)
                  END IF
  150          CONTINUE
            END IF
            GO TO 140
         END IF

C  ENTER A SIDE INTO THE DATABASE

      ELSE IF ((CIN(ICOM)(1:2) .EQ. 'SI') .OR.
     &   (CIN(ICOM)(1:2) .EQ. 'si')) THEN
         ICOM = ICOM + 1
         CALL MESAGE ('ENTER SIDE DATA IN THE FOLLOWING FORMAT:')
         CALL MESAGE ('[ SIDE NO., LINE 1, LINE 2, ..., LINE N ]')
         CALL MESAGE ('HIT RETURN TO END INPUT')
         ICOM = JCOM + 1
  160    CONTINUE
         CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, IFOUND, KIN, CIN, IIN,
     &      RIN)
         IF (IFOUND .GT. 0) THEN
            JJ = IIN(1)
            IF ((JJ .LE. 0) .OR. (JJ .GE. 10000)) THEN
               WRITE(*, 10040)JJ
               GO TO 160
            ELSE IF (IFOUND .LT. 2) THEN
               WRITE(*, 10070)JJ
               GO TO 160
            END IF
            CALL INSIDE(MS, N(3), N(4), N(20), JJ, IIN(2), IFOUND - 1,
     &         ISIDE, NLPS, IFLINE, ILLIST, LINKS, NHOLDS, IHOLDS,
     &         MERGE, NOROOM)
            IF (NOROOM) GO TO 400
            GO TO 160
         END IF

C  ENTER A HOLE INTO THE DATABASE

      ELSE IF ((CIN(ICOM)(1:2) .EQ. 'HO') .OR.
     &   (CIN(ICOM)(1:2) .EQ. 'ho')) THEN
         ICOM = ICOM + 1
         CALL MESAGE ('ENTER HOLE DATA IN THE FOLLOWING FORMAT:')
         CALL MESAGE ('[ REGION NO., HOLE 1, HOLE 2, ..., HOLE N ]')
         CALL MESAGE ('HIT RETURN TO END INPUT')
         ICOM = JCOM + 1
  170    CONTINUE
         CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, IFOUND, KIN, CIN, IIN,
     &      RIN)
         IF (IFOUND .GT. 0) THEN
            JJ = IIN(1)
            ADDLNK = .FALSE.
            CALL LTSORT (MR, LINKR, JJ, JJPNTR, ADDLNK)
            ADDLNK = .TRUE.
            IF (JJPNTR  .LE.  0) THEN
               WRITE(*, 10080) JJ
               GO TO 170
            ELSE IF (IFOUND  .LT.  1) THEN
               WRITE(*, 10090) JJ
               GO TO 170
            END IF
            CALL INHOLE (MR, N(7), N(29), JJPNTR, IIN(2), IFOUND  -  1,
     &         IFHOLE, NHPR, IHLIST, MERGE, NOROOM)
            IF (NOROOM) GO TO 400
            GO TO 170
         END IF

C  ENTER A BARSET INTO THE DATABASE

      ELSE IF ((CIN(ICOM)(1:2) .EQ. 'BA') .OR.
     &   (CIN(ICOM)(1:2) .EQ. 'ba')) THEN
         ICOM = ICOM + 1
         CALL MESAGE ('ENTER BAR SET DATA IN THE FOLLOWING FORMAT:')
         CALL MESAGE ('[ BAR SET NO., MAT NO., REFR. PNT., LINE 1, '//
     &      'LINE 2, ..., LINE N ]')
         CALL MESAGE ('HIT RETURN TO END INPUT')
         ICOM = JCOM + 1
  180    CONTINUE
         CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, IFOUND, KIN, CIN, IIN,
     &      RIN)
         IF (IFOUND .GT. 0) THEN
            JJ = IIN(1)
            IF ((JJ .LE. 0) .OR. (JJ .GE. 10000)) THEN
               WRITE(*, 10050)JJ
               GO TO 180
            ELSE IF (IFOUND .LT. 4) THEN
               WRITE(*, 10060)JJ
               GO TO 180
            END IF
            CALL INBRST(MS, MR, N(5), N(6), N(21), N(23), JJ, IIN(2),
     &         IIN(3), IIN(4), IFOUND - 3, IBARST, JMAT, JCENT, NLPB,
     &         JFLINE, JLLIST, LINKB, LINKM, NHOLDM, IHOLDM, NHOLDB,
     &         IHOLDB, MERGE, NOROOM)
            IF (NOROOM) GO TO 400
            GO TO 180
         END IF

C  INPUT A BODY DEFINITION INTO THE DATABASE

      ELSE IF ((CIN(ICOM)(1:3) .EQ. 'BOD') .OR.
     &   (CIN(ICOM)(1:3) .EQ. 'bod')) THEN
         ICOM = ICOM + 1
         CALL MESAGE ('NOTE: THE BODY CAN BE MADE UP OF ANY')
         CALL MESAGE ('      COMBINATION OF REGIONS AND BAR SETS.')
         CALL MESAGE ('      ENTER BAR SETS AS NEGATIVE REGIONS.')
         CALL MESAGE ('ENTER REGIONS (& BAR SETS) IN THE BODY IN THE '//
     &      'FOLLOWING FORMAT')
         CALL MESAGE ('[ REGION 1, REGION 2, ..., REGION N ]')
         CALL MESAGE ('HIT A RETURN TO END INPUT')
         ICOM = JCOM + 1
  190    CONTINUE
         CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, IFOUND, KIN, CIN, IIN,
     &      RIN)
         IF ((IFOUND .GT. 0).AND.(IFOUND .LE. NIN)) THEN
            IF (N(9) .GT. 0) THEN
               DO 200 J = 1, IFOUND
                  JIN(J) = IIN(J)
  200          CONTINUE
               CALL INTRUP('REPLACE THE CURRENT BODY DEFINITION', IANS,
     &            MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN)
               DO 210 J = 1, IFOUND
                  IIN(J) = JIN(J)
  210          CONTINUE
               IF (IANS) THEN
                  ADDOLD = .FALSE.
               ELSE
                  ADDOLD = .TRUE.
               END IF
            END IF
            CALL INBODY(MR, N(9), IIN, IFOUND, IRPB, ADDOLD, NOROOM)
            IF (NOROOM) GO TO 400
            GO TO 190
         ELSE IF (IFOUND .GT. NIN) THEN
            CALL MESAGE ('TOO MANY BODIES BEING INPUT A ONCE - TRY '//
     &         'AGAIN')
            GO TO 190
         END IF

C  SPAWN A PROCESS

      ELSE IF ((CIN(ICOM)(1:2) .EQ. 'SP') .OR.
     &   (CIN(ICOM)(1:2) .EQ. 'sp')) THEN
         ICOM = ICOM + 1
         CALL SPAWN(VAXVMS)

C  INPUT A SCHEME INTO THE DATABASE

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'S') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 's')) THEN
         ICOM = ICOM + 1
         CALL MESAGE ('NOTE: ENTER A DEFAULT SCHEME BY SPECIFYING')
         CALL MESAGE ('      THE REGION NUMBER AS ZERO')
         CALL MESAGE ('ENTER A SCHEME IN THE FOLLOWING FORMAT:')
         CALL MESAGE ('[ REGION NO., SCHEME ]')
         CALL MESAGE ('HIT RETURN TO END INPUT')
         ICOM = JCOM + 1
  220    CONTINUE
         CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, IFOUND, KIN, CIN, IIN,
     &      RIN)
         IF (IFOUND .GT. 0) THEN
            JJ = IIN(1)
            IF (JJ .GE. 10000) THEN
               WRITE(*, 10150)JJ
               GO TO 220
            ELSE IF (IFOUND .LT. 2) THEN
               WRITE(*, 10160)JJ
               GO TO 220
            END IF
            NOLD10 = N(10)
            NOLD24 = N(24)
            CALL INSCHM(MR, MSC, N(10), N(24), JJ, CIN(2), ISCHM,
     &         SCHEME, LINKSC, DEFSCH, NOROOM, DOLINK)
            IF (NOROOM) THEN
               N(10) = NOLD10
               N(24) = NOLD24
               CALL MESAGE ('************************************')
               CALL MESAGE ('NOT ENOUGH ROOM FOR SCHEME CARD')
               CALL MESAGE ('NO DYNAMIC DIMENSIONING INCREASES')
               CALL MESAGE ('AVAILABLE FOR CHARACTER STRINGS')
               CALL MESAGE ('SCHEME INPUT IS THUS IGNORED')
               CALL MESAGE ('************************************')
            END IF
            GO TO 220
         END IF

C  INPUT A BOUNDARY INTO THE DATABASE

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'B') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'b')) THEN
         ICOM = ICOM + 1
         IF (ICOM .GT. JCOM) THEN
            CALL MESAGE (' ')
            CALL MESAGE ('THE FOLLOWING BOUNDARY FLAGS ARE AVAILABLE:')
            CALL MESAGE ('        P*OINT FLAGS    - FOR NODES AT '//
     &         'POINTS')
            CALL MESAGE ('        N*ODE FLAGS     - FOR NODES ON A '//
     &         'BOUNDARY')
            CALL MESAGE ('        E*LEMENT FLAGS  - FOR ELEMENT '//
     &         'SIDES ON A BOUNDARY')
            CALL FREFLD (IZ, IZ, 'WHICH BOUNDARY FLAG WOULD YOU LIKE '//
     &         'TO ENTER: ', MCOM, IOSTAT, JCOM, KIN, CIN, IIN, RIN)
            ICOM = 1
         END IF

C  INPUT A POINT BOUNDARY INTO THE DATABASE

         IF ((CIN(ICOM)(1:1) .EQ. 'P') .OR.
     &      (CIN(ICOM)(1:1) .EQ. 'p')) THEN
            ICOM = ICOM + 1
            CALL MESAGE ('INPUT POINT BOUNDARY FLAG DATA IN THE '//
     &         'FOLLOWING FORMAT:')
            CALL MESAGE ('[ FLAG NO., POINT 1, POINT 2, ..., POINT N ]')
            CALL MESAGE ('HIT RETURN TO END INPUT')
            ICOM = JCOM + 1
  230       CONTINUE
            CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, IFOUND, KIN, CIN,
     &         IIN, RIN)
            IF (IFOUND .GT. 0) THEN
               ADDOLD = .FALSE.
               CALL INBOUN(MP, IIN(1), IFOUND - 1, IIN(2), N(25), N(11),
     &            N(12), NOLD(11), MERGE, NOROOM, NEWNUM, NHOLD1,
     &            IHOLD1, IPBF, NPPF, IFPB, LISTPB, LINKPB, IWTPBF,
     &            JHOLD, ADDOLD)
               IF (NOROOM) GO TO 400
               GO TO 230
            END IF

C  INPUT A NODE BOUNDARY INTO THE DATABASE

         ELSE IF ((CIN(ICOM)(1:1) .EQ. 'N') .OR.
     &      (CIN(ICOM)(1:1) .EQ. 'n')) THEN
            ICOM = ICOM + 1
            CALL MESAGE ('INPUT NODE BOUNDARY FLAG DATA IN THE '//
     &         'FOLLOWING FORMAT:')
            CALL MESAGE ('[ FLAG NO., LINE (OR NEG. SIDE) 1, LINE '//
     &         '(OR NEG. SIDE) 2, ...]')
            CALL MESAGE ('HIT RETURN TO END INPUT')
            ICOM = JCOM + 1
  240       CONTINUE
            CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, IFOUND, KIN, CIN,
     &         IIN, RIN)
            IF (IFOUND .GT. 0) THEN
               ADDOLD = .FALSE.
               CALL INBOUN(ML, IIN(1), IFOUND - 1, IIN(2), N(26), N(13),
     &            N(14), NOLD(13), MERGE, NOROOM, NEWNUM, NHOLD2,
     &            IHOLD2, ILBF, NLPF, IFLB, LISTLB, LINKLB, IWTLBF,
     &            JHOLD, ADDOLD)
               IF (NOROOM) GO TO 400
               GO TO 240
            END IF

C  INPUT AN ELEMENT BOUNDARY INTO THE DATABASE

         ELSE IF ((CIN(ICOM)(1:1) .EQ. 'E') .OR.
     &      (CIN(ICOM)(1:1) .EQ. 'e')) THEN
            ICOM = ICOM + 1
            CALL MESAGE ('INPUT ELEMENT BOUNDARY FLAG DATA IN THE '//
     &         'FOLLOWING FORMAT:')
            CALL MESAGE ('[ FLAG NO., LINE (OR NEG. SIDE) 1, LINE '//
     &         '(OR NEG. SIDE) 2, ...]')
            CALL MESAGE ('HIT RETURN TO END INPUT')
            ICOM = JCOM + 1
  250       CONTINUE
            CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, IFOUND, KIN, CIN,
     &         IIN, RIN)
            IF (IFOUND .GT. 0) THEN
               ADDOLD = .FALSE.
               CALL INBOUN(ML, IIN(1), IFOUND - 1, IIN(2), N(27), N(15),
     &            N(16), NOLD(15), MERGE, NOROOM, NEWNUM, NHOLD3,
     &            IHOLD3, ISBF, NSPF, IFSB, LISTSB, LINKSB, IWTSBF,
     &            JHOLD, ADDOLD)
               IF (NOROOM) GO TO 400
               GO TO 250
            END IF
         END IF

C  INPUT A BOUNDARY FLAG WEIGHTING INTO THE DATABASE

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'W') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'w')) THEN
         ICOM = ICOM + 1
         IF (ICOM .GT. JCOM) THEN
            CALL MESAGE (' ')
            CALL MESAGE ('THE FOLLOWING BOUNDARY FLAGS CAN BE '//
     &         'WEIGHTED:')
            CALL MESAGE ('        P*OINT FLAGS    - FOR NODES AT '//
     &         'POINTS')
            CALL MESAGE ('        N*ODE FLAGS     - FOR NODES ON A '//
     &         'BOUNDARY')
            CALL MESAGE ('        E*LEMENT FLAGS  - FOR ELEMENT '//
     &         'SIDES ON A BOUNDARY')
            CALL FREFLD (IZ, IZ, 'WHICH BOUNDARY FLAG WOULD YOU LIKE '//
     &         'TO WEIGHT: ', MCOM, IOSTAT, JCOM, KIN, CIN, IIN, RIN)
            ICOM = 1
         END IF

C  INPUT A POINT BOUNDARY FLAG WEIGHT INTO THE DATABASE

         IF ((CIN(ICOM)(1:1) .EQ. 'P') .OR.
     &      (CIN(ICOM)(1:1) .EQ. 'p')) THEN
            ICOM = ICOM + 1
            CALL MESAGE ('INPUT POINT BOUNDARY FLAG WEIGHTS IN THE '//
     &         'FOLLOWING FORMAT:')
            CALL MESAGE ('[ FLAG NO., WEIGHTING POINT, BOUNDARY '//
     &         'POINT ]')
            CALL MESAGE ('HIT RETURN TO END INPUT')
            ICOM = JCOM + 1
  260       CONTINUE
            CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, IFOUND, KIN, CIN,
     &         IIN, RIN)
            IF (IFOUND .GT. 0) THEN
               CALL LTSORT(MP, LINKPB, IIN(1), JJ, ADDLNK)
               IF (JJ .GT. 0) THEN
                  IWTPBF(1, JJ) = IIN(2)
                  IWTPBF(2, JJ) = IIN(3)
                  IWTPBF(3, JJ) = 0
               ELSE
                  WRITE(*, 10100)'POINT', IIN(1)
               END IF
               GO TO 260
            END IF

C  INPUT A NODE BOUNDARY WEIGHT INTO THE DATABASE

         ELSE IF ((CIN(ICOM)(1:1) .EQ. 'N') .OR.
     &      (CIN(ICOM)(1:1) .EQ. 'n')) THEN
            ICOM = ICOM + 1
            CALL MESAGE ('INPUT NODE BOUNDARY FLAG WEIGHTS IN THE '//
     &         'FOLLOWING FORMAT:')
            CALL MESAGE ('[ FLAG NO., WEIGHTING SIDE (OR NEG. LINE) '//
     &         'NO., BEGINNING POINT NO., ')
            CALL MESAGE ('  BEGINNING LINE NO. (OPTIONAL) ]')
            CALL MESAGE ('HIT RETURN TO END INPUT')
            ICOM = JCOM + 1
  270       CONTINUE
            CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, IFOUND, KIN, CIN,
     &         IIN, RIN)
            IF (IFOUND .GT. 0) THEN
               CALL LTSORT(ML, LINKLB, IIN(1), JJ, ADDLNK)
               IF (JJ .GT. 0) THEN
                  IWTLBF(1, JJ) = IIN(2)
                  IWTLBF(2, JJ) = IIN(3)
                  IWTLBF(3, JJ) = IIN(4)
               ELSE
                  WRITE(*, 10100)'NODE', IIN(1)
               END IF
               GO TO 270
            END IF

C  INPUT AN ELEMENT BOUNDARY WEIGHT INTO THE DATABASE

         ELSE IF ((CIN(ICOM)(1:1) .EQ. 'E') .OR.
     &      (CIN(ICOM)(1:1) .EQ. 'e')) THEN
            ICOM = ICOM + 1
            CALL MESAGE ('INPUT ELEMENT BOUNDARY FLAG WEIGHTS IN THE '//
     &         'FOLLOWING FORMAT:')
            CALL MESAGE ('[ FLAG NO., WEIGHTING SIDE (OR NEG. LINE) NO.,
     &         BEGINNING POINT NO., ')
            CALL MESAGE ('  BEGINNING LINE NO. (OPTIONAL) ]')
            CALL MESAGE ('HIT RETURN TO END INPUT')
            ICOM = JCOM + 1
  280       CONTINUE
            CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, IFOUND, KIN, CIN,
     &         IIN, RIN)
            IF (IFOUND .GT. 0) THEN
               CALL LTSORT(ML, LINKSB, IIN(1), JJ, ADDLNK)
               IF (JJ .GT. 0) THEN
                  IWTSBF(1, JJ) = IIN(2)
                  IWTSBF(2, JJ) = IIN(3)
                  IWTSBF(3, JJ) = IIN(4)
               ELSE
                  WRITE(*, 10100)'ELEMENT', IIN(1)
               END IF
               GO TO 280
            END IF
         END IF

C  TOGGLE THE BANDWIDTH OPTIMIZATION FLAG

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'O') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'o')) THEN
         ICOM = ICOM + 1
         IF (OPTIM) THEN
            OPTIM = .FALSE.
            CALL MESAGE ('BANDWIDTH OPTIMIZER DISABLED')
         ELSE
            OPTIM = .TRUE.
            CALL MESAGE ('BANDWIDTH OPTIMIZER ENABLED')
         END IF

C  FLAG THE BANDWIDTH OPTIMIZATION ROUTINES ON, AND READ A RENUM CARD

      ELSE IF ((CIN(ICOM)(1:3) .EQ. 'REN') .OR.
     &   (CIN(ICOM)(1:3) .EQ. 'ren')) THEN
         ICOM = ICOM + 1
         OPTIM = .TRUE.
         CALL MESAGE (' ')
         CALL MESAGE ('THE FOLLOWING RENUMBERING OPTIONS ARE '//
     &      'AVAILABLE:')
         CALL MESAGE ('        P*-L-P = POINT LINE POINT STARTING LIST')
         CALL MESAGE ('        X*, Y   = X, Y STARTING LOCATION')
         CALL MESAGE ('        N*UID  = NODE UNIQUE ID STARTING LIST')
         IF (ICOM .GT. JCOM) THEN
            CALL FREFLD (IZ, IZ, 'WHICH RENUMBER OPTION WOULD YOU '//
     &         'LIKE TO ENTER: ', MCOM, IOSTAT, JCOM, KIN, CIN, IIN,
     &         RIN)
            ICOM = 1
         END IF

C  ENTER A POINT-LINE-POINT RENUM CARD

         IF ((CIN(ICOM)(1:1) .EQ. 'P') .OR.
     &      (CIN(ICOM)(1:1) .EQ. 'p')) THEN
            ICOM = ICOM + 1
            HOLD = 'P-L-P'
            CALL MESAGE ('ENTER A POINT-LINE-POINT CARD IN THE '//
     &         'FOLLOWING FORMAT')
            CALL MESAGE ('[ ENTER POINT, LINE, POINT, LINE, ... ]')
            CALL MESAGE ('HIT A RETURN TO END INPUT')
            ICOM = JCOM + 1
  290       CONTINUE
            CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, IFOUND, KIN, CIN,
     &         IIN, RIN)
            IF (IFOUND .GE. 1) THEN
               CALL INRENM(MSC, N(28), HOLD, RIN, IIN, IFOUND, NUMBER,
     &            NOROOM)
               IF (NOROOM) GO TO 400
               GO TO 290
            END IF

C  ENTER A X, Y LOCATION RENUM CARD

         ELSE IF ((CIN(ICOM)(1:1) .EQ. 'X') .OR.
     &      (CIN(ICOM)(1:1) .EQ. 'x')) THEN
            ICOM = ICOM + 1
            HOLD = 'X-Y  '
            CALL MESAGE ('ENTER X, Y PAIRS IN THE FOLLOWING FORMAT')
            CALL MESAGE ('[ X , Y ]')
            CALL MESAGE ('HIT A RETURN TO END INPUT')
            ICOM = JCOM + 1
  300       CONTINUE
            CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, IFOUND, KIN, CIN,
     &         IIN, RIN)
            IF (IFOUND .GE. 1) THEN
               CALL INRENM(MSC, N(28), HOLD, RIN, IIN, IFOUND, NUMBER,
     &            NOROOM)
               IF (NOROOM) GO TO 400
               GO TO 300
            END IF

C  ENTER A NODE UNIQUE ID RENUM CARD

         ELSE IF ((CIN(ICOM)(1:1) .EQ. 'N') .OR.
     &      (CIN(ICOM)(1:1) .EQ. 'n')) THEN
            ICOM = ICOM + 1
            HOLD = 'NODE '
            CALL MESAGE ('NOTE: NODE UNIQUE ID NUMBERS (NUID) ARE')
            CALL MESAGE ('      DEFINED IN THE DOCUMENTATION')
            CALL MESAGE ('ENTER NUID CARDS IN THE FOLLOWING FORMAT:')
            CALL MESAGE ('[ NUID 1, NUID 2, ..., NUID N ]')
            CALL MESAGE ('HIT A RETURN TO END INPUT')
            ICOM = JCOM + 1
  310       CONTINUE
            CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, IFOUND, KIN, CIN,
     &         IIN, RIN)
            IF (IFOUND .GE. 1) THEN
               CALL INRENM(MSC, N(28), HOLD, RIN, IIN, IFOUND, NUMBER,
     &            NOROOM)
               IF (NOROOM) GO TO 400
               GO TO 310
            END IF
         END IF

C  ENTER A REGION INTO THE DATABASE

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'R') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'r')) THEN
         ICOM = ICOM + 1
         CALL INTRUP('ARE YOU USING SIDES IN DEFINING REGIONS', IANS,
     &      MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN)
         IF (IANS) THEN
            CALL MESAGE ('NOTE: ENTER ANY LINES NEEDED IN THE REGION '//
     &         'AS')
            CALL MESAGE ('      A NEGATIVE SIDE.')
            CALL MESAGE ('ENTER REGION DATA IN THE FOLLOWING FORMAT:')
            CALL MESAGE ('[ REGION NO., MATERIAL NO., SIDE 1, SIDE '//
     &         '2, ..., SIDE N ]')
         ELSE
            CALL MESAGE ('ENTER REGION DATA IN THE FOLLOWING FORMAT:')
            CALL MESAGE ('[ REGION NO., MATERIAL NO., LINE 1, LINE '//
     &         '2, ..., LINE N ]')
         END IF
         ICOM = JCOM + 1
         CALL MESAGE ('HIT RETURN TO END INPUT')
  320    CONTINUE
         CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, IFOUND, KIN, CIN, IIN,
     &      RIN)
         IF (IFOUND .GT. 0) THEN
            JJ = IIN(1)
            IF ((JJ .LE. 0) .OR. (JJ .GE. 10000)) THEN
               WRITE(*, 10110)JJ
               GO TO 320
            ELSE IF (IFOUND .LT. 3) THEN
               WRITE(*, 10120)JJ
               GO TO 320
            END IF
            IF (.NOT.IANS) THEN
               DO 330 J = 3, IFOUND
                  IIN(J) = -IIN(J)
  330          CONTINUE
            END IF
            CALL INREGN(MS, MR, N(7), N(8), N(22), N(23), JJ, IIN(2),
     &         IIN(3), IFOUND - 2, IREGN, IMAT, NSPR, IFSIDE, ISLIST,
     &         LINKR, LINKM, NHOLDR, IHOLDR, NHOLDM, IHOLDM, IRGFLG,
     &         MERGE, NOROOM)
            IF (NOROOM) GO TO 400
            ADDLNK = .FALSE.
            CALL LTSORT(MR, LINKR, IIN(1), JJPNTR, ADDLNK)
            RSIZE(JJPNTR) = 0.
            GO TO 320
         END IF

C  ENTER A GROUP OF REGIONS INTO THE DATABASE

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'G') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'g')) THEN
         ICOM = ICOM + 1
         CALL MESAGE
     &      ('ENTER GROUP OF REGIONS IN THE FOLLOWING FORMAT:')
         CALL MESAGE
     &      ('[ GROUP NO., REGION 1, REGION 2, ..., REGION N ]')
         ICOM = JCOM + 1
         CALL MESAGE ('HIT RETURN TO END INPUT')
  340    CONTINUE
         CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, IFOUND, KIN, CIN, IIN,
     &      RIN)
         IF (IFOUND .GT. 0) THEN
            JJ = IIN(1)
            IF ((JJ .LE. 0) .OR. (JJ .GE. 10000)) THEN
               WRITE(*, 10130)JJ
               GO TO 340
            ELSE IF (IFOUND .LT. 2) THEN
               WRITE(*, 10140)JJ
               GO TO 340
            END IF
            CALL INGRPN(MS, MR, N(7), N(8), N(22), JJ, IIN(2),
     &         IFOUND - 1, IREGN, NSPR, IFSIDE, ISLIST, LINKR, NHOLDR,
     &         IHOLDR, IRGFLG, MERGE, NOROOM)
            IF (NOROOM) GO TO 400
            GO TO 340
         END IF

C  ENTER A TITLE

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'T') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 't')) THEN
         ICOM = ICOM + 1
         IF (ICOM .GT. JCOM) THEN
            CALL GETINP(IUNIT, IDUMP, 'TITLE: ', TITLE, IOSTAT)
         END IF

C  ENTER LINE INTERVALS

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'I') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'i')) THEN
         ICOM = ICOM + 1
         IF (ICOM .GT. JCOM) THEN
            CALL MESAGE ('ENTER LINE INTERVALS IN THE FOLLOWING '//
     &         'FORMAT:')
            CALL MESAGE ('[ LINE NO. (OR NEG SIDE NO.), INTERVALS ]')
            CALL MESAGE ('HIT RETURN TO END INPUT')
         END IF
  350    CONTINUE
         IF (ICOM .GT. JCOM) THEN
            CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN, CIN, IIN,
     &         RIN)
            ICOM = 1
         END IF
         CALL GETI12(MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2, IFOUND)
         III(1) = I1
         IF (IFOUND .GT. 0) THEN
            CALL ININTR(ML, MS, 1, I2, III, N(19), N(20), NINT, NLPS,
     &         IFLINE, ILLIST, LINKL, LINKS, ADDLNK)
            GO TO 350
         END IF

C  ENTER LINE FACTORS

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'F') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'f')) THEN
         ICOM = ICOM + 1
         IF (ICOM .GT. JCOM) THEN
            CALL MESAGE ('ENTER LINE FACTORS IN THE FOLLOWING FORMAT:')
            CALL MESAGE ('[ LINE NO. (OR NEG. SIDE NO., ) FACTOR ]')
            CALL MESAGE ('HIT RETURN TO END INPUT')
         END IF
  360    CONTINUE
         IF (ICOM .GT. JCOM) THEN
            CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN, CIN, IIN,
     &         RIN)
            ICOM = 1
         END IF
         CALL GETI1R(MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN, I1, R1,
     &      IFOUND)
         III(1) = I1
         IF (IFOUND .GT. 0) THEN
            CALL INFACT(ML, MS, 1, R1, III(1), N(19), N(20), FACTOR,
     &         NLPS, IFLINE, ILLIST, LINKL, LINKS, ADDLNK)
            GO TO 360
         END IF

C  ENTER MATERIAL NUMBERS

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'M') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'm')) THEN
         ICOM = ICOM + 1
         CALL MESAGE ('ENTER REGION MATERIALS IN THE FOLLOWING FORMAT:')
         CALL MESAGE ('[ REGION NO., MATERIAL NO. ]')
         CALL MESAGE ('HIT RETURN TO END INPUT')
         ICOM = JCOM + 1
  370    CONTINUE
         CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, IFOUND, KIN, CIN, IIN,
     &      RIN)
         IF (IFOUND .GT. 0) THEN
            JJ = IIN(1)
            CALL LTSORT(MR, LINKR, JJ, IPNTR, ADDLNK)
            IF ((JJ .GT. N(22)) .OR. (IPNTR .LE. 0)) THEN
               WRITE(*, 10170)JJ
            ELSE
               IMAT(IPNTR) = IIN(2)
            END IF
            GO TO 370
         END IF

C  FLAG THE THREE-NODE ELEMENT OPTION

      ELSE IF ((CIN(ICOM)(1:2) .EQ. 'TH') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'th')) THEN
         ICOM = ICOM + 1
         IF (THREE) THEN
            THREE = .FALSE.
            CALL MESAGE ('THREE NODE ELEMENT GENERATION - OFF')
         ELSE
            THREE = .TRUE.
            CALL MESAGE ('THREE NODE ELEMENT GENERATION - ON')
         END IF

C  FLAG THE EIGHT-NODE ELEMENT OPTION

      ELSE IF ((CIN(ICOM)(1:2) .EQ. 'EI') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'ei')) THEN
         ICOM = ICOM + 1
         IF (EIGHT) THEN
            EIGHT = .FALSE.
            CALL MESAGE ('EIGHT NODE ELEMENT GENERATION - OFF')
         ELSE
            EIGHT = .TRUE.
            NINE = .FALSE.
            CALL MESAGE ('EIGHT NODE ELEMENT GENERATION - ON')
         END IF

C  FLAG THE NINE-NODE ELEMENT OPTION

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'N') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'n')) THEN
         ICOM = ICOM + 1
         IF (NINE) THEN
            NINE = .FALSE.
            CALL MESAGE ('NINE NODE ELEMENT GENERATION - OFF')
         ELSE
            NINE = .TRUE.
            EIGHT = .FALSE.
            CALL MESAGE ('NINE NODE ELEMENT GENERATION - ON')
         END IF

C  EXIT OPTION - EXITS FASTQ

      ELSE IF ((CIN(ICOM)(1:2) .EQ. 'EX') .OR.
     &   (CIN(ICOM)(1:2) .EQ. 'ex')) THEN
         ICOM = ICOM + 1
         CALL STRLNG (CIN(ICOM), LEN)
         IF (((LEN .GT. 1) .AND. (CIN(ICOM)(2:2) .NE. 'X')) .OR.
     &      ((LEN .GT. 1) .AND. (CIN(ICOM)(2:2) .NE. 'x'))) THEN
            CALL HELP_FQ(9)
         ELSE
            CALL FEXIT(WROTE, MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN,
     &         TIME1, BATCH, VERSN)
         ENDIF
         GO TO 110

C  LINK ALL NEW DATA AS NEEDED, AND RETURN FROM THE KEYIN OPTION

      ELSE IF (CIN(ICOM)(1:1) .EQ. ' ') THEN
         ICOM = ICOM + 1

C  LINK UP THE POINTS AND LINES TO THEIR ASSOCIATED FLAGS

         SIDEOK = .FALSE.
         CALL LINKBC(MP, MS, NOLD(11) + 1, N(11), N(1), N(25), N(11),
     &      N(12), N(20), IPBF, IFPB, NPPF, LISTPB, NLPS, IFLINE,
     &      ILLIST, IPBOUN, LINKPB, IWTPBF, LINKP, LINKS, SIDEOK,
     &      NOROOM)
         IF (NOROOM) GO TO 400
         SIDEOK = .TRUE.
         CALL LINKBC(ML, MS, NOLD(13) + 1, N(13), N(2), N(26), N(13),
     &      N(14), N(20), ILBF, IFLB, NLPF, LISTLB, NLPS, IFLINE,
     &      ILLIST, ILBOUN, LINKLB, IWTLBF, LINKL, LINKS, SIDEOK,
     &      NOROOM)
         IF (NOROOM) GO TO 400
         CALL LINKBC(ML, MS, NOLD(15) + 1, N(15), N(2), N(27), N(15),
     &      N(16), N(20), ISBF, IFSB, NSPF, LISTSB, NLPS, IFLINE,
     &      ILLIST, ISBOUN, LINKSB, IWTSBF, LINKL, LINKS, SIDEOK,
     &      NOROOM)
         IF (NOROOM) GO TO 400

C  IF NO BODY CARDS HAVE BEEN READ, ASSUME THE BODY IS ALL THE
C  REGIONS AND ALL THE BAR SETS

         IF (N(9) .EQ. NOLD(9)) THEN
            ADDOLD = .TRUE.
            IFOUND = 1
            DO 380 I = NOLD(5) + 1, N(5)
               CALL LTSORT(MS, LINKB, IBARST(I), IPNTR, ADDLNK)
               IF (IPNTR .EQ. I) THEN
                  III(1) = -IBARST(I)
                  CALL INBODY(MR, N(9), III, IFOUND, IRPB, ADDOLD,
     &               NOROOM)
                  IF (NOROOM) GO TO 400
               END IF
  380       CONTINUE
            DO 390 I = NOLD(7) + 1, N(7)
               CALL LTSORT(MR, LINKR, IREGN(I), IPNTR, ADDLNK)
               IF (IPNTR .EQ. I) THEN
                  CALL INBODY(MR, N(9), IREGN(I), IFOUND, IRPB, ADDOLD,
     &               NOROOM)
                  IF (NOROOM) GO TO 400
               END IF
  390       CONTINUE
         END IF
         RETURN

      ELSE
         ICOM = ICOM + 1
         CALL HELP_FQ(9)
      END IF
      GO TO 110

C  MORE ROOM IN DIMENSIONS NEEDED

  400 CONTINUE
      CALL MESAGE (' ')
      CALL MESAGE ('DIMENSIONS MUST BE INCREASED - PLEASE WAIT')
      DO 410 I = 1, 29
         N(I) = NOLD(I)
  410 CONTINUE
      NOROOM = .TRUE.
      RETURN

10000 FORMAT(' A POINT NO. OF:', I7, ' IS NOT ALLOWED', /,
     &   ' THIS POINT WILL NOT BE INPUT INTO DATABASE')
10010 FORMAT(' A LINE NO. OF:', I7, ' IS NOT ALLOWED', /,
     &   ' THIS LINE WILL NOT BE INPUT INTO DATABASE')
10020 FORMAT(' FOR LINE NO.:', I5, ' NOT ENOUGH INFORMATION IS '//
     &   'SUPPLIED',/, ' THIS LINE WILL NOT BE INPUT INTO DATABASE')
10030 FORMAT(' REGION NO:', I5, ' IS NOT IN THE DATABASE', /,
     &   ' THUS NO SIZE CAN BE ENTERED')
10040 FORMAT(' A SIDE NO. OF:', I7, ' IS NOT ALLOWED', /,
     &   ' THIS SIDE WILL NOT BE INPUT INTO DATABASE')
10050 FORMAT(' A BAR SET NO. OF:', I7, ' IS NOT ALLOWED', /,
     &   ' THIS BAR SET WILL NOT BE INPUT INTO DATABASE')
10060 FORMAT(' FOR BAR SET NO.:', I5,
     &   /, ' NOT ENOUGH INFORMATION IS SUPPLIED',
     &   /, ' THIS BAR SET WILL NOT BE INPUT INTO DATABASE')
10070 FORMAT(' FOR SIDE NO.:', I5, ' NOT ENOUGH INFORMATION IS '//
     &   'SUPPLIED',/, ' THIS SIDE WILL NOT BE INPUT INTO DATABASE')
10080 FORMAT(' REGION NO:', I5, ' IS NOT IN THE DATABASE', /,
     &   ' THUS NO HOLE CAN BE ENTERED')
10090 FORMAT(' REGION:', I5, ' HAS LESS THAN ONE HOLE', /,
     &   ' THE HOLE FOR THIS REGION WILL NOT BE INPUT INTO DATABASE')
10100 FORMAT(' WEIGHTING OF ', A, ' FLAG:', I5, ' NOT POSSIBLE - '//
     &   'FLAG NOT FOUND')
10110 FORMAT(' A REGION NO. OF:', I7, ' IS NOT ALLOWED', /,
     &   ' THIS REGION WILL NOT BE INPUT INTO DATABASE')
10120 FORMAT(' FOR REGION NO.:', I5, ' NOT ENOUGH INFORMATION IS '//
     &   'SUPPLIED', /, ' THIS REGION WILL NOT BE INPUT INTO DATABASE')
10130 FORMAT(' A GROUP NO. OF:', I7, ' IS NOT ALLOWED', /,
     &   ' THIS GROUP WILL NOT BE INPUT INTO DATABASE')
10140 FORMAT(' FOR GROUP NO.:', I5, ' NOT ENOUGH INFORMATION IS '//
     &   'SUPPLIED', /, ' THIS GROUP WILL NOT BE INPUT INTO DATABASE')
10150 FORMAT(' A REGION NO. OF:', I7, ' IS NOT ALLOWED', /,
     &   ' THE SCHEME FOR THIS REGION WILL NOT BE INPUT INTO DATABASE')
10160 FORMAT(' FOR THE SCHEME FOR REGION NO.:', I5,
     &   ' NOT ENOUGH INFORMATION IS SUPPLIED',
     &   /, ' THIS SCHEME WILL NOT BE INPUT INTO DATABASE')
10170 FORMAT(' REGION NO:', I5, ' IS NOT IN THE DATABASE', /,
     &   ' THUS NO MATERIAL NUMBER CAN BE ENTERED')

      END
