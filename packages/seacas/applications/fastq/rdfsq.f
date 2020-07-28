C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE RDFSQ (MP, ML, MS, MR, MSNAP, MSC, MA, IUNIT, IDUMP, N,
     &   IPOINT, COOR, IPBOUN, ILINE, LTYPE, NINT, FACTOR, LCON, ILBOUN,
     &   ISBOUN, ISIDE, NLPS, IFLINE, ILLIST, IBARST, JMAT, JCENT, NLPB,
     &   JFLINE, JLLIST, IREGN, IMAT, NSPR, IFSIDE, ISLIST, IRPB, IPBF,
     &   NPPF, IFPB, LISTPB, ILBF, NLPF, IFLB, LISTLB, ISBF, NSPF, IFSB,
     &   LISTSB, ATTRIB, LINKP, LINKL, LINKS, LINKB, LINKR, LINKM,
     &   LINKSC, LINKPB, LINKLB, LINKSB, IHOLDP, IHOLDL, IHOLDS, IHOLDB,
     &   IHOLDR, IHOLDM, IHOLD1, IHOLD2, IHOLD3, IWTPBF, IWTLBF, IWTSBF,
     &   RSIZE, IFHOLE, NHPR, IHLIST, IRGFLG, ISCHM, SCHEME, NUMBER,
     &   DEFSCH, DEFSIZ, TITLE, OPTIM, MERGE, THREE, EIGHT, NINE,
     &   SNAP, SNAPDX, NSNAP, RATIO, NOROOM, EXODUSII)
C***********************************************************************

C  SUBROUTINE RDFSQ =  READS AND/OR MERGES FASTQ CARD FILE(S)

C***********************************************************************

C  SUBROUTINE CALLED BY:
C     FASTQ  =  A PROGRAM TO QUICKLY PREPARE FASTQ INPUT

C***********************************************************************

C  VARIABLES USED:
C     TITLE   =  MESH TITLE
C     IHOLDP  =  AN ARRAY TO HOLD THE POINTS DURING RENUMBERING
C     IHOLDL  =  AN ARRAY TO HOLD THE LINES DURING RENUMBERING
C     IHOLDS  =  AN ARRAY TO HOLD THE SIDES DURING RENUMBERING
C     IHOLDR  =  AN ARRAY TO HOLD THE REGIONS DURING RENUMBERING
C     DUMB     =  DUMMY VARIABLE WHERE THE DATA IS READ IN
C     OPTIM   =  .TRUE. IF THE MESH IS TO BE OPTIMIZED
C     MERGE   =  .TRUE. IF THE DATA IS TO BE MERGED WITH EXISTING DATA
C     NOROOM  =  .TRUE. IF THE AMOUNT OF DATA EXCEEDS DIMENSIONED LIMITS
C     NODATA  =  .TRUE. IF NO DATA HAS BEEN READ FROM THE FILE

C***********************************************************************

      PARAMETER (NIN = 1000)

      DIMENSION IPOINT(MP), COOR(2, MP), IPBOUN(MP)
      DIMENSION ILINE(ML), LTYPE(ML), NINT(ML), FACTOR(ML), LCON(3, ML)
      DIMENSION ILBOUN(ML), ISBOUN(ML)
      DIMENSION ISIDE(MS), NLPS(MS), IFLINE(MS), ILLIST(MS*3)
      DIMENSION IBARST(MS), JMAT(MS), JCENT(MS), NLPB(MS), JFLINE(MS)
      DIMENSION JLLIST(3*MS)
      DIMENSION IREGN(MR), IMAT(MR), NSPR(MR), IFSIDE(MR), ISLIST(MR*4)
      DIMENSION IRPB(MR), ISCHM(MSC), SCHEME(MSC), RSIZE(MR)
      DIMENSION IPBF(MP), NPPF(MP), IFPB(MP), LISTPB(2, MP)
      DIMENSION ILBF(ML), NLPF(ML), IFLB(ML), LISTLB(2, ML)
      DIMENSION ISBF(ML), NSPF(ML), IFSB(ML), LISTSB(2, ML)
      DIMENSION IWTPBF(3, MP), IWTLBF(3, ML), IWTSBF(3, ML)
      DIMENSION ATTRIB((MR + MS)*MA)
      DIMENSION LINKP(2, MP), LINKL(2, ML), LINKS(2, MS), LINKB(2, MS)
      DIMENSION LINKR(2, MR), LINKM(2, (MS + MR)), LINKSC(2, MR)
      DIMENSION LINKPB(2, MP), LINKLB(2, ML), LINKSB(2, ML)
      DIMENSION IFHOLE(MR), NHPR(MR), IHLIST(MR*2), IRGFLG(MR)
      DIMENSION NUMBER(MSC)
      DIMENSION IHOLDP(2, MP), IHOLDL(2, ML), IHOLDS(2, MS)
      DIMENSION IHOLDB(2, MS), IHOLDR(2, MR), IHOLDM(2, (MS + MR))
      DIMENSION IHOLD1(2, MP), IHOLD2(2, ML), IHOLD3(2, ML)
      DIMENSION N(29), NOLD(29), III(1)
      DIMENSION KIN(NIN), CIN(NIN), IIN(NIN), RIN(NIN)
      DIMENSION SNAPDX(2, MSNAP), NSNAP(2)

      CHARACTER*72 SCHEME, DEFSCH
      CHARACTER*72 TITLE, HOLD, NUMBER*80, CIN*72

      LOGICAL OPTIM, MERGE, NEWNUM, NOROOM, ADDOLD, DOLINK, ERR
      LOGICAL NODATA, ADDLNK, THREE, EIGHT, NINE, SIDEOK, SNAP
      LOGICAL EXODUSII

C  SET UP THE INITIALIZATION OF VARIABLES

      DO 100 I = 1, 29
         NOLD(I) = N(I)
  100 CONTINUE
      LHOLD  = 0
      NHOLDP = 0
      NHOLDL = 0
      NHOLDS = 0
      NHOLDB = 0
      NHOLDR = 0
      NHOLDM = 0
      NHOLD1 = 0
      NHOLD2 = 0
      NHOLD3 = 0
      THREE = .FALSE.
      EIGHT = .FALSE.
      NINE = .FALSE.
      OPTIM = .FALSE.
      NOROOM = .FALSE.
      ADDOLD = .TRUE.
      NODATA = .TRUE.
      ADDLNK = .TRUE.
      DEFSCH = 'M'

C  READ THE INPUT CARDS AND SORT AS NEEDED

      DO 130 I = 1, MP + 2*ML + MS + 2*MR
         CALL FREFLD (IUNIT, IDUMP, ' ', NIN, IOSTAT, IFOUND,
     &      KIN, CIN, IIN, RIN)

C  CHECK FOR THE END OF THE FILE OR FOR AN ERROR

         IF (IOSTAT .LT. 0) GO TO 140

C  INPUT THE TITLE

         IF (CIN(1)(1:5) .EQ. 'TITLE') THEN
            NODATA = .FALSE.
            IF (MERGE) THEN
               CALL STRLNG (TITLE, LEN)
               IF (LEN .LE. 70) THEN
                  LHOLD = LEN
                  TITLE(LEN + 1:LEN + 1) = ':'
                  CALL GETINP (IUNIT, IDUMP, ' ', HOLD, IOSTAT)
                  CALL STRCUT (HOLD)
                  LEND = 71 - LEN
                  TITLE(LEN + 2:72) = HOLD(1:LEND)
               END IF
            ELSE
               CALL GETINP (IUNIT, IDUMP, ' ', TITLE, IOSTAT)
            END IF
            CALL STRCUT (TITLE)
            CALL STRLNG (TITLE, LEN)

C  INPUT A POINT INTO THE DATABASE

         ELSE IF (CIN(1)(1:5) .EQ. 'POINT') THEN
            NODATA = .FALSE.
            JJ = IIN(2)
            IF ((JJ .LE. 0) .OR. (JJ .GE. 10000)) THEN
               WRITE(*, 10000) JJ
               GO TO 120
            END IF
            CALL INPOIN (MP, N(1), N(18), JJ, RIN(3), RIN(4), NHOLDP,
     &         IHOLDP, IPOINT, COOR, IPBOUN, LINKP, MERGE, NOROOM)
            IF (NOROOM) GO TO 310

C  INPUT A LINE INTO THE DATABASE

         ELSE IF (CIN(1)(1:5) .EQ. 'LINE ')THEN
            NODATA = .FALSE.
            JJ = IIN(2)
            IF ((JJ .LE. 0) .OR. (JJ .GE. 10000)) THEN
               WRITE(*, 10010)JJ
               GO TO 120
            ELSE IF (IFOUND .LT. 5) THEN
               WRITE(*, 10020) JJ
               GO TO 120
            END IF
            IF ((CIN(3)(1:4) .EQ. 'STR ') .OR.
     &         (CIN(3)(1:4) .EQ. '    ')) THEN
               LTYP = 1
            ELSE IF (CIN(3)(1:4) .EQ. 'CIRC') THEN
               LTYP = 3
            ELSE IF (CIN(3)(1:4) .EQ. 'CIRM') THEN
               LTYP = 4
            ELSE IF (CIN(3)(1:4) .EQ. 'CIRR') THEN
               LTYP = 6
            ELSE IF (CIN(3)(1:4) .EQ. 'PARA') THEN
               LTYP = 5
            ELSE IF (CIN(3)(1:4) .EQ. 'CORN') THEN
               LTYP = 2
            ELSE IF (CIN(3)(1:4) .EQ. 'ELIP') THEN
               LTYP = 7
            ELSE IF (CIN(3)(1:4) .EQ. 'ELPR') THEN
               LTYP = 8
            ELSE
               LTYP = 1
               WRITE(*, 10030) CIN(3)(1:4), JJ
            END IF
            CALL INLINE (ML, N(2), N(19), JJ, LTYP, IIN(4), IIN(5),
     &         IIN(6), IIN(7), RIN(8), NHOLDL, IHOLDL, ILINE, LTYPE,
     &         NINT, FACTOR, LCON, ILBOUN, ISBOUN, LINKL, MERGE, NOROOM)
            IF (NOROOM) GO TO 310

C  INPUT A SIDE INTO THE DATABASE

         ELSE IF (CIN(1)(1:5) .EQ. 'SIDE ') THEN
            NODATA = .FALSE.
            JJ = IIN(2)
            IF ((JJ .LE. 0) .OR. (JJ .GE. 10000)) THEN
               WRITE(*, 10040) JJ
               GO TO 120
            ELSE IF (IFOUND .LT. 3) THEN
               WRITE(*, 10050) JJ
               GO TO 120
            END IF
            CALL INSIDE (MS, N(3), N(4), N(20), JJ, IIN(3), IFOUND - 2,
     &         ISIDE, NLPS, IFLINE, ILLIST, LINKS, NHOLDS, IHOLDS,
     &         MERGE, NOROOM)
            IF (NOROOM) GO TO 310

C  INPUT A BAR SET INTO THE DATABASE

         ELSE IF (CIN(1)(1:6) .EQ. 'BARSET') THEN
            NODATA = .FALSE.
            JJ = IIN(2)
            IF ((JJ .LE. 0) .OR. (JJ .GE. 10000)) THEN
               WRITE(*, 10040) JJ
               GO TO 120
            ELSE IF (IFOUND .LT. 5) THEN
               WRITE(*, 10050) JJ
               GO TO 120
            END IF
            CALL INBRST (MS, MR, N(5), N(6), N(21), N(23), JJ, IIN(3),
     &         IIN(4), IIN(5), IFOUND - 4, IBARST, JMAT, JCENT, NLPB,
     &         JFLINE, JLLIST, LINKB, LINKM, NHOLDM, IHOLDM, NHOLDB,
     &         IHOLDB, MERGE, NOROOM)
            IF (NOROOM) GO TO 310

C  INPUT A REGION INTO THE DATABASE

         ELSE IF (CIN(1)(1:6) .EQ. 'REGION') THEN
            NODATA = .FALSE.
            JJ = IIN(2)
            IF ((JJ .LE. 0) .OR. (JJ .GE. 10000)) THEN
               WRITE(*, 10060) JJ
               GO TO 120
            ELSE IF (IFOUND .LT. 4) THEN
               WRITE(*, 10080) JJ
               GO TO 120
            END IF
            CALL INREGN (MS, MR, N(7), N(8), N(22), N(23), JJ, IIN(3),
     &         IIN(4), IFOUND - 3, IREGN, IMAT, NSPR, IFSIDE, ISLIST,
     &         LINKR, LINKM, NHOLDR, IHOLDR, NHOLDM, IHOLDM, IRGFLG,
     &         MERGE, NOROOM)
            ADDLNK = .FALSE.
            CALL LTSORT (MR, LINKR, IIN(2), JJPNTR, ADDLNK)
            ADDLNK = .TRUE.
            RSIZE(JJPNTR) = 0.
            NHPR(JJPNTR) = 0
            IF (NOROOM) GO TO 310

C  INPUT A GROUP INTO THE DATABASE

         ELSE IF (CIN(1)(1:6) .EQ. 'GROUP ') THEN
            NODATA = .FALSE.
            JJ = IIN(2)
            IF ((JJ .LE. 0) .OR. (JJ .GE. 10000)) THEN
               WRITE(*, 10060) JJ
               GO TO 120
            ELSE IF (IFOUND .LT. 3) THEN
               WRITE(*, 10080) JJ
               GO TO 120
            END IF
            CALL INGRPN (MS, MR, N(7), N(8), N(22), JJ, IIN(3),
     &         IFOUND - 2, IREGN, NSPR, IFSIDE, ISLIST, LINKR, NHOLDR,
     &         IHOLDR, IRGFLG, MERGE, NOROOM)
            IF (NOROOM) GO TO 310

C  INPUT A REGION'S HOLES INTO THE DATABASE

         ELSE IF (CIN(1)(1:6) .EQ. 'HOLE  ') THEN
            NODATA = .FALSE.
            JJ = IIN(2)
            ADDLNK = .FALSE.
            CALL LTSORT (MR, LINKR, JJ, JJPNTR, ADDLNK)
            IF (MERGE .AND. JJPNTR .GT. 0) THEN
               CALL LTSORT (MR, IHOLDR, JJ, JJNEW, ADDLNK)
               IF (JJNEW .GT. 0) THEN
                  CALL LTSORT (MR, LINKR, JJNEW, JJPNTR, ADDLNK)
               END IF
            END IF
            ADDLNK = .TRUE.
            IF (JJPNTR .LE. 0) THEN
               WRITE(*, 10090) JJ
               GO TO 120
            ELSE IF (IFOUND .LT. 3) THEN
               WRITE(*, 10100) JJ
               GO TO 120
            END IF
            CALL INHOLE (MR, N(7), N(29), JJPNTR, IIN(3), IFOUND - 2,
     &         IFHOLE, NHPR, IHLIST, MERGE, NOROOM)
            IF (NOROOM) GO TO 310

C  INPUT A SCHEME INTO THE DATABASE

         ELSE IF (CIN(1)(1:6) .EQ. 'SCHEME') THEN
            NODATA = .FALSE.
            JJ = IIN(2)
            IF (JJ .GE. 10000) THEN
               WRITE(*, 10110) JJ
               GO TO 120
            ELSE IF (IFOUND .LT. 3) THEN
               IF (JJ .EQ. 0) THEN
                  WRITE(*, 10120)
               ELSE
                  WRITE(*, 10130) JJ
               END IF
               GO TO 120
            END IF
            DOLINK = .FALSE.
            NOLD10 = N(10)
            NOLD24 = N(24)
            CALL INSCHM (MR, MSC, N(10), N(24), JJ, CIN(3), ISCHM,
     &         SCHEME, LINKSC, DEFSCH, NOROOM, DOLINK)
            IF (NOROOM) THEN
               N(10) = NOLD10
               N(24) = NOLD24
               CALL MESAGE ('************************************')
               CALL MESAGE ('NOT ENOUGH ROOM FOR SCHEME CARD')
               CALL MESAGE ('NO DYNAMIC DIMENSIONING INCREASES')
               CALL MESAGE ('AVAILABLE FOR CHARACTER STRINGS')
               CALL MESAGE ('SCHEME CARD IS THUS IGNORED')
               CALL MESAGE ('************************************')
            END IF

C  INPUT INTERVALS FOR SIDES OR LINES

         ELSE IF (CIN(1)(1:3) .EQ. 'INT') THEN
            NODATA = .FALSE.
            JJ = IIN(3)
            IF (IFOUND .LT. 3) THEN
               IF (JJ .LT. 0) THEN
                  WRITE(*, 10140) -JJ
               ELSE
                  WRITE(*, 10150) JJ
               END IF
               GO TO 120
            END IF
            ADDLNK = .FALSE.
            CALL ININTR (ML, MS, IFOUND - 2, IIN(2), IIN(3), N(19),
     &         N(20), NINT, NLPS, IFLINE, ILLIST, LINKL, LINKS, ADDLNK)
            ADDLNK = .TRUE.

C  INPUT FACTORS FOR SIDES OR LINES

         ELSE IF (CIN(1)(1:3) .EQ. 'FAC') THEN
            NODATA = .FALSE.
            JJ = IIN(3)
            IF (IFOUND .LT. 3) THEN
               IF (JJ .LT. 0) THEN
                  WRITE(*, 10160) -JJ
               ELSE
                  WRITE(*, 10170) JJ
               END IF
               GO TO 120
            END IF
            ADDLNK = .FALSE.
            CALL INFACT (ML, MS, IFOUND - 2, RIN(2), IIN(3), N(19),
     &         N(20), FACTOR, NLPS, IFLINE, ILLIST, LINKL, LINKS,
     &         ADDLNK)
            ADDLNK = .TRUE.

C  INPUT A REGION'S INTERVAL SIZE INTO THE DATABASE

         ELSE IF (CIN(1)(1:5) .EQ. 'SIZE ') THEN
            NODATA = .FALSE.
            IF (IFOUND .LT. 3) THEN
               DEFSIZ = AMAX1(RIN(2), 0.)
            ELSE
               ADDLNK = .FALSE.
               DO 110 IRSZ = 3, IFOUND
                  CALL LTSORT (MR, LINKR, IIN(IRSZ), JJ, ADDLNK)
                  IF (JJ .GT. 0) THEN
                     RSIZE(JJ) = AMAX1(RIN(2), 0.)
                  ELSE
                     WRITE(*, 10070) IIN(IRSZ)
                  END IF
  110          CONTINUE
               ADDLNK = .TRUE.
            END IF

C  INPUT A BODY DEFINITION INTO THE DATABASE

         ELSE IF (CIN(1)(1:4) .EQ. 'BODY') THEN
            NODATA = .FALSE.
            IF (IFOUND .GT. 1) THEN
               CALL INBODY (MR, N(9), IIN(2), IFOUND - 1, IRPB, ADDOLD,
     &            NOROOM)
               IF (NOROOM) GO TO 310
            END IF

C  INPUT A POINT BOUNDARY INTO THE DATABASE

         ELSE IF (CIN(1)(1:6) .EQ. 'POINBC') THEN
            IF (IFOUND .LT. 3) THEN
               WRITE(*, 10200) IIN(2)
               GO TO 120
            END IF
            NODATA = .FALSE.
            CALL INBOUN (MP, IIN(2), IFOUND - 2, IIN(3), N(25), N(11),
     &         N(12), NOLD(11), MERGE, NOROOM, NEWNUM, NHOLD1, IHOLD1,
     &         IPBF, NPPF, IFPB, LISTPB, LINKPB, IWTPBF, JHOLD, ADDOLD)
            IF (NOROOM) GO TO 310
            IF (NEWNUM) THEN
               ADDLNK = .FALSE.
               CALL LTSORT (MP, IHOLD1, JHOLD, IPNTR, ADDLNK)
               WRITE(*, 10210)JHOLD, IPNTR
               ADDLNK = .TRUE.
            END IF

C  INPUT A LINE BOUNDARY INTO THE DATABASE

         ELSE IF ((CIN(1)(1:6) .EQ. 'NODEBC')  .OR.
     &      (CIN(1)(1:6) .EQ. 'LINEBC')) THEN
            IF (IFOUND .LT. 3) THEN
               WRITE(*, 10220) IIN(2)
               GO TO 120
            END IF
            NODATA = .FALSE.
            CALL INBOUN (ML, IIN(2), IFOUND - 2, IIN(3), N(26), N(13),
     &         N(14), NOLD(13), MERGE, NOROOM, NEWNUM, NHOLD2, IHOLD2,
     &         ILBF, NLPF, IFLB, LISTLB, LINKLB, IWTLBF, JHOLD, ADDOLD)
            IF (NOROOM) GO TO 310
            IF (NEWNUM) THEN
               ADDLNK = .FALSE.
               CALL LTSORT (ML, IHOLD2, JHOLD, IPNTR, ADDLNK)
               WRITE(*, 10230)JHOLD, IPNTR
               ADDLNK = .TRUE.
            END IF

C  INPUT A SIDE BOUNDARY INTO THE DATABASE

         ELSE IF ((CIN(1)(1:6) .EQ. 'ELEMBC')  .OR.
     &      (CIN(1)(1:6) .EQ. 'SIDEBC')) THEN
            IF (IFOUND .LT. 3) THEN
               WRITE(*, 10240) IIN(2)
               GO TO 120
            END IF
            NODATA = .FALSE.
            CALL INBOUN (ML, IIN(2), IFOUND - 2, IIN(3), N(27), N(15),
     &         N(16), NOLD(15), MERGE, NOROOM, NEWNUM, NHOLD3, IHOLD3,
     &         ISBF, NSPF, IFSB, LISTSB, LINKSB, IWTSBF, JHOLD, ADDOLD)
            IF (NOROOM) GO TO 310
            IF (NEWNUM) THEN
               ADDLNK = .FALSE.
               CALL LTSORT (ML, IHOLD3, JHOLD, IPNTR, ADDLNK)
               WRITE(*, 10250)JHOLD, IPNTR
               ADDLNK = .TRUE.
            END IF

C  INPUT A FLAG WEIGHTING INTO THE DATABASE

         ELSE IF (CIN(1)(1:6) .EQ. 'WEIGHT') THEN
            ADDLNK = .FALSE.
            NODATA = .FALSE.

C  GET THE FLAG TYPE

            IF (CIN(2)(1:1) .EQ. 'P') THEN
               CALL LTSORT (MP, LINKPB, IIN(3), JJ, ADDLNK)
               IF (JJ .GT. 0) THEN
                  IWTPBF(1, JJ) = IIN(4)
                  IWTPBF(2, JJ) = IIN(5)
                  IWTPBF(3, JJ) = 0
               ELSE
                  WRITE(*, 10260) 'POINT', IIN(3)
               END IF
            ELSE IF (CIN(2)(1:1) .EQ. 'L') THEN
               CALL LTSORT (ML, LINKLB, IIN(3), JJ, ADDLNK)
               IF (JJ .GT. 0) THEN
                  IWTLBF(1, JJ) = IIN(4)
                  IWTLBF(2, JJ) = IIN(5)
                  IWTLBF(3, JJ) = IIN(6)
               ELSE
                  WRITE(*, 10270) 'LINE', IIN(3)
               END IF
            ELSE IF (CIN(2)(1:1) .EQ. 'S') THEN
               CALL LTSORT (ML, LINKSB, IIN(3), JJ, ADDLNK)
               IF (JJ .GT. 0) THEN
                  IWTSBF(1, JJ) = IIN(4)
                  IWTSBF(2, JJ) = IIN(5)
                  IWTSBF(3, JJ) = IIN(6)
               ELSE
                  WRITE(*, 10270) 'SIDE', IIN(3)
               END IF
            ELSE

C  NO FLAG TYPE HAS BEEN RECOGNIZED

               WRITE(*, 10280) CIN(2)(1:5)
            END IF
            ADDLNK = .TRUE.

C  FLAG THE BANDWIDTH OPTIMIZATION ROUTINES ON, AND READ A RENUM CARD

         ELSE IF (CIN(1)(1:5) .EQ. 'RENUM') THEN
            NODATA = .FALSE.
            OPTIM = .TRUE.
            IF (IFOUND .GT. 2) THEN
               HOLD = CIN(2)
               IDUM = IFOUND - 2
               CALL INRENM (MSC, N(28), HOLD, RIN(3), IIN(3), IDUM,
     &            NUMBER, NOROOM)
               IF (NOROOM) THEN
                  CALL MESAGE ('************************************')
                  CALL MESAGE ('NOT ENOUGH ROOM FOR RENUMBERING CARD')
                  CALL MESAGE ('NO DYNAMIC DIMENSIONING INCREASES')
                  CALL MESAGE ('AVAILABLE FOR CHARACTER STRINGS')
                  CALL MESAGE ('RENUMBERING CARD IS THUS IGNORED')
                  CALL MESAGE ('************************************')
               END IF
            ELSE IF (IFOUND .EQ. 2) THEN
               CALL MESAGE ('RENUM CARD READ WITHOUT ANY DATA')
               CALL MESAGE ('DEFAULT RENUM PROCESSING WILL BE USED')
            END IF

C  Write Database in exodusII format

         ELSE IF (CIN(1)(1:2) .EQ. 'X2') THEN
            NODATA = .FALSE.
            EXODUSII = .TRUE.

C  Write Database in exodusI/genesis format

         ELSE IF (CIN(1)(1:2) .EQ. 'X1') THEN
            NODATA = .FALSE.
            EXODUSII = .FALSE.

C  FLAG THE GENERATION OF THREE NODE ELEMENTS

         ELSE IF (CIN(1)(1:5) .EQ. 'THREE') THEN
            NODATA = .FALSE.
            THREE = .TRUE.

C  FLAG THE GENERATION OF EIGHT NODE ELEMENTS

         ELSE IF (CIN(1)(1:5) .EQ. 'EIGHT') THEN
            NODATA = .FALSE.
            EIGHT = .TRUE.
            NINE = .FALSE.

C  FLAG THE GENERATION OF NINE NODE ELEMENTS

         ELSE IF (CIN(1)(1:4) .EQ. 'NINE') THEN
            NODATA = .FALSE.
            NINE = .TRUE.
            EIGHT = .FALSE.

C  INPUT SNAP-TO-GRID FLAG

         ELSE IF (CIN(1)(1:4) .EQ. 'SNAP') THEN
            NODATA = .FALSE.
            SNAP = CIN(2)(1:2) .EQ. 'ON'

C  INPUT X-GRID LINES

         ELSE IF (CIN(1)(1:4) .EQ. 'XGRI') THEN
            IF (IFOUND .LT. 2) THEN
               WRITE(*, 10290) 'XGRID'
               GO TO 120
            END IF
            NODATA = .FALSE.
            CALL INGRID (MSNAP, SNAPDX, NSNAP, 1, RIN(2), IFOUND - 1,
     &         ERR)

C  INPUT Y-GRID LINES

         ELSE IF (CIN(1)(1:4) .EQ. 'YGRI') THEN
            IF (IFOUND .LT. 2) THEN
               WRITE(*, 10290) 'YGRID'
               GO TO 120
            END IF
            NODATA = .FALSE.
            CALL INGRID (MSNAP, SNAPDX, NSNAP, 2, RIN(2), IFOUND - 1,
     &         ERR)

C  END OF DATA

         ELSE IF (CIN(1)(1:4) .EQ. 'EXIT') THEN
            NODATA = .FALSE.
            GO TO 150
         END IF
  120    CONTINUE
  130 CONTINUE

  140 CONTINUE
      CALL MESAGE ('FILE END ENCOUNTERED BEFORE -EXIT- CARD WAS FOUND')
      CALL MESAGE ('POSSIBLE PROBLEM IN FILE')

C  RENUMBER THE CARDS IF MERGING

  150 CONTINUE
      ADDLNK = .FALSE.
      IF (MERGE) THEN

C  RENUMBER THE POINTS CONTAINED IN THE LINE, AND POINT BOUNDARY CARDS

         DO 170 I = NOLD(2) + 1, N(2)
            DO 160 J = 1, 3
               CALL LTSORT (MP, IHOLDP, LCON(J, I), IPNTR, ADDLNK)
               IF ((LCON(J, I) .LE. NHOLDP) .AND. (IPNTR .GT. 0))
     &            LCON(J, I) = IPNTR
  160       CONTINUE
  170    CONTINUE
         DO 180 I = NOLD(12) + 1, N(12)
            CALL LTSORT (MP, IHOLDP, LISTPB(1, I), IPNTR, ADDLNK)
            IF ((LISTPB(1, I) .LE. NHOLDP) .AND. (IPNTR .GT. 0))
     &         LISTPB(1, I) = IPNTR
  180    CONTINUE

C  RENUMBER THE LINES CONTAINED IN THE SIDE, BAR SET, REGION,
C  LINE BOUNDARY, AND SIDE BOUNDARY CARDS

         DO 190 I =  NOLD(4) + 1, N(4)
            CALL LTSORT (ML, IHOLDL, ILLIST(I), IPNTR, ADDLNK)
            IF ((ILLIST(I) .LE. NHOLDL) .AND. (IPNTR .GT. 0))
     &         ILLIST(I) = IPNTR
  190    CONTINUE
         DO 200 I =  NOLD(6) + 1, N(6)
            CALL LTSORT (ML, IHOLDL, JLLIST(I), IPNTR, ADDLNK)
            IF ((JLLIST(I) .LE. NHOLDL) .AND. (IPNTR .GT. 0))
     &         JLLIST(I) = IPNTR
  200    CONTINUE
         DO 210 I = NOLD(8) + 1, N(8)
            IF (ISLIST(I) .LT. 0)THEN
               KK = ABS(ISLIST(I))
               CALL LTSORT (ML, IHOLDL, KK, IPNTR, ADDLNK)
               IF ((KK .LE. NHOLDL) .AND. (IPNTR .GT. 0))
     &            ISLIST(I) = -IPNTR
            END IF
  210    CONTINUE
         DO 220 I = NOLD(14) + 1, N(14)
            CALL LTSORT (ML, IHOLDL, LISTLB(1, I), IPNTR, ADDLNK)
            IF ((LISTLB(1, I) .LE. NHOLDL) .AND. (IPNTR .GT. 0))
     &         LISTLB(1, I) = IPNTR
  220    CONTINUE
         DO 230 I = NOLD(16) + 1, N(16)
            CALL LTSORT (ML, IHOLDL, LISTSB(1, I), IPNTR, ADDLNK)
            IF ((LISTSB(1, I) .LE. NHOLDL) .AND. (IPNTR .GT. 0))
     &         LISTSB(1, I) = IPNTR
  230    CONTINUE

C  RENUMBER THE SIDES CONTAINED IN THE REGION CARDS

         DO 240 I = NOLD(8) + 1, N(8)
            IF (ISLIST(I) .GT. 0) THEN
               CALL LTSORT (MS, IHOLDS, ISLIST(I), IPNTR, ADDLNK)
               IF ((ISLIST(I) .LE. NHOLDS) .AND. (IPNTR .GT. 0))
     &            ISLIST(I) = IPNTR
            END IF
  240    CONTINUE

C  RENUMBER THE REGIONS CONTAINED IN THE HOLE CARDS

         DO 250 I = NOLD(29) + 1, N(29)
            IF (IHLIST(I) .GT. 0) THEN
               CALL LTSORT (MR, IHOLDR, IHLIST(I), IPNTR, ADDLNK)
               IF ((IHLIST(I) .LE. NHOLDR) .AND. (IPNTR .GT. 0))
     &            IHLIST(I) = IPNTR
            END IF
  250    CONTINUE

C  RENUMBER THE BAR SETS AND REGIONS CONTAINED IN THE BODY,
C  AND THE REGIONS CONTAINED IN THE SCHEME CARDS

         DO 260 I = NOLD(9) + 1, N(9)
            IF (IRPB(I) .GT. 0) THEN
               CALL LTSORT (MR, IHOLDR, IRPB(I), IPNTR, ADDLNK)
               IF ((IRPB(I) .LE. NHOLDR) .AND. (IPNTR .GT. 0))
     &            IRPB(I) = IPNTR
            ELSE IF (IRPB(I) .LT. 0) THEN
               CALL LTSORT (MS, IHOLDB, IABS(IRPB(I)), IPNTR, ADDLNK)
               IF ((IABS(IRPB(I)) .LE. NHOLDB) .AND. (IPNTR .GT. 0))
     &            IRPB(I) = -IPNTR
            END IF
  260    CONTINUE
         DO 270 I = NOLD(10) + 1, N(10)
            CALL LTSORT (MR, IHOLDR, ISCHM(I), IPNTR, ADDLNK)
            IF ((ISCHM(I) .LE. NHOLDR) .AND. (IPNTR .GT. 0))
     &         ISCHM(I) = IPNTR
  270    CONTINUE
      END IF

C  LINK THE SCHEME CARDS

      ADDLNK = .TRUE.
      DO 280 I = NOLD(10) + 1, N(10)
         IF (ISCHM(I) .GT. N(24)) N(24) = ISCHM(I)
         CALL LTSORT (MR, LINKSC, ISCHM(I), I, ADDLNK)
  280 CONTINUE

C  LINK UP THE POINTS AND LINES TO THEIR ASSOCIATED FLAGS

      SIDEOK = .FALSE.
      CALL LINKBC (MP, MS, NOLD(11) + 1, N(11), N(1), N(25), N(11),
     &   N(12), N(20), IPBF, IFPB, NPPF, LISTPB, NLPS, IFLINE, ILLIST,
     &   IPBOUN, LINKPB, IWTPBF, LINKP, LINKS, SIDEOK, NOROOM)
      IF (NOROOM) GO TO 310
      SIDEOK = .TRUE.
      CALL LINKBC (ML, MS, NOLD(13) + 1, N(13), N(2), N(26), N(13),
     &   N(14), N(20), ILBF, IFLB, NLPF, LISTLB, NLPS, IFLINE, ILLIST,
     &   ILBOUN, LINKLB, IWTLBF, LINKL, LINKS, SIDEOK, NOROOM)
      IF (NOROOM) GO TO 310
      CALL LINKBC (ML, MS, NOLD(15) + 1, N(15), N(2), N(27), N(15),
     &   N(16), N(20), ISBF, IFSB, NSPF, LISTSB, NLPS, IFLINE, ILLIST,
     &   ISBOUN, LINKSB, IWTSBF, LINKL, LINKS, SIDEOK, NOROOM)
      IF (NOROOM) GO TO 310

C  IF NO BODY CARDS HAVE BEEN READ, ASSUME THE BODY IS ALL THE REGIONS
C  AND ALL THE BAR SETS

      ADDLNK = .FALSE.
      IF (N(9) .EQ. NOLD(9)) THEN
         IFOUND = 1
         DO 290 I = NOLD(5) + 1, N(5)
            CALL LTSORT (MS, LINKB, IBARST(I), IPNTR, ADDLNK)
            IF (IPNTR .EQ. I) THEN
               III(1) = -IBARST(I)
               CALL INBODY (MR, N(9), III(1), IFOUND, IRPB, ADDOLD,
     &            NOROOM)
               IF (NOROOM) GO TO 310
            END IF
  290    CONTINUE
         DO 300 I = NOLD(7) + 1, N(7)
            CALL LTSORT (MR, LINKR, IREGN(I), IPNTR, ADDLNK)
            IF (IPNTR .EQ. I) THEN
               CALL INBODY (MR, N(9), IREGN(I), IFOUND, IRPB, ADDOLD,
     &            NOROOM)
               IF (NOROOM) GO TO 310
            END IF
  300    CONTINUE
      END IF

C  SUCCESSFUL COMPLETION - WRITE SUMMARY OF SUCCESSFUL READS

      IF (NODATA) THEN
         CALL MESAGE (' ')
         CALL MESAGE (' *----------------------------------------- - *')
         CALL MESAGE (' NO FASTQ DATA HAS BEEN FOUND IN CURRENT FILE')
         CALL MESAGE (' *----------------------------------------- - *')
         RETURN
      END IF
      CALL MESAGE (' ')
      CALL MESAGE ('FILE SUCCESSFULLY READ')
      CALL STRLNG (TITLE, LEN)
      WRITE(*, 10300) TITLE(1:LEN)
      IF (MERGE) THEN
         IF (N(1) .GT. 0) WRITE(*, 10310) N(1) - NOLD(1), N(1)
         IF (N(2) .GT. 0) WRITE(*, 10320) N(2) - NOLD(2), N(2)
         IF (N(3) .GT. 0) WRITE(*, 10330) N(3) - NOLD(3), N(3)
         IF (N(5) .GT. 0) WRITE(*, 10340) N(5) - NOLD(5), N(5)
         IF (N(7) .GT. 0) WRITE(*, 10350) N(7) - NOLD(7), N(7)
         IF (N(10) .GT. 0) WRITE(*, 10360) N(10) - NOLD(10), N(10)
         IF (N(11) .GT. 0) WRITE(*, 10370) N(11) - NOLD(11), N(11)
         IF (N(13) .GT. 0) WRITE(*, 10380) N(13) - NOLD(13), N(13)
         IF (N(15) .GT. 0) WRITE(*, 10390) N(15) - NOLD(15), N(15)
         IF (N(28) .GT. 0) WRITE(*, 10400) N(28) - NOLD(28), N(28)
         IF (N(29) .GT. 0) WRITE(*, 10410) N(29) - NOLD(29), N(29)
      ELSE
         IF (N(1) .GT. 0) WRITE(*, 10420) N(1)
         IF (N(2) .GT. 0) WRITE(*, 10430) N(2)
         IF (N(3) .GT. 0) WRITE(*, 10440) N(3)
         IF (N(5) .GT. 0) WRITE(*, 10450) N(5)
         IF (N(7) .GT. 0) WRITE(*, 10460) N(7)
         IF (N(10) .GT. 0) WRITE(*, 10480) N(10)
         IF (N(11) .GT. 0) WRITE(*, 10490) N(11)
         IF (N(13) .GT. 0) WRITE(*, 10500) N(13)
         IF (N(15) .GT. 0) WRITE(*, 10510) N(15)
         IF (N(28) .GT. 0) WRITE(*, 10520) N(28)
         IF (N(29) .GT. 0) WRITE(*, 10470) N(29)
      END IF
      RETURN

C  MORE ROOM IN DIMENSIONS NEEDED

  310 CONTINUE
      CALL MESAGE (' ')
      CALL MESAGE ('DIMENSIONS MUST BE INCREASED - PLEASE WAIT')
      DO 320 I = 1, 29
         N(I) = NOLD(I)
  320 CONTINUE
      IF (LHOLD .EQ. 0) THEN
         TITLE = ' '
      ELSE
         TITLE(LHOLD + 1:) = ' '
      END IF
      NOROOM = .TRUE.

C  FIND OUT HOW MUCH ROOM IS NEEDED

      REWIND IUNIT
      NEWMP = 0
      NEWML = 0
      NEWMS = 0
      NEWMR = 0
  330 CONTINUE
      CALL FREFLD (IUNIT, IDUMP, ' ', NIN, IOSTAT, IFOUND, KIN, CIN,
     &   IIN, RIN)

C  CHECK FOR THE END OF THE FILE OR FOR AN ERROR

      IF (IOSTAT .LT. 0) GO TO 340

C  COUNT THE CARDS FOR NEEDED DIMENSIONING

      IF (CIN(1)(1:5) .EQ. 'POINT') THEN
         NEWMP = NEWMP + 1
      ELSE IF (CIN(1)(1:5) .EQ. 'LINE ') THEN
         NEWML = NEWML + 1
      ELSE IF (CIN(1)(1:5) .EQ. 'SIDE ') THEN
         NEWMS = NEWMS + 1
      ELSE IF (CIN(1)(1:6) .EQ. 'REGION') THEN
         NEWMR = NEWMR + 1
      END IF
      GO TO 330

C  GET THE LARGEST RATIO OF NEEDED/CURRENT

  340 CONTINUE
      RATIO = AMAX1(DBLE(NEWMP)/DBLE(MP), DBLE(NEWML)/DBLE(ML),
     &   DBLE(NEWMS)/DBLE(MS), DBLE(NEWMR)/DBLE(MR), 1.5000001)*1.1
      RETURN

10000 FORMAT (' A POINT NO. OF:', I7, ' IS NOT ALLOWED', /,
     &   ' THIS POINT WILL NOT BE INPUT INTO DATABASE')
10010 FORMAT (' A LINE NO. OF:', I7, ' IS NOT ALLOWED', /,
     &   ' THIS LINE WILL NOT BE INPUT INTO DATABASE')
10020 FORMAT (' FOR LINE NO.:', I5, ' NOT ENOUGH INFORMATION IS ',
     &   'SUPPLIED',/, ' THIS LINE WILL NOT BE INPUT INTO DATABASE')
10030 FORMAT (' LINE TYPE "', A4, '" FOR LINE:', I5, ' IS ILLEGAL.', /,
     &   ' LINE TYPE CHANGED TO "STRAIGHT"')
10040 FORMAT (' A SIDE NO. OF:', I7, ' IS NOT ALLOWED', /,
     &   ' THIS SIDE WILL NOT BE INPUT INTO DATABASE')
10050 FORMAT (' FOR SIDE NO.:', I5, ' NOT ENOUGH INFORMATION IS ',
     &   'SUPPLIED',/, ' THIS SIDE WILL NOT BE INPUT INTO DATABASE')
10060 FORMAT (' A REGION NO. OF:', I7, ' IS NOT ALLOWED', /,
     &   ' THIS REGION WILL NOT BE INPUT INTO DATABASE')
10070 FORMAT (' REGION NO:', I5, ' IS NOT IN THE DATABASE', /,
     &   ' THUS NO SIZE CAN BE ENTERED')
10080 FORMAT (' FOR REGION NO.:', I5, ' NOT ENOUGH INFORMATION IS ',
     &   'SUPPLIED'/, ' THIS REGION WILL NOT BE INPUT INTO DATABASE')
10090 FORMAT (' A REGION NO. OF:', I7, ' IS NOT IN THE DATABASE', /,
     &   ' THE HOLES FOR THIS REGION WILL NOT BE INPUT INTO DATABASE')
10100 FORMAT (' FOR HOLE REGION NO.:', I5, ' NOT ENOUGH INFORMATION ',
     &   'IS SUPPLIED'/,' THE HOLES FOR THIS REGION WILL NOT BE INPUT ',
     &   'INTO DATABASE')
10110 FORMAT (' A REGION NO. OF:', I7, ' IS NOT ALLOWED', /,
     &   ' THE SCHEME FOR THIS REGION WILL NOT BE INPUT INTO DATABASE')
10120 FORMAT (' THE DEFAULT SCHEME HAS NOT BEEN SPECIFIED ADEQUATELY',
     &   /, ' THIS DEFAULT WILL NOT BE INPUT INTO THE DATABASE')
10130 FORMAT (' FOR THE SCHEME FOR REGION NO.:', I5,
     &   ' NOT ENOUGH INFORMATION IS SUPPLIED',
     &   /, ' THIS SCHEME WILL NOT BE INPUT INTO DATABASE')
10140 FORMAT (' FOR INTERVALS TO BE INPUT ON LINES IN SIDE: ', I5, /,
     &   ' NOT ENOUGH INFORMATION IS SUPPLIED', /,
     &   ' THESE INTERVALS WILL NOT BE ASSIGNED')
10150 FORMAT (' FOR INTERVALS TO BE INPUT ON LINE: ', I5, /,
     &   ' NOT ENOUGH INFORMATION IS SUPPLIED', /,
     &   ' THIS INTERVAL WILL NOT BE ASSIGNED')
10160 FORMAT (' FOR FACTORS TO BE INPUT ON LINES IN SIDE: ', I5, /,
     &   ' NOT ENOUGH INFORMATION IS SUPPLIED', /,
     &   ' THESE INTERVALS WILL NOT BE ASSIGNED')
10170 FORMAT (' FOR FACTORS TO BE INPUT ON LINE: ', I5, /,
     &   ' NOT ENOUGH INFORMATION IS SUPPLIED', /,
     &   ' THIS INTERVAL WILL NOT BE ASSIGNED')
10200 FORMAT (' FOR POINBC NO.:', I5, ' NOT ENOUGH INFORMATION IS ',
     &   'SUPPLIED',/, ' THIS POINBC WILL NOT BE INPUT INTO DATABASE')
10210 FORMAT (' OLD POINBC NO:', I5, ' TO NEW POINBC NO:', I5)
10220 FORMAT (' FOR NODEBC NO.:', I5, ' NOT ENOUGH INFORMATION IS ',
     &   'SUPPLIED',/, ' THIS NODEBC WILL NOT BE INPUT INTO DATABASE')
10230 FORMAT (' OLD NODEBC NO:', I5, ' TO NEW NODEBC NO:', I5)
10240 FORMAT (' FOR ELEMBC NO.:', I5, ' NOT ENOUGH INFORMATION IS ',
     &   'SUPPLIED',/, ' THIS ELEMBC WILL NOT BE INPUT INTO DATABASE')
10250 FORMAT (' OLD ELEMBC NO:', I5, ' TO NEW ELEMBC NO:', I5)
10260 FORMAT (' WEIGHTING OF ', A5, ' FLAG:', I5, ' NOT POSSIBLE', /,
     &   'FLAG NOT FOUND')
10270 FORMAT (' WEIGHTING OF ', A4, ' FLAG:', I5, ' NOT POSSIBLE', /,
     &   'FLAG NOT FOUND')
10280 FORMAT (' WEIGHTING TYPE OF ', A5, ' CANNOT BE ENTERED')
10290 FORMAT (' NOT ENOUGH INFORMATION SUPPLIED WITH CARD:  ', A)
10300 FORMAT (' TITLE: ', A)
10310 FORMAT (' ', I5, ' NEW POINTS READ   -    TOTAL POINTS:', I5)
10320 FORMAT (' ', I5, ' NEW LINES READ    -     TOTAL LINES:', I5)
10330 FORMAT (' ', I5, ' NEW SIDES READ    -     TOTAL SIDES:', I5)
10340 FORMAT (' ', I5, ' NEW BAR SETS READ -  TOTAL BAR SETS:', I5)
10350 FORMAT (' ', I5, ' NEW REGIONS READ  -   TOTAL REGIONS:', I5)
10360 FORMAT (' ', I5, ' NEW SCHEMES READ  -   TOTAL SCHEMES:', I5)
10370 FORMAT (' ', I5, ' NEW POINBCS READ  -   TOTAL POINBCS:', I5)
10380 FORMAT (' ', I5, ' NEW NODEBCS READ  -   TOTAL NODEBCS:', I5)
10390 FORMAT (' ', I5, ' NEW ELEMBCS READ  -   TOTAL ELEMBCS:', I5)
10400 FORMAT (' ', I5, ' NEW RENUMS READ   -    TOTAL RENUMS:', I5)
10410 FORMAT (' ', I5, ' NEW HOLES READ    -     TOTAL HOLES:', I5)
10420 FORMAT ('   NUMBER OF POINTS READ:', I5)
10430 FORMAT ('    NUMBER OF LINES READ:', I5)
10440 FORMAT ('    NUMBER OF SIDES READ:', I5)
10450 FORMAT (' NUMBER OF BAR SETS READ:', I5)
10460 FORMAT ('  NUMBER OF REGIONS READ:', I5)
10470 FORMAT ('    NUMBER OF HOLES READ:', I5)
10480 FORMAT ('  NUMBER OF SCHEMES READ:', I5)
10490 FORMAT ('  NUMBER OF POINBCS READ:', I5)
10500 FORMAT ('  NUMBER OF NODEBCS READ:', I5)
10510 FORMAT ('  NUMBER OF ELEMBCS READ:', I5)
10520 FORMAT ('   NUMBER OF RENUMS READ:', I5)

      END
