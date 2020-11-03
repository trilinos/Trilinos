C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE MESH (A, IA, MP, ML, MS, MR, MSC, MA, MCOM, ICOM, JCOM,
     &   CIN, RIN, IIN, KIN, IDUMP, N, IPOINT, COOR, IPBOUN, ILINE,
     &   LTYPE, NINT, FACTOR, LCON, ILBOUN, ISBOUN, ISIDE, NLPS, IFLINE,
     &   ILLIST, IBARST, JMAT, JCENT, NLPB, JFLINE, JLLIST, IREGN, IMAT,
     &   NSPR, IFSIDE, ISLIST, IRPB, IPBF, NPPF, IFPB, LISTPB, ILBF,
     &   NLPF, IFLB, LISTLB, ISBF, NSPF, IFSB, LISTSB, LINKP, LINKL,
     &   LINKS, LINKB, LINKR, LINKM, LINKSC, LINKPB, LINKLB, LINKSB,
     &   IWTPBF, IWTLBF, IWTSBF, RSIZE, IFHOLE, NHPR, IHLIST, IRGFLG,
     &   ISCHM, SCHEME, NUMBER, DEFSCH, DEFSIZ, TITLE, OPTIM, IDEV,
     &   ALPHA, DEV1, THREE, EIGHT, NINE, BATCH, VAXVMS, VERSN, AXIS,
     &   AREACG, LABN, LABE, LABO, LABNB, LABSB, LABM, LABW, WROTE,
     &   TIME1, HARDPL, EXODUSII)
C***********************************************************************

C  SUBROUTINE MESH = PROCESSES THE MESH AND THEN GRAPHICALLY DISPLAYS
C                    THE MESH

C***********************************************************************

C  SUBROUTINE CALLED BY:
C     FASTQ  = A PROGRAM TO QUICKLY PREPARE QMESH INPUT

C***********************************************************************

C  SUBROUTINES CALLED:
C     QMESH  = GENERATES THE MESH
C     RENUM  = ASSIGNS NODAL AND ELEMENT NUMBERS TO THE MESH

C***********************************************************************

C  VARIABLES USED:
C     IANS   = LOGICAL RESPONSE FROM YES-NO QUESTION
C     TITLE  = MESH TITLE
C     TALL   = CHARACTER HEIGHT SETTING PARAMETER
C     STEP   = .TRUE. IF GENERATION TO BE STEPPED THROUGH INTERACTIVELY
C     DONE   = .TRUE. IF MESH HAS BEEN PROCESSED

C***********************************************************************

      PARAMETER (MAXNAM = 40, MLINK = 55)

      include 'exodusII.inc'

C  DIMENSIONS FOR MESH DEFINING ENTITIES  (I.E. POINTS,  LINES,  ETC.)

      DIMENSION A(1), IA(1)
      DIMENSION IPOINT(MP), COOR(2, MP), IPBOUN(MP)
      DIMENSION ILINE(ML), LTYPE(ML), NINT(ML), FACTOR(ML)
      DIMENSION LCON(3, ML)
      DIMENSION ILBOUN(ML), ISBOUN(ML)
      DIMENSION ISIDE(MS), NLPS(MS), IFLINE(MS), ILLIST(MS*3)
      DIMENSION IREGN(MR), IMAT(MR)
      DIMENSION NSPR(MR), IFSIDE(MR), ISLIST(MR*4), IRGFLG(MR)
      DIMENSION RSIZE(MR), IFHOLE(MR), NHPR(MR), IHLIST(MR*2)
      DIMENSION IBARST(MS), JMAT(MS), JCENT(MS)
      DIMENSION NLPB(MS), JFLINE(MS), JLLIST(MS*3)
      DIMENSION IRPB(MR), ISCHM(MSC), SCHEME(MSC)
      DIMENSION IPBF(MP), NPPF(MP), IFPB(MP)
      DIMENSION LISTPB(2, MP), IWTPBF(3, MP)
      DIMENSION ILBF(ML), NLPF(ML), IFLB(ML)
      DIMENSION LISTLB(2, ML), IWTLBF(3, ML)
      DIMENSION ISBF(ML), NSPF(ML), IFSB(ML)
      DIMENSION LISTSB(2, ML), IWTSBF(3, ML)
      DIMENSION LINKP(2, MP), LINKL(2, ML), LINKS(2, MS)
      DIMENSION LINKB(2, MS)
      DIMENSION LINKR(2, MR), LINKM(2, (MS + MR)), LINKSC(2, MR)
      DIMENSION LINKPB(2, MP), LINKLB(2, ML), LINKSB(2, ML)
      DIMENSION NUMBER(MSC)
      DIMENSION N(29), K(41), IDEV(2), III(1)
      DIMENSION KIN(MCOM), IIN(MCOM), RIN(MCOM)

      CHARACTER*72 SCHEME, DEFSCH, TITLE, DEV1*3, CIN(MCOM)
      CHARACTER*2048 FNAME
      CHARACTER*80 NUMBER, HOLD, VERSN*10

      LOGICAL OPTIM, STEP, ERR, ALPHA, THREE, EIGHT, NINE
      LOGICAL AXIS, AREACG, LABE, LABO, LABN, LABNB, LABSB, LABM, LABW
      LOGICAL ADDLNK, BATCH, VAXVMS, WROTE, HARDPL, LGROUP
      LOGICAL REMESH, LONG
      LOGICAL EXODUSII, ISBARS

      CHARACTER*8 CDUMH, CDUMS
      INTEGER CMPSIZ

      NPREGN = 0
      IZ = 0
      ADDLNK = .FALSE.
      REMESH = .FALSE.

C  ENTER THE MESH OPTION

  100 CONTINUE
      IF (ICOM.GT.JCOM) THEN
         CALL MESAGE (' ')
         CALL FREFLD (IZ, IZ, 'ENTER MESH OPTION: ', MCOM, IOSTAT, JCOM,
     &      KIN, CIN, IIN, RIN)
         ICOM = 1
      END IF

C  RETURN FROM MESHING AFTER DELETING THE MESH

      IF (CIN(ICOM)(1:1) .EQ. ' ') THEN
         ICOM = ICOM + 1
         IF (NPREGN.GT.0) THEN
            CALL MDDEL ('IPART')
            CALL MDDEL ('LSTNBC')
            CALL MDDEL ('LSTSBC')
            CALL MDDEL ('NNFLG')
            CALL MDDEL ('NNPTR')
            CALL MDDEL ('NNLEN')
            CALL MDDEL ('NSFLG')
            CALL MDDEL ('NSPTR')
            CALL MDDEL ('NSLEN')
            CALL MDDEL ('NVPTR')
            CALL MDDEL ('NVLEN')
            CALL MDDEL ('NSIDEN')
            CALL MDDEL ('NUID')
            CALL MDDEL ('XN')
            CALL MDDEL ('YN')
            CALL MDDEL ('NXK')
            CALL MDDEL ('MAT')
            CALL MDDEL ('KXN')
            CALL MDDEL ('LIST')
            CALL MDDEL ('LA')
            CALL MDDEL ('LB')
            CALL MDDEL ('CENTK')
            CALL MDDEL ('MATMAP')
            CALL MDDEL ('LISTN')
            CALL MDDEL ('WTNODE')
            CALL MDDEL ('WTSIDE')
            CALL MDDEL ('WTHOLD')
            CALL MDDEL ('IHERE')
            CALL MDDEL ('ILIST')
            CALL MDDEL ('XLIST')
            CALL MDDEL ('AMESUR')
            CALL MDDEL ('XNOLD')
            CALL MDDEL ('YNOLD')
            CALL MDDEL ('NXKOLD')
            CALL MDDEL ('MMPOLD')
            CALL MDDEL ('LINKEG')
            CALL MDDEL ('LISTEG')
            CALL MDDEL ('BMESUR')
            CALL MDDEL ('CMESUR')
            CALL MDDEL ('DMESUR')
            CALL MDSTAT (NERR, MUSED)
            IF (NERR.GT.0) THEN
               CALL MDEROR (6)
               STOP ' '
            END IF
            KKK = 0
         END IF
         RETURN

C  SPECIFY ExodusI (X1) or ExodusII (X2) database format
      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'X') .OR.
     &        (CIN(ICOM)(1:1) .EQ. 'x')) THEN
         if (cin(icom)(2:2) .eq. '2') then
            exodusII = .TRUE.
            call mesage ('Writing EXODUSII Format')
         else if (cin(icom)(2:2) .eq. '1') then
            exodusII = .FALSE.
            call mesage ('Writing EXODUSI/GENESIS Format')
         end if
         ICOM = ICOM + 1

C  ENTER THE MESH GRAPHICS OPTION

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'G') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'g')) THEN
         ICOM = ICOM + 1
         IF (KKK.LE.0) THEN
            CALL MESAGE (' ')
            CALL MESAGE ('****************************************')
            CALL MESAGE ('*    NO MESH HAS BEEN PROCESSED        *')
            CALL MESAGE ('*        NO PLOTTING POSSIBLE          *')
            CALL MESAGE ('****************************************')
            CALL MESAGE (' ')
         ELSE
            CALL GMESH (NPNODE, NPELEM, MXNFLG, MXSFLG, NPNBC, NPSBC,
     &         MAXKXN, MR, NPREGN, MCOM, ICOM, JCOM, CIN, RIN, IIN, KIN,
     &         NNN, KKK, NUMMAT, NNXK, IA(K(1)), IA(K(2)), IA(K(4)),
     &         IA(K(5)), IA(K(7)), IA(K(10)), IA(K(11)), IA(K(12)),
     &         IA(K(13)), A(K(14)), A(K(15)), IA(K(16)),
     &         IA(K(17)), IA(K(18)), IA(K(19)), A(K(22)),
     &         IA(K(23)), A(K(25)), A(K(26)), NBCNOD, NNLIST,
     &         NBCSID, NVLIST, TITLE, IDUMP, AXIS, AREACG, LABE, LABO,
     &         LABN, LABNB, LABSB, LABM, LABW, IDEV, ALPHA, DEV1, EIGHT,
     &         NINE, VAXVMS, VERSN, WROTE, TIME1, HARDPL, BATCH)
         END IF

C  TOGGLE OPTIMIZATION

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'O') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'o')) THEN
         ICOM = ICOM + 1
         IF (OPTIM) THEN
            OPTIM = .FALSE.
            CALL MESAGE (' ')
            CALL MESAGE ('BANDWIDTH OPTIMIZATION - DISABLED')
         ELSE
            OPTIM = .TRUE.
            CALL MESAGE (' ')
            CALL MESAGE ('BANDWIDTH OPTIMIZATION - ENABLED')
         END IF

C  TOGGLE THREE NODE BAR ELEMENTS

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'T') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 't')) THEN
         ICOM = ICOM + 1
         IF (THREE) THEN
            THREE = .FALSE.
            CALL MESAGE (' ')
            CALL MESAGE ('THREE NODE BAR GENERATION - DISABLED')
         ELSE
            THREE = .TRUE.
            CALL MESAGE (' ')
            CALL MESAGE ('THREE NODE BAR GENERATION - ENABLED')
         END IF
         KKK = 0

C  TOGGLE EIGHT NODE ELEMENTS

      ELSE IF ((CIN(ICOM)(1:2) .EQ. 'EI') .OR.
     &   (CIN(ICOM)(1:2) .EQ. 'ei')) THEN
         ICOM = ICOM + 1
         IF (EIGHT) THEN
            EIGHT = .FALSE.
            CALL MESAGE (' ')
            CALL MESAGE ('EIGHT NODE QUAD GENERATION - DISABLED')
         ELSE
            EIGHT = .TRUE.
            NINE = .FALSE.
            CALL MESAGE (' ')
            CALL MESAGE ('EIGHT NODE QUAD GENERATION - ENABLED')
         END IF
         KKK = 0

C  TOGGLE NINE NODE ELEMENT GENERATION

      ELSE IF ((CIN(ICOM)(1:2) .EQ. 'NI') .OR.
     &   (CIN(ICOM)(1:2) .EQ. 'ni')) THEN
         ICOM = ICOM + 1
         IF (NINE) THEN
            NINE = .FALSE.
            CALL MESAGE (' ')
            CALL MESAGE ('NINE NODE QUAD GENERATION - DISABLED')
         ELSE
            NINE = .TRUE.
            EIGHT = .FALSE.
            CALL MESAGE (' ')
            CALL MESAGE ('NINE NODE QUAD GENERATION - ENABLED')
         END IF
         KKK = 0

C  EXIT OPTION - EXITS FASTQ

      ELSE IF ((CIN(ICOM)(1:2) .EQ. 'EX') .OR.
     &   (CIN(ICOM)(1:2) .EQ. 'ex')) THEN
         CALL STRLNG (CIN(ICOM), LEN)
         IF (LEN .GT. 1) THEN
            IF ( (CIN(ICOM)(2:2) .EQ. 'X') .OR.
     &         (CIN(ICOM)(2:2) .EQ. 'x')) THEN
               ICOM = ICOM + 1
               CALL FEXIT (WROTE, MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN,
     &            TIME1, BATCH, VERSN)
            ELSE
               ICOM = ICOM + 1
               CALL HELP_FQ(12)
            ENDIF
         ELSE
            ICOM = ICOM + 1
            CALL FEXIT (WROTE, MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN,
     &         TIME1, BATCH, VERSN)
         ENDIF
         GOTO 100

C  ENTER LINE INTERVALS

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'I') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'i')) THEN
         ICOM = ICOM + 1
         IF (ICOM.GT.JCOM) THEN
            CALL MESAGE ('ENTER LINE INTERVALS IN THE FOLLOWING '//
     &         'FORMAT:')
            CALL MESAGE ('[ LINE NO. (OR NEG SIDE NO.),  INTERVALS ]')
            CALL MESAGE ('HIT RETURN TO END INPUT')
         END IF
  110    CONTINUE
         IF (ICOM.GT.JCOM) THEN
            CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN, CIN, IIN,
     &         RIN)
            ICOM = 1
         END IF
         CALL GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2, IFOUND)
         IF (IFOUND.GT.0) THEN
            III(1) = I1
            CALL ININTR (ML, MS, 1, I2, III(1), N(19), N(20), NINT,
     &         NLPS, IFLINE, ILLIST, LINKL, LINKS, ADDLNK)
            GOTO 110
         END IF

C  SPAWN A PROCESS

      ELSE IF ((CIN(ICOM)(1:2) .EQ. 'SP') .OR.
     &   (CIN(ICOM)(1:2) .EQ. 'sp')) THEN
         ICOM = ICOM + 1
         CALL SPAWN (VAXVMS)

C  ENTER A NEW SIZE FOR A REGION

      ELSE IF ((CIN(ICOM)(1:2)  .EQ.  'SI') .OR.
     &   (CIN(ICOM)(1:2)  .EQ.  'si')) THEN
         ICOM = ICOM + 1
         CALL MESAGE ('NOTE: ENTER A DEFAULT SIZE BY SPECIFYING')
         CALL MESAGE ('      A SIZE WITH NO REGIONS.')
         CALL MESAGE ('ENTER REGION SIZE DATA IN THE FOLLOWING FORMAT:')
         CALL MESAGE ('[ SIZE, REGION 1, REGION 2, ..., REGION N ]')
         CALL MESAGE ('HIT RETURN TO END INPUT')
         ICOM = JCOM + 1
  120    CONTINUE
         CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, IFOUND, KIN, CIN, IIN,
     &      RIN)
         IF (IFOUND .GT. 0) THEN
            IF (IFOUND .LT. 2) THEN
               DEFSIZ = RIN(1)
            ELSE
               ADDLNK = .FALSE.
               DO 130 IRSZ = 2, IFOUND
                  CALL LTSORT(MR, LINKR, IIN(IRSZ), JJ, ADDLNK)
                  IF (JJ .GE. 0) THEN
                     RSIZE(JJ) = AMAX1(RIN(1), 0.)
                  ELSE
                     WRITE(*, 10000)IIN(IRSZ)
                  END IF
  130          CONTINUE
            END IF
            GO TO 120
         END IF

C  GENERATE THE MESH

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'P') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'p') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'S') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 's')) THEN
         DO 140 I = 1, 2
            IF (N(I).LE.0) THEN
               CALL MESAGE ('***************************************')
               IF (I .EQ. 1) THEN
                  CALL MESAGE
     &               ('*    NO POINT CARDS IN DATABASE       *')
               ELSE
                  CALL MESAGE
     &               ('*     NO LINE CARDS IN DATABASE       *')
               END IF
               CALL MESAGE ('*     NO MESH GENERATION POSSIBLE     *')
               CALL MESAGE ('***************************************')
               ICOM = ICOM + 1
               GOTO 100
            END IF
  140    CONTINUE
         IF ((CIN(ICOM)(1:1) .EQ. 'S') .OR.
     &      (CIN(ICOM)(1:1) .EQ. 's')) THEN
            STEP = .TRUE.
         ELSE
            STEP = .FALSE.
         END IF
         ICOM = ICOM + 1

C  OPEN THE TEMPORARY FILE

         IUNIT = 99
         OPEN (UNIT = IUNIT, STATUS = 'scratch', FORM = 'unformatted',
     &      ACCESS = 'sequential')
         REWIND IUNIT

C  GENERATE THE MESH

         IF (NPREGN.GT.0) THEN
            CALL MDDEL ('IPART')
            CALL MDDEL ('LSTNBC')
            CALL MDDEL ('LSTSBC')
            CALL MDDEL ('NNFLG')
            CALL MDDEL ('NNPTR')
            CALL MDDEL ('NNLEN')
            CALL MDDEL ('NSFLG')
            CALL MDDEL ('NSPTR')
            CALL MDDEL ('NSLEN')
            CALL MDDEL ('NVPTR')
            CALL MDDEL ('NVLEN')
            CALL MDDEL ('NSIDEN')
            CALL MDDEL ('NUID')
            CALL MDDEL ('XN')
            CALL MDDEL ('YN')
            CALL MDDEL ('NXK')
            CALL MDDEL ('MAT')
            CALL MDDEL ('KXN')
            CALL MDDEL ('LIST')
            CALL MDDEL ('LA')
            CALL MDDEL ('LB')
            CALL MDDEL ('CENTK')
            CALL MDDEL ('MATMAP')
            CALL MDDEL ('LISTN')
            CALL MDDEL ('WTNODE')
            CALL MDDEL ('WTSIDE')
            CALL MDDEL ('WTHOLD')
            CALL MDDEL ('IHERE')
            CALL MDDEL ('ILIST')
            CALL MDDEL ('XLIST')
            IF (.NOT. REMESH) THEN
               CALL MDDEL ('AMESUR')
               CALL MDDEL ('XNOLD')
               CALL MDDEL ('YNOLD')
               CALL MDDEL ('NXKOLD')
               CALL MDDEL ('MMPOLD')
               CALL MDDEL ('LINKEG')
               CALL MDDEL ('LISTEG')
               CALL MDDEL ('BMESUR')
               CALL MDDEL ('CMESUR')
               CALL MDDEL ('DMESUR')
            END IF
            CALL MDSTAT (NERR, MUSED)
            IF (NERR.GT.0) THEN
               CALL MDEROR (6)
               STOP ' '
            END IF
            NPREGN = 0

C  SET UP THE ARRAYS FOR ADAPTIVE REMESHING

         END IF
         IF (.NOT. REMESH) THEN
            CALL MDRSRV ('AMESUR', K(31), 1)
            CALL MDRSRV ('XNOLD',  K(32), 1)
            CALL MDRSRV ('YNOLD',  K(33), 1)
            CALL MDRSRV ('NXKOLD', K(34), 4)
            CALL MDRSRV ('MMPOLD', K(35), 3)
            CALL MDRSRV ('LINKEG', K(36), 2)
            CALL MDRSRV ('LISTEG', K(37), 2)
            CALL MDRSRV ('BMESUR', K(38), 1)
            CALL MDRSRV ('CMESUR', K(39), 1)
            CALL MDRSRV ('DMESUR', K(40), 1)
            NPROLD = 1
            NPNOLD = 1
            NPEOLD = 1
            NNXK = 1
            CALL MDSTAT (NERR, MUSED)
            IF (NERR.GT.0) THEN
               CALL MDEROR (6)
               STOP ' '
            END IF
         END IF

         LGROUP = .FALSE.
         DO 150 I = 1, N(7)
            IF (IRGFLG(I) .GE. 1) THEN
               LGROUP = .TRUE.
               GO TO 160
            END IF
  150    CONTINUE
  160    CONTINUE
         CALL QMESH (A, IA, MP, ML, MS, MR, MSC, MCOM, ICOM, JCOM, CIN,
     &      RIN, IIN, KIN, IUNIT, IDUMP, N, IPOINT, COOR, IPBOUN,
     &      ILINE, LTYPE, NINT, FACTOR, LCON, ILBOUN, ISBOUN, ISIDE,
     &      NLPS, IFLINE, ILLIST, IBARST, JMAT, JCENT, NLPB, JFLINE,
     &      JLLIST, IREGN, IMAT, NSPR, IFSIDE, ISLIST, IRPB, IPBF,
     &      NPPF, IFPB, LISTPB, ILBF, NLPF, IFLB, LISTLB, ISBF, NSPF,
     &      IFSB, LISTSB, LINKP, LINKL, LINKS, LINKB, LINKR, LINKSC,
     &      LINKPB, LINKLB, LINKSB, RSIZE, IFHOLE, NHPR, IHLIST,
     &      IRGFLG, ISCHM, SCHEME, DEFSCH, DEFSIZ, NPREGN, NPNBC, NPSBC,
     &      NPNODE, NPELEM, MAXKXN, STEP, DEV1, THREE, EIGHT, NINE,
     &      LGROUP, BATCH, A(K(31)), A(K(32)), A(K(33)), A(K(34)),
     &      A(K(35)), A(K(36)), A(K(37)), A(K(38)), MLINK, NPROLD,
     &      NPNOLD, NPEOLD, NNXK, REMESH, REXMIN, REXMAX, REYMIN,
     &      REYMAX, IDIVIS, SIZMIN, EMAX, EMIN)

         IF (REMESH) REMESH = .FALSE.
         IF (NPREGN.GT.0) THEN

C  SET UP THE NECESSARY DIMENSIONING FOR NUMBERING

C  A(K(1))  =  IPART   =  ARRAY OF BEGINNING/ENDING ELEMENT NUMBERS/REGION
C  A(K(2))  =  LSTNBC  =  LIST OF NODAL BOUNDARY CONDITIONS
C                      (REORDERED TO NODES)
C  A(K(3))  =  LSTSBC  =  LIST OF SIDE BOUNDARY CONDITIONS
C                      (REORDERED TO NELEMS)
C  A(K(4))  =  NNFLG   =  LIST OF NODAL FLAGS
C  A(K(5))  =  NNPTR   =  POINTERS INTO THE NODES LIST FOR EACH NODE FLAG
C  A(K(6))  =  NNLEN   =  NUMBER OF NODES IN THE LIST FOR EACH NODE FLAG
C  A(K(7))  =  NSFLG   =  LIST OF SIDE FLAGS
C  A(K(8))  =  NSPTR   =  POINTERS INTO THE SIDES LIST FOR EACH SIDE FLAG
C  A(K(9))  =  NSLEN   =  NUMBER OF SIDES IN THE LIST FOR EACH SIDE FLAG
C  A(K(10)) =  NVPTR   =  POINTERS INTO THE NSIDEN LIST FOR EACH SIDE FLAG
C  A(K(11)) =  NVLEN   =  NUMBER OF NODES IN THE NSIDEN LIST FOR EACH
C                     SIDE FLAG
C  A(K(12)) =  NSIDEN  =  LIST OF SIDE NODES ASSOCIATED WITH THE SIDE LIST
C  A(K(13)) =  NUID    =  NODE UNIQUE IDENTIFIER ARRAY
C                     (LATER - MAPDXG ARRAY)
C  A(K(14)) =  XN      =  X COORDINATE OF THE NODES
C  A(K(15)) =  YN      =  Y COORDINATE OF THE NODES
C  A(K(16)) =  NXK     =  NODES PER ELEMENT (CONNECTIVITY) ARRAY
C  A(K(17)) =  MAT     =  MATERIAL NUMBER FOR EACH ELEMENT
C  A(K(18)) =  KXN     =  ELEMENTS ATTACHED TO EACH NODE
C                     (LATER - LOOKUP TABLE)
C  A(K(19)) =  LIST    =  WORKING ARRAY FOR NODE NUMBERING
C                     (LATER - MAPGXD ARRAY)
C  A(K(20)) =  LA      =  WORKING ARRAY FOR NODE NUMBERING
C  A(K(21)) =  LB      =  WORKING ARRAY FOR NODE NUMBERING
C  A(K(22)) =  CENTK   =  ARRAY CONTAINING X, Y LOCATION OF
C                           ELEMENT CENTER
C  A(K(23)) =  MATMAP  =  ARRAY OF BEGIN/END ELEMENT NUMBERS/REGION IN
C                     MAPGXD
C  A(K(24)) =  LISTN   =  LIST FOR RENUMBERING THE NODES FOR
C                     BANDWIDTH OPTIMIZATION
C  A(K(25)) =  WTNODE  =  NODAL BOUNDARY FLAG WEIGHTING FACTOR ARRAY
C  A(K(26)) =  WTSIDE  =  ELEMENT SIDE BOUNDARY FLAG WEIGHTING
C                           FACTOR ARRAY
C  A(K(27)) =  WTHOLD  =  TEMPORARY ARRAY USED FOR MIDSIDE NODE
C                           PROCESSING
C  A(K(28)) =  IHERE   =  A WORK ARRAY FOR SORTING BOUNDARY FLAGS
C  A(K(29)) =  ILIST   =  LINE LIST FOR WEIGHTED BOUNDARIES
C  A(K(30)) =  XLIST   =  WORKING LIST FOR WEIGHTED BOUNDARIES
C  A(K(31)) =  AMESUR  =  FIRST ADAPTIVE MESHING VARIABLE
C  A(K(32)) =  XNOLD   =  OLD XN ARRAY USED DURING ADAPTIVE MESHING
C  A(K(33)) =  YNOLD   =  OLD YN ARRAY USED DURING ADAPTIVE MESHING
C  A(K(34)) =  NXKOLD  =  OLD NXK ARRAY USED DURING ADAPTIVE MESHING
C  A(K(35)) =  MMPOLD  =  OLD MATMAP ARRAY USED DURING ADAPTIVE MESHING
C  A(K(36)) =  LINKEG  =  ELEMENT GRID LINK FOR ADAPTIVE MESHING
C  A(K(37)) =  LISTEG  =  ELEMENT GRID LIST FOR ADAPTIVE MESHING
C  A(K(38)) =  BMESUR  =  AMESUR VALUES AVERAGED AT THE NODES
C  A(K(31)) =  CMESUR  =  SECOND ADAPTIVE MESHING VARIABLE
C  A(K(40)) =  DMESUR  =  CMESUR VALUES AVERAGED AT THE NODES

            IF (EIGHT) THEN
               NPNBC = NPNBC*2
               NPSBC = NPSBC*2
               NNXK = 8
            ELSE IF (NINE) THEN
               NPNBC = NPNBC*2
               NPSBC = NPSBC*2
               NNXK = 9
            ELSE
               NNXK = 4
            END IF
            NPWTS = MAX0 (NPNBC, NPSBC)
            MXNFLG = N(11) + N(13) + 1
            MXSFLG = N(15) + 1
            NNUID = MAX0 (NPNODE, NPNBC, NPSBC, NPELEM)
            MXLPS = 2
            DO 170 I = 1, N(3)
               MXLPS = MAX0(NLPS(I) + 1, MXLPS)
  170       CONTINUE
            CALL MDRSRV ('IPART', K(1), 3*NPREGN)
            CALL MDRSRV ('LSTNBC', K(2), NPNBC)
            CALL MDRSRV ('LSTSBC', K(3), NPSBC)
            CALL MDRSRV ('NNFLG', K(4), MXNFLG)
            CALL MDRSRV ('NNPTR', K(5), MXNFLG)
            CALL MDRSRV ('NNLEN', K(6), MXNFLG)
            CALL MDRSRV ('NSFLG', K(7), MXSFLG)
            CALL MDRSRV ('NSPTR', K(8), MXSFLG)
            CALL MDRSRV ('NSLEN', K(9), MXSFLG)
            CALL MDRSRV ('NVPTR', K(10), MXSFLG)
            CALL MDRSRV ('NVLEN', K(11), MXSFLG)
            CALL MDRSRV ('NSIDEN', K(12), NPSBC)
            CALL MDRSRV ('NUID', K(13), NNUID)
            CALL MDRSRV ('XN', K(14), NPNODE)
            CALL MDRSRV ('YN', K(15), NPNODE)
            CALL MDRSRV ('NXK', K(16), NPELEM*NNXK)
            CALL MDRSRV ('MAT', K(17), NPELEM)
C ... The kxn array is used for a work array in renum and needs to hold npsbc items.
            if (npsbc .gt. maxkxn) maxkxn = npsbc
            CALL MDRSRV ('KXN', K(18), MAXKXN*NNXK)
            CALL MDRSRV ('LIST', K(19), NNUID)
            CALL MDRSRV ('LA', K(20), NPNODE)
            CALL MDRSRV ('LB', K(21), NPNODE)
            CALL MDRSRV ('CENTK', K(22), NPELEM*2)
            CALL MDRSRV ('MATMAP', K(23), NPREGN*3)
            CALL MDRSRV ('LISTN', K(24), NNUID)
            CALL MDRSRV ('WTNODE', K(25), NPNBC)
            CALL MDRSRV ('WTSIDE', K(26), NPSBC)
            CALL MDRSRV ('WTHOLD', K(27), NPWTS)
            CALL MDRSRV ('IHERE',  K(28), NNUID)
            CALL MDRSRV ('ILIST',  K(29), MXLPS)
            CALL MDRSRV ('XLIST',  K(30), MXLPS)
            CALL MDSTAT (NERR, MUSED)
            IF (NERR.GT.0) THEN
               CALL MDEROR (6)
               STOP ' '
            END IF

            if (n(5) .gt. 0) then
              isbars = .true.
            else
              isbars = .false.
            end if

            CALL RENUM (NPNODE, NPELEM, MXNFLG, MXSFLG, NPNBC, NPSBC,
     &         NPWTS, NPREGN, MP, ML, MS, MR, MSC, MAXKXN, NNUID, NNXK,
     &         MXLPS, IUNIT, NNN, KKK, N(28), N(2), N(11), N(13),
     &         N(15), IA(K(1)), IA(K(2)), IA(K(3)), IA(K(4)),
     &         IA(K(5)), IA(K(6)), IA(K(7)), IA(K(8)), IA(K(9)),
     &         IA(K(10)), IA(K(11)), IA(K(12)), IA(K(13)),
     &         A(K(14)), A(K(15)), IA(K(16)), IA(K(17)),
     &         IA(K(18)), IA(K(19)), IA(K(20)), IA(K(21)),
     &         IA(K(23)), IA(K(24)), A(K(25)), A(K(26)),
     &         A(K(27)), IA(K(28)), IA(K(29)), A(K(30)), NUMMAT,
     &         NBCNOD, NNLIST, NBCSID, NSLIST, NVLIST, COOR, ILINE,
     &         LTYPE, LCON, ISIDE, NLPS, IFLINE, ILLIST, LINKP, LINKL,
     &         LINKS, LINKR, IMAT, LINKB, JMAT, IPBF, NPPF, IFPB,
     &         LISTPB, IWTPBF, ILBF, NLPF, IFLB, LISTLB, IWTLBF, ISBF,
     &         NSPF, IFSB, LISTSB, IWTSBF, LINKPB, LINKLB, LINKSB,
     &         NUMBER, THREE, EIGHT, NINE, OPTIM, ISBARS)

C  AN ERROR HAS OCCURRED AND THUS THE DUMMY REMESH ARRAYS MUST BE
C  DELETED IF NOT REMESHING

         ELSEIF (.NOT. REMESH) THEN
            CALL MDDEL ('AMESUR')
            CALL MDDEL ('XNOLD')
            CALL MDDEL ('YNOLD')
            CALL MDDEL ('NXKOLD')
            CALL MDDEL ('MMPOLD')
            CALL MDDEL ('LINKEG')
            CALL MDDEL ('LISTEG')
            CALL MDDEL ('BMESUR')
            CALL MDDEL ('CMESUR')
            CALL MDDEL ('DMESUR')
            CALL MDSTAT (NERR, MUSED)
            IF (NERR.GT.0) THEN
               CALL MDEROR (6)
               STOP ' '
            END IF
         END IF

         CLOSE (IUNIT)

C  ENTER LINE FACTORS

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'F') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'f')) THEN
         ICOM = ICOM + 1
         IF (ICOM.GT.JCOM) THEN
            CALL MESAGE ('ENTER LINE FACTORS IN THE FOLLOWING FORMAT:')
            CALL MESAGE ('[ LINE NO. (OR NEG. SIDE NO.),  FACTOR ]')
            CALL MESAGE ('HIT RETURN TO END INPUT')
         END IF
  180    CONTINUE
         IF (ICOM.GT.JCOM) THEN
            CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN, CIN, IIN,
     &         RIN)
            ICOM = 1
         END IF
         CALL GETI1R (MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN, I1, R1,
     &      IFOUND)
         IF (IFOUND.GT.0) THEN
            III(1) = I1
            CALL INFACT (ML, MS, 1, R1, III(1), N(19), N(20), FACTOR,
     &         NLPS, IFLINE, ILLIST, LINKL, LINKS, ADDLNK)
            GOTO 180
         END IF

C  CLEAR OUT ALL THE INTERVALS ASSIGNED TO LINES BY REGION

      ELSE IF ((CIN(ICOM)(1:2) .EQ. 'CI') .OR.
     &   (CIN(ICOM)(1:2) .EQ. 'ci')) THEN
         ICOM = ICOM + 1
         CALL MESAGE ('CLEAR INTERVALS FROM REGIONS FROM <I1> TO <I2>')
         CALL MESAGE ('HIT RETURN TO END INPUT')
  200    CONTINUE
         IF (ICOM.GT.JCOM)THEN
            CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN, CIN,
     &         IIN, RIN)
            ICOM = 1
         ENDIF
         CALL GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2, IFOUND)

         IF (IFOUND .GT. 0) THEN
            IF (I1 .GT. 0) THEN
               CALL CHECK (I1, I2, N (22))

C  REMOVE INTERVALS ON LINES ASSOCIATED WITH THE REGIONS

               DO 230 I = I1, I2
                  CALL LTSORT (MR, LINKR, I, II, ADDLNK)
                  IF (II .GT. 0) THEN
                     DO 220 J = IFSIDE (II), IFSIDE (II) + NSPR (II)-1

C  FIRST REMOVE INTERVALS OF LINES ON SIDE DATA

                        CALL LTSORT (MS, LINKS, ISLIST (J), JJ, ADDLNK)
                        IF ((ISLIST (J) .GT. 0) .AND. (JJ .GT. 0)) THEN
                           DO 210 KKK = IFLINE (JJ), IFLINE (JJ) +
     &                        NLPS (JJ)-1
                              CALL LTSORT (ML, LINKL, ILLIST (KKK), KK,
     &                           ADDLNK)
                              IF (KK .GT. 0) NINT (KK) = 0
  210                      CONTINUE

C  NEXT REMOVE INTERVALS ON LINES ALONE

                        ELSE
                           JJ = IABS (ISLIST (J))
                           CALL LTSORT (ML, LINKL, JJ, KK, ADDLNK)
                           IF (KK.GT.0) NINT (KK) = 0
                        ENDIF
  220                CONTINUE
                  ENDIF
  230          CONTINUE
               GOTO 200
            ENDIF
         ENDIF

C     ADJUST THE MESH

      ELSE IF ((CIN(ICOM)(1:2) .EQ. 'AD') .OR.
     &   (CIN(ICOM)(1:2) .EQ. 'ad')) THEN
         ICOM = ICOM + 1
         IF (KKK.LE.0) THEN
            CALL MESAGE (' ')
            CALL MESAGE ('****************************************')
            CALL MESAGE ('*      NO ELEMENTS IN DATABASE         *')
            CALL MESAGE ('*    NO MESH ADJUSTMENT POSSIBLE       *')
            CALL MESAGE ('****************************************')
            CALL MESAGE (' ')
            GOTO 100
         END IF
         CALL ADJMSH (MS, MR, NPNODE, NPELEM, MXNFLG, MXSFLG, NPREGN,
     &      NPNBC, NPSBC, MCOM, ICOM, JCOM, CIN, RIN, IIN, KIN, NNN,
     &      KKK, NNXK, IA(K(2)), IA(K(3)), IA(K(4)), IA(K(5)), IA(K(6)),
     &      IA(K(7)), IA(K(8)), IA(K(9)), IA(K(10)), IA(K(11)),
     &      IA(K(12)), IA(K(13)), A(K(14)), A(K(15)), IA(K(16)),
     &      IA(K(17)), IA(K(19)), IA(K(23)), A(K(25)), A(K(26)), NBCNOD,
     &      NNLIST, NBCSID, NSLIST, NVLIST, NUMMAT, LINKM, TITLE, ERR,
     &      EIGHT, NINE, VERSN)
         IF (ERR) THEN
            CALL MESAGE (' ')
            CALL MESAGE ('******************************************')
            CALL MESAGE ('*         ERROR ADJUSTING MESH           *')
            CALL MESAGE ('******************************************')
            CALL MESAGE (' ')
         END IF

C     CALCULATE A DISTORTION INDEX FOR THE REGIONS

      ELSE IF (((CIN(ICOM)(1:1) .EQ. 'D') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'd')) .AND.
     &   (CIN(ICOM)(2:2) .NE. 'N') .AND.
     &   (CIN(ICOM)(2:2) .NE. 'n')) THEN
         ICOM = ICOM + 1
         IF (KKK.LE.0) THEN
            CALL MESAGE (' ')
            CALL MESAGE ('****************************************')
            CALL MESAGE ('*      NO ELEMENTS IN DATABASE         *')
            CALL MESAGE ('*   NO DISTORTION INDEX CALCULATED     *')
            CALL MESAGE ('****************************************')
            CALL MESAGE (' ')
            GOTO 100
         END IF
         CALL RGDSTR (NPNODE, NPELEM, KKK, NNXK, IA(K(14)),
     &      A(K(15)), A(K(16)) )

C     Write out the mesh data into the genesis or exodusII data base
      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'W') .OR.
     &         (CIN(ICOM)(1:1) .EQ. 'w')) THEN
         ICOM = ICOM + 1
         IF (KKK.LE.0) THEN
            CALL MESAGE (' ')
            CALL MESAGE ('****************************************')
            CALL MESAGE ('*       NO ELEMENTS IN DATABASE        *')
            CALL MESAGE ('*  NO EXODUSII FILE WRITING POSSIBLE   *')
            CALL MESAGE ('****************************************')
            CALL MESAGE (' ')
            GOTO 100
         END IF
         IUNIT = 9
  240    CONTINUE
         IF (BATCH) THEN
            CALL EXNAME (IUNIT, FNAME, LN)
         ELSE
            IF (ICOM.LE.JCOM) THEN
               FNAME = CIN(ICOM)
               LN = lenstr(FNAME)
               ICOM = ICOM + 1
            ELSE
               CALL INQSTR ('OUTPUT DATABASE FILE NAME: ',
     &            FNAME)
               LN = lenstr(FNAME)
            END IF
         END IF
C ... Create the exodusII file
         if (EXODUSII) THEN
            CMPSIZ = 0
            IOWS   = 0

C ... See if user-specified output word size EXT05
            call exname (-5, hold, llen)
            if (llen .lt. 1) goto 25
            read(hold,'(i1)',ERR=25)iows
            goto 26
 25         continue

            call exparm (cdumh, cdums, idum, iows, idum, idum)
 26         continue

C ... One final check to make sure we have a valid iows.
            if (iows .ne. 4 .and. iows .ne. 8) then
              iows = 4
            endif

            iunit = excre(FNAME(:LN), EXCLOB, CMPSIZ, IOWS, IERR)
            if (ierr .lt. 0) then
               call exopts (EXVRBS, ierr)
               call exerr('fastq', 'Error from excre', ierr)
               if (batch) then
                 stop 'Exodus Error'
               else
                 go to 240
               endif
            endif
            CALL WREX2 (MS, MR, NPNODE, NPELEM, MXNFLG, MXSFLG, NPREGN,
     &           NPNBC, NPSBC, IUNIT, NNN, KKK, NNXK, IA(K(2)),
     &           IA(K(3)), IA(K(4)), IA(K(5)), IA(K(6)), IA(K(7)),
     &           IA(K(8)), IA(K(9)), IA(K(10)), IA(K(11)), IA(K(12)),
     &           IA(K(13)), A(K(14)), A(K(15)), IA(K(16)), IA(K(17)),
     &           IA(K(19)), IA(K(23)), A(K(25)), A(K(26)),
     &           NBCNOD, NNLIST, NBCSID, NSLIST, NVLIST, NUMMAT, LINKM,
     &           TITLE, ERR, EIGHT, NINE, VERSN, A, IA, FNAME(:LN))
            IF (ERR) THEN
               IF (BATCH) THEN
                  CALL MESAGE (' ')
                  CALL MESAGE
     &                 ('******************************************')
                  CALL MESAGE
     &                 ('*  ERROR WRITING GENESIS DATABASE FILE   *')
                  CALL MESAGE
     &                 ('*         NO OUTPUT FILE SAVED           *')
                  CALL MESAGE
     &                 ('******************************************')
                  CALL MESAGE (' ')
                  STOP
               END IF
            END IF
            call exclos(iunit, ierr)
         else
            OPEN (UNIT = IUNIT, FILE = FNAME(:LN), STATUS = 'UNKNOWN',
     &           FORM = 'UNFORMATTED', ACCESS = 'SEQUENTIAL', ERR = 240)
            REWIND IUNIT
            CALL WRGENS (MS, MR, NPNODE, NPELEM, MXNFLG, MXSFLG, NPREGN,
     &           NPNBC, NPSBC, IUNIT, NNN, KKK, NNXK, IA(K(2)),
     &           IA(K(3)), IA(K(4)), IA(K(5)), IA(K(6)), IA(K(7)),
     &           IA(K(8)), IA(K(9)), IA(K(10)), IA(K(11)), IA(K(12)),
     &           IA(K(13)), A(K(14)), A(K(15)), IA(K(16)), IA(K(17)),
     &           IA(K(19)), IA(K(23)), A(K(25)), A(K(26)),
     &           NBCNOD, NNLIST, NBCSID, NSLIST, NVLIST, NUMMAT, LINKM,
     &           TITLE, ERR, EIGHT, NINE, VERSN)
            IF (ERR) THEN
               CLOSE (UNIT = IUNIT, STATUS = 'DELETE')
               IF (BATCH) THEN
                  CALL MESAGE (' ')
                  CALL MESAGE
     &                 ('******************************************')
                  CALL MESAGE
     &                 ('*  ERROR WRITING GENESIS DATABASE FILE   *')
                  CALL MESAGE
     &                 ('*         NO OUTPUT FILE SAVED           *')
                  CALL MESAGE
     &                 ('******************************************')
                  CALL MESAGE (' ')
                  STOP
               END IF
            ELSE
               CLOSE (UNIT = IUNIT, STATUS = 'KEEP')
            end if
         end if

C     WRITE OUT THE MESH DATA INTO THE ERROR CODE DATA FORMAT

      ELSE IF ((CIN (ICOM) (1:1) .EQ. 'J') .OR.
     &   (CIN (ICOM) (1:1) .EQ. 'j')) THEN
         ICOM = ICOM + 1
         IF (KKK.LE.0) THEN
            CALL MESAGE (' ')
            CALL MESAGE ('****************************************')
            CALL MESAGE ('*      NO ELEMENTS IN DATABASE         *')
            CALL MESAGE ('*     NO FILE WRITING POSSIBLE        *')
            CALL MESAGE ('****************************************')
            CALL MESAGE (' ')
            GOTO 100
         END IF
         IUNIT = 99
  250    CONTINUE
         IF (ICOM.LE.JCOM) THEN
            FNAME = CIN (ICOM)
            ICOM = ICOM + 1
         ELSE
            CALL INQSTR ('JOE''S ERROR OUTPUT FILE NAME: ', FNAME)
         END IF
         OPEN (UNIT = IUNIT, FILE = FNAME, STATUS = 'NEW', ERR = 250)
         REWIND IUNIT
         CALL WRJERR (MS, MR, NPNODE, NPELEM, MXNFLG, MXSFLG, NPREGN,
     &      NPNBC, NPSBC, IUNIT, NNN, KKK, NNXK, IA(K(2)), IA(K(3)),
     &      IA(K(4)), IA(K(5)), IA(K(6)), IA(K(7)), IA(K(8)),
     &      IA(K(9)), IA(K(10)), IA(K(11)), IA(K(12)), IA(K(13)),
     &      A(K(14)), A(K(15)), IA(K(16)), IA(K(17)), IA(K(19)),
     &      IA(K(23)), NBCNOD, NNLIST, NBCSID, NSLIST, NVLIST, NUMMAT,
     &      LINKM, TITLE, ERR, EIGHT, NINE)
         IF(ERR) THEN
            CLOSE (UNIT = IUNIT, STATUS = 'DELETE')
         ELSE
            CLOSE (UNIT = IUNIT, STATUS = 'KEEP')
         END IF

C     WRITE OUT THE MESH DATA INTO THE ABAQUS DATA FORMAT

      ELSE IF ((CIN (ICOM) (1:1) .EQ. 'A') .OR.
     &   (CIN (ICOM) (1:1) .EQ. 'a')) THEN
         ICOM = ICOM + 1
         IF (KKK.LE.0) THEN
            CALL MESAGE (' ')
            CALL MESAGE ('****************************************')
            CALL MESAGE ('*      NO ELEMENTS IN DATABASE         *')
            CALL MESAGE ('*     NO FILE WRITING POSSIBLE        *')
            CALL MESAGE ('****************************************')
            CALL MESAGE (' ')
            GOTO 100
         END IF
         IUNIT = 99
  260    CONTINUE
         IF (ICOM.LE.JCOM) THEN
            FNAME = CIN (ICOM)
            ICOM = ICOM + 1
         ELSE
            CALL INQSTR ('ABAQUS OUTPUT FILE NAME: ', FNAME)
         END IF
         OPEN (UNIT = IUNIT, FILE = FNAME, STATUS = 'UNKNOWN',
     &        ERR = 260)
         REWIND IUNIT
         CALL WRABQS (MS, MR, NPNODE, NPELEM, MXNFLG, MXSFLG, NPREGN,
     &      NPNBC, NPSBC, IUNIT, NNN, KKK, NNXK, IA(K(2)), IA(K(3)),
     &      IA(K(4)), IA(K(5)), IA(K(6)), IA(K(7)), IA(K(8)),
     &      IA(K(9)), IA(K(10)), IA(K(11)), IA(K(12)), IA(K(13)),
     &      A(K(14)), A(K(15)), IA(K(16)), IA(K(17)), IA(K(19)),
     &      IA(K(23)), NBCNOD, NNLIST, NBCSID, NSLIST, NVLIST, NUMMAT,
     &      LINKM, TITLE, ERR, EIGHT, NINE)
         IF(ERR) THEN
            CLOSE (UNIT = IUNIT, STATUS = 'DELETE')
         ELSE
            CLOSE (UNIT = IUNIT, STATUS = 'KEEP')
         END IF

C     WRITE OUT THE MESH DATA INTO THE NASTRAN DATA FORMAT

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'N') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'n') .OR. (CIN(ICOM)(1:2) .EQ. 'DN') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'DN')) THEN
         ICOM = ICOM + 1
         IF ((CIN(ICOM)(1:1) .EQ. 'N') .OR.
     &      (CIN(ICOM)(1:1) .EQ. 'n')) THEN
            LONG = .TRUE.
         ELSE
            LONG = .FALSE.
         ENDIF
         IF (KKK.LE.0) THEN
            CALL MESAGE (' ')
            CALL MESAGE ('****************************************')
            CALL MESAGE ('*      NO ELEMENTS IN DATABASE         *')
            CALL MESAGE ('*     NO FILE WRITING POSSIBLE        *')
            CALL MESAGE ('****************************************')
            CALL MESAGE (' ')
            GOTO 100
         END IF
         IUNIT = 99
  270    CONTINUE
         IF (ICOM.LE.JCOM) THEN
            FNAME = CIN(ICOM)
            ICOM = ICOM + 1
         ELSE
            CALL INQSTR ('NASTRAN OUTPUT FILE NAME: ', FNAME)
         END IF
         OPEN (UNIT = IUNIT, FILE = FNAME, STATUS = 'UNKNOWN',
     &        ERR = 270)
         REWIND IUNIT
         CALL WRNAST (MS, MR, NPNODE, NPELEM, MXNFLG, MXSFLG, NPREGN,
     &      NPNBC, NPSBC, IUNIT, NNN, KKK, NNXK, IA(K(2)), IA(K(3)),
     &      IA(K(4)), IA(K(5)), IA(K(6)), IA(K(7)), IA(K(8)),
     &      IA(K(9)), IA(K(10)), IA(K(11)), IA(K(12)), IA(K(13)),
     &      A(K(14)), A(K(15)), IA(K(16)), IA(K(17)), IA(K(19)),
     &      IA(K(23)), NBCNOD, NNLIST, NBCSID, NSLIST, NVLIST, NUMMAT,
     &      LINKM, TITLE, ERR, EIGHT, NINE, LONG)
         IF (ERR) THEN
            CLOSE (UNIT = IUNIT, STATUS = 'DELETE')
         ELSE
            CLOSE (UNIT = IUNIT, STATUS = 'KEEP')
         END IF
      ELSE
         ICOM = ICOM + 1
         CALL HELP_FQ (12)
      END IF
      GOTO 100

10000 FORMAT (' REGION NO:', I5, ' IS NOT IN THE DATABASE', /,
     &   ' THUS NO SIZE CAN BE ENTERED')

      END
