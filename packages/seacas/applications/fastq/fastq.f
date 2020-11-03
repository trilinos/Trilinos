C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      PROGRAM FASTQ
C***********************************************************************

C  FASTQ = A PROGRAM TO QUICKLY GENERATE QUADRALATERAL MESHES

C***********************************************************************

C                       WRITTEN AND MAINTAINED BY
C                             TED D. BLACKER
C                             DIVISION 1523
C                              VERSION 1.4X

C***********************************************************************

C                          USES WORK PREVIOUSLY
C                              COMPLETED BY
C                            RONDALL E. JONES
C                             DIVISION 2644
C                     (QMESH, RENUM, AND QNUM CODES)

C***********************************************************************

C  NOTE:  FASTQ CALLS SEVERAL GRAPHICS ROUTINES FROM THE PLT PLOT
C         PACKAGE, AS WELL AS A NUMBER OF UTILITY ROUTINES FROM
C         THE DEPARTMENT 1520 SUPES LIBRARY OF ROUTINES.  OF PRIME
C         USE IS THE FREE FIELD READER ROUTINES AND THE DYNAMIC
C         MEMORY ALLOCATION ROUTINES.

C***********************************************************************

C  VARIABLES USED:
C     IANS   = LOGICAL RESPONSE FROM YES-NO QUESTION
C     TITLE  = MESH TITLE
C     DRWTAB = .TRUE. IF THE DIGITIZER HAS BEEN INITIALIZED TO A DRAWING
C     TBZOOM = .TRUE. IF THE DIGITIZER ZOOM HAS BEEN SET
C     WROTE  = .FALSE. IF THE DATA HAS BEEN CHANGED SINCE THE LAST WRITE
C     K      = POINTER ARRAY TO DIMENSION A
C     A      = LARGE ARRAY FOR VARIABLE DIMENSIONING
C              NOTE: FOR DYNAMIC VARIABLE DIMENSIONING, THE ROUTINES
C              WILL WORK AS CONFIGURED - THE SWITCH TO NORMAL VARIABLE
C              DIMENSIONING IS NOTED IN THE CODE AS IT IS NEEDED.
C              DYNAMIC VARIABLE DIMENSIONING IS NOT MACHINE INDEPENDENT.
C     MERGE  = .TRUE. IF TWO DATA FILES ARE TO BE MERGED
C     NOROOM = .TRUE. IF MORE ROOM IS NEEDED TO INPUT THE DATA
C     BATCH  = .TRUE. IF THE PROGRAM IS BEING RUN IN BATCH MODE
C     START  = .TRUE. IF THE PROGRAM IS JUST STARTING - TRY A READ FIRST
C     VAXVMS = .TRUE. IF THE INSTALLATION IS ON A VAX/VMS MACHINE
C              (IT IS ASSUMED THAT VAXVMS WILL HAVE MULTIPLE VIRTUAL
C              DEVICE DRIVERS AVAILABLE - MVDI)

C***********************************************************************

      LOGICAL IANS, DRWTAB, WROTE, OPTIM, MERGE, NOROOM, TBZOOM
      LOGICAL LABP, LABL, LABR, AXISD, LABMD, LABI, LABF, LABPB, LABLB,
     &   LABSBD
      LOGICAL FULL, LABSC, LABSZ, AREACG, AXIS, AXIST
      LOGICAL LABE, LABO, LABN, LABNB, LABSB, LABM, LABW
      LOGICAL BATCH, VAXVMS, START, ALPHA, HARDPL, SNAP
      LOGICAL THREE, EIGHT, NINE, REGWRT, BARWRT
      LOGICAL EXODUSII

      PARAMETER (MSC = 60, MA = 4, MCOM = 50, MSNAP = 100)

C  NOTE:  IF DYNAMIC VARIABLE DIMENSIONING IS NOT BEING USED, THIS
C         PARAMETER STATEMENT SHOULD BE EXPANDED TO THE FORM:
C      PARAMETER (MP = 100, ML = 100, MS = 50, MR = 30, MSC = 30, MA = 4)
C         IF LARGER DIMENSIONS ARE DESIRED, MP, ML, MS, AND MR CAN
C         BE INCREASED ACCORDINGLY.
C         ALSO, THE VARIABLE A SHOULD BE DIMENSIONED AS:
C            DIMENSION A(MP*17 + ML*31 + MS*10 + MR*15)

      CHARACTER  DEV1*3, DEV2*3, VERSN*10, NUMBER*80, DATATYPE*8
      CHARACTER*8  HARD, SOFT, DATE, TIME
      CHARACTER*72 SCHEME, DEFSCH, TITLE, CIN(MCOM)
      CHARACTER*8  MEMDBG
      CHARACTER*2048 FNAME

      DIMENSION K(67), N(29), ISCHM(MSC), SCHEME(MSC), NUMBER(MSC)
      DIMENSION IDEV(2), SNAPDX(2,MSNAP), NSNAP(2)
      DIMENSION KIN(MCOM), IIN(MCOM), RIN(MCOM)
      DIMENSION A(1), IA(1)

      EQUIVALENCE (A, IA)

C  INITIALIZE VARIABLES

C ... By default, write exodusII format. If the environment variable
C     EXT04 is set to 'EXODUSII', write exodusII format.  If EXT04
C     is set to 'GENESIS', write exodusI format
      call exname(-4, datatype, klen)
      CALL EXUPCS (datatype(:klen))
      if (datatype(:8) .eq. 'EXODUSII') then
         exodusII = .TRUE.
      else if (datatype(:7) .eq. 'GENESIS') then
         exodusII = .FALSE.
      else
         exodusII = .TRUE.
      end if

      CALL EXCPUS (TIME1)
      TITLE = ' '
      DEFSCH = 'M'
      VERSN = 'FASTQ 3.22'
      DEFSIZ = 0.
      SNAP = .TRUE.
      TBZOOM = .FALSE.
      MERGE = .FALSE.
      NSNAP(1) = 0
      NSNAP(2) = 0
      DO 100 I = 1, MSNAP
         SNAPDX(1, I) = 0.0
         SNAPDX(2, I) = 0.0
  100 CONTINUE

      three = .false.
      eight = .false.
      nine  = .false.
      optim = .false.

C  GET THE CURRENT SYSTEM PARAMETERS AND SET MODE FOR RUNNING

      CALL EXPARM (HARD, SOFT, MODE, KCSU, KNSU, IDAU)

C**               FILE MODE

      call exname (-3, dev2, klen)
      if (dev2(:3) .eq. 'BAT') THEN
         BATCH = .TRUE.
         ALPHA = .TRUE.
         START = .FALSE.
         CIN(1) = 'READ'
         CIN(2) = 'MESH'
         CIN(3) = 'PROCESS'
         CIN(4) = 'WRITE'
         CIN(5) = 'EXIT'
         ICOM = 1
         JCOM = 5
      ELSE
         BATCH = .FALSE.
         ALPHA = .FALSE.
         START = .TRUE.
         CIN(1) = 'READ'
         ICOM = 1
         JCOM = 1
      END IF

      IF (SOFT(1:3) .EQ. 'VMS') THEN
         VAXVMS = .TRUE.
      ELSE
         VAXVMS = .FALSE.
      ENDIF

      IF (HARD(1:4) .EQ. 'CRAY' .AND. SOFT(1:3) .NE. 'UNI') THEN

C  WE MUST NOW INPUT THE DEVICE IN AN ADHOCK MANNER FROM THE CRAY

         CALL EXNAME (75, FNAME, LEN)
         DEV1 = FNAME (1:3)
         CALL EXUPCS (DEV1)
      END IF

C***********************************************************************

      TITLE = ' '
      DEFSCH = 'M'
      WROTE = .TRUE.
      DRWTAB = .FALSE.
      HARDPL = .FALSE.
      DO 110 I = 1, 29
         N(I) = 0
  110 CONTINUE

C  SET UP THE DEFAULT LABELING FOR DATA PLOTTING

      AREACG = .FALSE.
      AXIS = .FALSE.
      AXISD = .FALSE.
      AXIST = .FALSE.
      FULL = .FALSE.
      LABP = .TRUE.
      LABL = .TRUE.
      LABR = .TRUE.
      LABMD = .FALSE.
      LABI = .FALSE.
      LABF = .FALSE.
      LABPB = .FALSE.
      LABLB = .FALSE.
      LABSBD = .FALSE.
      LABSC = .FALSE.
      LABSZ = .FALSE.
      LABN = .FALSE.
      LABE = .FALSE.
      LABO = .FALSE.
      LABNB = .FALSE.
      LABSB = .FALSE.
      LABM = .FALSE.
      LABW = .FALSE.

C  PRINT GREETING AND TRACE

      CALL MESAGE (' ')
      CALL MESAGE ('WELCOME TO FASTQ:')
      CALL EXDATE (DATE)
      CALL EXTIME (TIME)
      WRITE (*, *) '            DATE: ', DATE
      WRITE (*, *) '            TIME: ', TIME
      WRITE (*, *) '         VERSION: ', VERSN
      if (exodusII) then
         write (*,*) '   Output Format: ExodusII'
      else
         write (*,*) '   Output Format: Genesis/ExodusI'
      end if

      CALL MESAGE (' ')
      WRITE (*, *)
     *  '+++                Copyright 2014 NTESS                   +++'
      WRITE (*, *)
     *  '+++ The U.S. Government retains a limited license in this +++'
      WRITE (*, *)
     *  '+++    software as prescribed in AL 88-1 and AL 91-7.     +++'
      WRITE (*, *)
     *  '+++ Export of this program may require a license from the +++'
      WRITE (*, *)
     *  '+++               United States Government                +++'

C  IF THE CODE IS BEING RUN ON THE VAX INTERACTIVELY,
C  GET WHICH DEVICE IS BEING USED
C  AND SET UP THE MULTIPLE DEVICE OUTPUT ROUTINES

      IF ((VAXVMS) .AND. (.NOT.BATCH)) THEN
         CALL EXNAME (-1, DEV1, LEN)
         CALL EXNAME (-2, DEV2, LEN)
         CALL VDIQES (10001, KAVAL1)
         CALL VDIQES (10002, KAVAL2)
         IF (KAVAL1.NE.1) THEN
            ALPHA = .TRUE.
            CALL MESAGE ('TERMINAL PLOTTING DEVICE NOT AVAILABLE')
         ELSE
            ALPHA = .FALSE.
         END IF
         IF (KAVAL2.NE.1) CALL MESAGE ('HARDCOPY DEVICE NOT AVAILABLE')
      END IF
      IF ((.NOT.BATCH) .AND. (.NOT.ALPHA)) THEN
         CALL VDESCP (10003, 0, 0)
         CALL PLTINT
         CALL VDESCP (10001, 0, 0)
         CALL PLTSTV (2, 160.)
      END IF

C  SET UP THE DUMP LOCATION FOR THE LOG FILE

      IDUMP = 0

C-----------------------------------------------------------------------

C  THE NEXT SERIES OF STATEMENTS MUST BE TAKEN OUT IF NOT USING
C  DYNAMIC VARIABLE DIMENSIONING

C  SET UP THE INITIAL POINTER ARRAY SYSTEM

      MP = 1000
      ML = 1000
      MS = 1000
      MR = 1000

C  INITIALIZE THE DYNAMIC DIMENSIONING ROUTINES

      CALL MDINIT (A)
      CALL MDFILL(0)

C ... See if supes memory debugging desired
C     If EXT99 Environment variable set, turn on supes memory debugging
C     The numeric value of the variable is used as the unit to write
C     debug information to.
      CALL EXNAME (-99, MEMDBG, L)
      IF (L .GE. 1) THEN
        READ(MEMDBG(:L), '(I8)', ERR=20) IUNIT
        CALL MDDEBG(IUNIT)
      END IF
 20   CONTINUE

C  GET INITIAL SPACE IN ARRAY A

      CALL MDRSRV ('IPOINT', K(1), MP)
      CALL MDRSRV ('COOR', K(2), MP*2)
      CALL MDRSRV ('IPBOUN', K(3), MP)
      CALL MDRSRV ('ILINE', K(4), ML)
      CALL MDRSRV ('LTYPE', K(5), ML)
      CALL MDRSRV ('NINT', K(6), ML)
      CALL MDRSRV ('FACTOR', K(7), ML)
      CALL MDRSRV ('LCON', K(8), ML*3)
      CALL MDRSRV ('ILBOUN', K(9), ML)
      CALL MDRSRV ('ISBOUN', K(10), ML)
      CALL MDRSRV ('ISIDE', K(11), MS)
      CALL MDRSRV ('NLPS', K(12), MS)
      CALL MDRSRV ('IFLINE', K(13), MS)
      CALL MDRSRV ('ILLIST', K(14), MS*3)
      CALL MDRSRV ('IBARST', K(15), MS)
      CALL MDRSRV ('JMAT', K(16), MS)
      CALL MDRSRV ('JCENT', K(17), MS)
      CALL MDRSRV ('NLPB', K(18), MS)
      CALL MDRSRV ('JFLINE', K(19), MS)
      CALL MDRSRV ('JLLIST', K(20), MS*3)
      CALL MDRSRV ('IREGN', K(21), MR)
      CALL MDRSRV ('IMAT', K(22), MR)
      CALL MDRSRV ('NSPR', K(23), MR)
      CALL MDRSRV ('IFSIDE', K(24), MR)
      CALL MDRSRV ('ISLIST', K(25), MR*4)
      CALL MDRSRV ('IRPB', K(26), MR)
      CALL MDRSRV ('IPBF', K(27), MP)
      CALL MDRSRV ('NPPF', K(28), MP)
      CALL MDRSRV ('IFPB', K(29), MP)
      CALL MDRSRV ('LISTPB', K(30), MP*2)
      CALL MDRSRV ('ILBF', K(31), ML)
      CALL MDRSRV ('NLPF', K(32), ML)
      CALL MDRSRV ('IFLB', K(33), ML)
      CALL MDRSRV ('LISTLB', K(34), ML*2)
      CALL MDRSRV ('ISBF', K(35), ML)
      CALL MDRSRV ('NSPF', K(36), ML)
      CALL MDRSRV ('IFSB', K(37), ML)
      CALL MDRSRV ('LISTSB', K(38), ML*2)
      CALL MDRSRV ('ATTRIB', K(39), MA*(MR + MS))
      CALL MDRSRV ('LINKP', K(40), MP*2)
      CALL MDRSRV ('LINKL', K(41), ML*2)
      CALL MDRSRV ('LINKS', K(42), MS*2)
      CALL MDRSRV ('LINKB', K(43), MS*2)
      CALL MDRSRV ('LINKR', K(44), MR*2)
      CALL MDRSRV ('LINKM', K(45), (MS + MR)*2)
      CALL MDRSRV ('LINKSC', K(46), MR*2)
      CALL MDRSRV ('LINKPB', K(47), MP*2)
      CALL MDRSRV ('LINKLB', K(48), ML*2)
      CALL MDRSRV ('LINKSB', K(49), ML*2)
      CALL MDRSRV ('REXTRM', K(50), MR*4)
      CALL MDRSRV ('IHOLDP', K(51), MP*2)
      CALL MDRSRV ('IHOLDL', K(52), ML*2)
      CALL MDRSRV ('IHOLDS', K(53), MS*2)
      CALL MDRSRV ('IHOLDB', K(54), MS*2)
      CALL MDRSRV ('IHOLDR', K(55), MR*2)
      CALL MDRSRV ('IHOLDM', K(56), (MS + MR)*2)
      CALL MDRSRV ('IHOLD1', K(57), MP*2)
      CALL MDRSRV ('IHOLD2', K(58), ML*2)
      CALL MDRSRV ('IHOLD3', K(59), ML*2)
      CALL MDRSRV ('IWTPBF', K(60), MP*3)
      CALL MDRSRV ('IWTLBF', K(61), ML*3)
      CALL MDRSRV ('IWTSBF', K(62), ML*3)
      CALL MDRSRV ('RSIZE', K(63), MR)
      CALL MDRSRV ('IFHOLE', K(64), MR)
      CALL MDRSRV ('NHPR', K(65), MR)
      CALL MDRSRV ('IHLIST', K(66), MR)
      CALL MDRSRV ('IRGFLG', K(67), MR)
      CALL MDSTAT (NERR, MUSED)
      IF (NERR .GT. 0) THEN
         CALL MDEROR (6)
         STOP' '
      END IF

C  THIS ENDS THE SECTION THAT NEEDS TO BE REMOVED IF NOT USING
C  DYNAMIC VARIABLE DIMENSIONING.  AS A REPLACEMENT, THE POINTERS
C  MUST BE HARD WIRED INTO THE PROGRAM.  THIS WOULD BE HANDLED IN THE
C  FOLLOWING PATTERN OF STATEMENTS:
C        K(1) = 1
C        K(2) = K(1) + MP      !NOTE - MP IS THE DIMENSION FOR IPOINT, ETC.
C        K(3) = K(2) + MP*2
C        K(4) = K(3) + MP
C        K(5) = K(4) + ML
C         ....
C        K(67) = K(66) + MR

C-----------------------------------------------------------------------

C  ZERO THE LINK ARRAYS

      CALL LTNEW (MP, IA(K(40)))
      CALL LTNEW (ML, IA(K(41)))
      CALL LTNEW (MS, IA(K(42)))
      CALL LTNEW (MS, IA(K(43)))
      CALL LTNEW (MR, IA(K(44)))
      CALL LTNEW (MS + MR, IA(K(45)))
      CALL LTNEW (MR, IA(K(46)))
      CALL LTNEW (MP, IA(K(47)))
      CALL LTNEW (ML, IA(K(48)))
      CALL LTNEW (ML, IA(K(49)))

C  ENTER FASTQ MAIN OPTION

      IZ = 0
  120 CONTINUE
      IF ((.NOT.BATCH) .AND. (ICOM .GT. JCOM)) THEN
         CALL MESAGE (' ')
         CALL FREFLD (IZ, IZ, 'ENTER OPTION: ', MCOM, IOSTAT, JCOM,
     &      KIN, CIN, IIN, RIN)
         ICOM = 1
      END IF

C  GRAPHICS OPTION - PLOTS FASTQ DATA

      IF ((CIN(ICOM)(1:1) .EQ. 'G') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'g')) THEN
         ICOM = ICOM + 1
         CALL GDATA (MP, ML, MS, MR, MSC, MCOM, ICOM, JCOM, CIN, RIN,
     &      IIN, KIN, IDUMP, N, IA(K(1)), A(K(2)), IA(K(3)), IA(K(4)),
     &      IA(K(5)), IA(K(6)), A(K(7)), IA(K(8)), IA(K(9)), IA(K(10)),
     &      IA(K(11)), IA(K(12)), IA(K(13)), IA(K(14)), IA(K(15)),
     &      IA(K(16)), IA(K(17)), IA(K(18)), IA(K(19)), IA(K(20)),
     &      IA(K(21)), IA(K(22)), IA(K(23)), IA(K(24)), IA(K(25)),
     &      IA(K(40)), IA(K(41)), IA(K(42)), IA(K(43)), IA(K(44)),
     &      IA(K(46)), A(K(50)), A(K(63)), SCHEME, DEFSCH, DEFSIZ,
     &      TITLE, LABP, LABL, LABR, AXISD, LABMD, LABI, LABF, LABPB,
     &      LABLB, LABSBD, LABSC, LABSZ, FULL, IDEV, ALPHA, DEV1,
     &      VAXVMS, VERSN, WROTE, TIME1, HARDPL, BATCH)

C  DELETE OPTION - DELETES FASTQ DATA

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'D') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'd')) THEN
         ICOM = ICOM + 1
         CALL DELFSQ (MP, ML, MS, MR, MSC, MCOM, ICOM, JCOM, CIN, RIN,
     &      IIN, KIN, N, IA(K(3)), IA(K(9)), IA(K(10)), IA(K(12)),
     &      IA(K(13)), IA(K(14)), IA(K(23)), IA(K(24)), IA(K(25)),
     &      IA(K(26)), IA(K(27)), IA(K(28)), IA(K(29)), IA(K(30)),
     &      IA(K(31)), IA(K(32)), IA(K(33)), IA(K(34)), IA(K(35)),
     &      IA(K(36)), IA(K(37)), IA(K(38)), IA(K(40)), IA(K(41)),
     &      IA(K(42)), IA(K(43)), IA(K(44)), IA(K(46)), IA(K(47)),
     &      IA(K(48)), IA(K(49)), IA(K(60)), IA(K(61)), IA(K(62)),
     &      IA(K(64)), IA(K(65)), IA(K(66)), IA(K(67)), NUMBER, DEFSCH,
     &      OPTIM, VAXVMS, WROTE, TIME1, BATCH, VERSN)
         WROTE = .FALSE.

C  FLUSH OPTION - ERASES ALL DATA

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'F') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'f')) THEN
         ICOM = ICOM + 1
         CALL INTRUP ('THIS OPTION ERASES ALL DATA - CONTINUE ANYWAY',
     &      IANS, MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN)
         IF (IANS) THEN
            TITLE = ' '
            DO 130 I = 1, 29
               N(I) = 0
  130       CONTINUE
            NSNAP(1) = 0
            NSNAP(2) = 0
            DO 140 I = 1, MSNAP
               SNAPDX(1, I) = 0.0
               SNAPDX(2, I) = 0.0
  140       CONTINUE
            TBZOOM = .FALSE.
            CALL LTNEW (MP, IA(K(40)))
            CALL LTNEW (ML, IA(K(41)))
            CALL LTNEW (MS, IA(K(42)))
            CALL LTNEW (MS, IA(K(43)))
            CALL LTNEW (MR, IA(K(44)))
            CALL LTNEW (MS + MR, IA(K(45)))
            CALL LTNEW (MR, IA(K(46)))
            CALL LTNEW (MP, IA(K(47)))
            CALL LTNEW (ML, IA(K(48)))
            CALL LTNEW (ML, IA(K(49)))
         END IF

C  MESH OPTION - BEGINS MESH PROCESSING

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'M') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'm')) THEN
         ICOM = ICOM + 1
         CALL MESH (A, IA, MP, ML, MS, MR, MSC, MA, MCOM, ICOM, JCOM,
     &      CIN, RIN, IIN, KIN, IDUMP, N, IA(K(1)), A(K(2)), IA(K(3)),
     &      IA(K(4)), IA(K(5)), IA(K(6)), A(K(7)), IA(K(8)), IA(K(9)),
     &      IA(K(10)), IA(K(11)), IA(K(12)), IA(K(13)), IA(K(14)),
     &      IA(K(15)), IA(K(16)), IA(K(17)), IA(K(18)), IA(K(19)),
     &      IA(K(20)), IA(K(21)), IA(K(22)), IA(K(23)), IA(K(24)),
     &      IA(K(25)), IA(K(26)), IA(K(27)), IA(K(28)), IA(K(29)),
     &      IA(K(30)), IA(K(31)), IA(K(32)), IA(K(33)), IA(K(34)),
     &      IA(K(35)), IA(K(36)), IA(K(37)), IA(K(38)), IA(K(40)),
     &      IA(K(41)), IA(K(42)), IA(K(43)), IA(K(44)), IA(K(45)),
     &      IA(K(46)), IA(K(47)), IA(K(48)), IA(K(49)), IA(K(60)),
     &      IA(K(61)), IA(K(62)), A(K(63)), IA(K(64)), IA(K(65)),
     &      IA(K(66)), IA(K(67)), ISCHM, SCHEME, NUMBER, DEFSCH, DEFSIZ,
     &      TITLE, OPTIM, IDEV, ALPHA, DEV1, THREE, EIGHT, NINE, BATCH,
     &      VAXVMS, VERSN, AXIS, AREACG, LABN, LABE, LABO, LABNB,
     &      LABSB, LABM, LABW, WROTE, TIME1, HARDPL, EXODUSII)

C  SPAWN A PROCESS

      ELSE IF ((CIN(ICOM)(1:2) .EQ. 'SP') .OR.
     &   (CIN(ICOM)(1:2) .EQ. 'sp')) THEN
         ICOM = ICOM + 1
         CALL SPAWN (VAXVMS)

C  STRAIGHTEN OPTION - STRAIGHTEN LINES IN THE X OR Y DIRECTION

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'S') .OR.
     &   (CIN(ICOM)(1:2) .EQ. 's')) THEN
         ICOM = ICOM + 1
         CALL STRAIT (MP, ML, MCOM, ICOM, JCOM, CIN, RIN, IIN, KIN,
     &      IDUMP, N, A(K(2)), IA(K(8)), IA(K(40)), IA(K(41)))

C  TABLET DIGITIZE OPTION - DIGITIZE THE GEOMETRY

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'T') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 't')) THEN
         ICOM = ICOM + 1
         MERGE = .FALSE.
  150    CONTINUE
         CALL TABLET (MP, ML, MS, MR, MSNAP, MCOM, ICOM, JCOM, CIN, RIN,
     &      IIN, KIN, IDUMP, N, IA(K(1)), A(K(2)), IA(K(3)), IA(K(4)),
     &      IA(K(5)), IA(K(6)), A(K(7)), IA(K(8)), IA(K(9)), IA(K(10)),
     &      IA(K(11)), IA(K(12)), IA(K(13)), IA(K(14)), IA(K(21)),
     &      IA(K(22)), IA(K(23)), IA(K(24)), IA(K(25)), IA(K(26)),
     &      IA(K(27)), IA(K(28)), IA(K(29)), IA(K(30)), IA(K(31)),
     &      IA(K(32)), IA(K(33)), IA(K(34)), IA(K(35)), IA(K(36)),
     &      IA(K(37)), IA(K(38)), IA(K(40)), IA(K(41)), IA(K(42)),
     &      IA(K(44)), IA(K(45)), IA(K(47)), IA(K(48)), IA(K(49)),
     &      A(K(50)), IA(K(51)), IA(K(52)), IA(K(53)), IA(K(55)),
     &      IA(K(56)), IA(K(58)), IA(K(59)), IA(K(60)), IA(K(61)),
     &      IA(K(62)), IA(K(67)), TITLE, NOROOM, DRWTAB, XX1, YY1,
     &      SCALE, CT, ST, X1, X2, Y1, Y2, ALPHA, DEV1, SNAP, SNAPDX,
     &      NSNAP, VAXVMS, TBZOOM, AXIST, WROTE, BATCH, VERSN, TIME1)

C  EXTEND THE MEMORY AND CONTINUE IF USING DYNAMIC VARIABLE DIMENSIONING.
C  IN CONVERTING TO NORMAL VARIABLE DIMENSIONING, THE EXTEND MEMORY LINES
C  MUST BE TAKEN OUT, AND AN EXIT OF THE PROGRAM INPUT.  THEN THE
C  PARAMETER STATEMENT CONTAINING MP, ML, MS, AND MR MUST BE INCREASED TO
C  INCREASE DIMESIONING.

         IF (NOROOM) THEN
            MPOLD = MP
            MLOLD = ML
            MSOLD = MS
            MROLD = MR
            MP = NINT(DBLE(MP)*1.5000001)
            ML = NINT(DBLE(ML)*1.5000001)
            MS = NINT(DBLE(MS)*1.5000001)
            MR = NINT(DBLE(MR)*1.5000001)
            CALL MDLONG ('IPOINT', K(1), MP)
            CALL MDLONG ('COOR', K(2), MP*2)
            CALL MDLONG ('IPBOUN', K(3), MP)
            CALL MDLONG ('ILINE', K(4), ML)
            CALL MDLONG ('LTYPE', K(5), ML)
            CALL MDLONG ('NINT', K(6), ML)
            CALL MDLONG ('FACTOR', K(7), ML)
            CALL MDLONG ('LCON', K(8), ML*3)
            CALL MDLONG ('ILBOUN', K(9), ML)
            CALL MDLONG ('ISBOUN', K(10), ML)
            CALL MDLONG ('ISIDE', K(11), MS)
            CALL MDLONG ('NLPS', K(12), MS)
            CALL MDLONG ('IFLINE', K(13), MS)
            CALL MDLONG ('ILLIST', K(14), MS*3)
            CALL MDLONG ('IBARST', K(15), MS)
            CALL MDLONG ('JMAT', K(16), MS)
            CALL MDLONG ('JCENT', K(17), MS)
            CALL MDLONG ('NLPB', K(18), MS)
            CALL MDLONG ('JFLINE', K(19), MS)
            CALL MDLONG ('JLLIST', K(20), MS*3)
            CALL MDLONG ('IREGN', K(21), MR)
            CALL MDLONG ('IMAT', K(22), MR)
            CALL MDLONG ('NSPR', K(23), MR)
            CALL MDLONG ('IFSIDE', K(24), MR)
            CALL MDLONG ('ISLIST', K(25), MR*4)
            CALL MDLONG ('IRPB', K(26), MR)
            CALL MDLONG ('IPBF', K(27), MP)
            CALL MDLONG ('NPPF', K(28), MP)
            CALL MDLONG ('IFPB', K(29), MP)
            CALL MDLONG ('LISTPB', K(30), MP*2)
            CALL MDLONG ('ILBF', K(31), ML)
            CALL MDLONG ('NLPF', K(32), ML)
            CALL MDLONG ('IFLB', K(33), ML)
            CALL MDLONG ('LISTLB', K(34), ML*2)
            CALL MDLONG ('ISBF', K(35), ML)
            CALL MDLONG ('NSPF', K(36), ML)
            CALL MDLONG ('IFSB', K(37), ML)
            CALL MDLONG ('LISTSB', K(38), ML*2)
            CALL MDLONG ('ATTRIB', K(39), MA*(MR + MS))
            CALL MDLONG ('LINKP', K(40), MP*2)
            CALL MDLONG ('LINKL', K(41), ML*2)
            CALL MDLONG ('LINKS', K(42), MS*2)
            CALL MDLONG ('LINKB', K(43), MS*2)
            CALL MDLONG ('LINKR', K(44), MR*2)
            CALL MDLONG ('LINKM', K(45), (MS + MR)*2)
            CALL MDLONG ('LINKSC', K(46), MR*2)
            CALL MDLONG ('LINKPB', K(47), MP*2)
            CALL MDLONG ('LINKLB', K(48), ML*2)
            CALL MDLONG ('LINKSB', K(49), ML*2)
            CALL MDLONG ('REXTRM', K(50), MR*4)
            CALL MDLONG ('IHOLDP', K(51), MP*2)
            CALL MDLONG ('IHOLDL', K(52), ML*2)
            CALL MDLONG ('IHOLDS', K(53), MS*2)
            CALL MDLONG ('IHOLDB', K(54), MS*2)
            CALL MDLONG ('IHOLDR', K(55), MR*2)
            CALL MDLONG ('IHOLDM', K(56), (MS + MR)*2)
            CALL MDLONG ('IHOLD1', K(57), MP*2)
            CALL MDLONG ('IHOLD2', K(58), ML*2)
            CALL MDLONG ('IHOLD3', K(59), ML*2)
            CALL MDLONG ('IWTPBF', K(60), MP*3)
            CALL MDLONG ('IWTLBF', K(61), ML*3)
            CALL MDLONG ('IWTSBF', K(62), ML*3)
            CALL MDLONG ('RSIZE', K(63), MR)
            CALL MDLONG ('IFHOLE', K(64), MR)
            CALL MDLONG ('NHPR', K(65), MR)
            CALL MDLONG ('IHLIST', K(66), MR)
            CALL MDLONG ('IRGFLG', K(67), MR)
            CALL MDSTAT (NERR, MUSED)
            IF (NERR .GT. 0) THEN
               CALL MDEROR (6)
               STOP' '
            END IF

C  RESORT THE LINK ARRAYS

            CALL LTNEW (ML, IA(K(51)))
            CALL LTADD (ML, MLOLD, N(1), IA(K(40)), IA(K(51)))
            CALL LTNEW (ML, IA(K(51)))
            CALL LTNEW (ML, IA(K(52)))
            CALL LTADD (ML, MLOLD, N(2), IA(K(41)), IA(K(52)))
            CALL LTNEW (ML, IA(K(52)))
            CALL LTNEW (MS, IA(K(53)))
            CALL LTADD (MS, MSOLD, N(3), IA(K(42)), IA(K(53)))
            CALL LTNEW (MS, IA(K(53)))
            CALL LTNEW (MS, IA(K(54)))
            CALL LTADD (MS, MSOLD, N(5), IA(K(43)), IA(K(54)))
            CALL LTNEW (MS, IA(K(54)))
            CALL LTNEW (MR, IA(K(55)))
            CALL LTADD (MR, MROLD, N(7), IA(K(44)), IA(K(55)))
            CALL LTNEW (MS + MR, IA(K(56)))
            CALL LTADD (MS + MR, MSOLD + MROLD, MSOLD + MROLD,
     &         IA(K(45)), IA(K(56)))
            CALL LTNEW (MS + MR, IA(K(56)))
            CALL LTNEW (MR, IA(K(55)))
            CALL LTADD (MR, MROLD, N(8), IA(K(46)), IA(K(55)))
            CALL LTNEW (MR, IA(K(55)))
            CALL LTNEW (MP, IA(K(57)))
            CALL LTADD (MP, MPOLD, N(11), IA(K(47)), IA(K(57)))
            CALL LTNEW (MP, IA(K(57)))
            CALL LTNEW (ML, IA(K(58)))
            CALL LTADD (ML, MLOLD, N(13), IA(K(48)), IA(K(58)))
            CALL LTNEW (ML, IA(K(58)))
            CALL LTNEW (ML, IA(K(59)))
            CALL LTADD (ML, MLOLD, N(15), IA(K(49)), IA(K(59)))
            CALL MESAGE('DIGITIZATION CAN NOW BE CONTINUED')
            GO TO 150
         END IF
         WROTE = .FALSE.

C  KEY-IN OPTION - TYPE IN THE DATA FROM THE KEYBOARD

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'K') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'k')) THEN
         ICOM = ICOM + 1
         MERGE = .FALSE.
         WROTE = .FALSE.
  160    CONTINUE
         CALL KEYIN (MP, ML, MS, MR, MSC, MA, MCOM, ICOM, JCOM, CIN,
     &      RIN, IIN, KIN, IDUMP, N, IA(K(1)), A(K(2)), IA(K(3)),
     &      IA(K(4)), IA(K(5)), IA(K(6)), A(K(7)), IA(K(8)), IA(K(9)),
     &      IA(K(10)), IA(K(11)), IA(K(12)), IA(K(13)), IA(K(14)),
     &      IA(K(15)), IA(K(16)), IA(K(17)), IA(K(18)), IA(K(19)),
     &      IA(K(20)), IA(K(21)), IA(K(22)), IA(K(23)), IA(K(24)),
     &      IA(K(25)), IA(K(26)), IA(K(27)), IA(K(28)), IA(K(29)),
     &      IA(K(30)), IA(K(31)), IA(K(32)), IA(K(33)), IA(K(34)),
     &      IA(K(35)), IA(K(36)), IA(K(37)), IA(K(38)), IA(K(40)),
     &      IA(K(41)), IA(K(42)), IA(K(43)), IA(K(44)), IA(K(45)),
     &      IA(K(46)), IA(K(47)), IA(K(48)), IA(K(49)), IA(K(51)),
     &      IA(K(52)), IA(K(53)), IA(K(54)), IA(K(55)), IA(K(56)),
     &      IA(K(57)), IA(K(58)), IA(K(59)), IA(K(60)), IA(K(61)),
     &      IA(K(62)), A(K(63)), IA(K(64)), IA(K(65)), IA(K(66)),
     &      IA(K(67)), ISCHM, SCHEME, NUMBER, DEFSCH, DEFSIZ, TITLE,
     &      OPTIM, THREE, EIGHT, NINE, NOROOM, VAXVMS, WROTE, TIME1,
     &      VERSN, BATCH)
         IF (NOROOM) THEN
            MPOLD = MP
            MLOLD = ML
            MSOLD = MS
            MROLD = MR
            MP = NINT(DBLE(MP)*1.5000001)
            ML = NINT(DBLE(ML)*1.5000001)
            MS = NINT(DBLE(MS)*1.5000001)
            MR = NINT(DBLE(MR)*1.5000001)
            CALL MDLONG ('IPOINT', K(1), MP)
            CALL MDLONG ('COOR', K(2), MP*2)
            CALL MDLONG ('IPBOUN', K(3), MP)
            CALL MDLONG ('ILINE', K(4), ML)
            CALL MDLONG ('LTYPE', K(5), ML)
            CALL MDLONG ('NINT', K(6), ML)
            CALL MDLONG ('FACTOR', K(7), ML)
            CALL MDLONG ('LCON', K(8), ML*3)
            CALL MDLONG ('ILBOUN', K(9), ML)
            CALL MDLONG ('ISBOUN', K(10), ML)
            CALL MDLONG ('ISIDE', K(11), MS)
            CALL MDLONG ('NLPS', K(12), MS)
            CALL MDLONG ('IFLINE', K(13), MS)
            CALL MDLONG ('ILLIST', K(14), MS*3)
            CALL MDLONG ('IBARST', K(15), MS)
            CALL MDLONG ('JMAT', K(16), MS)
            CALL MDLONG ('JCENT', K(17), MS)
            CALL MDLONG ('NLPB', K(18), MS)
            CALL MDLONG ('JFLINE', K(19), MS)
            CALL MDLONG ('JLLIST', K(20), MS*3)
            CALL MDLONG ('IREGN', K(21), MR)
            CALL MDLONG ('IMAT', K(22), MR)
            CALL MDLONG ('NSPR', K(23), MR)
            CALL MDLONG ('IFSIDE', K(24), MR)
            CALL MDLONG ('ISLIST', K(25), MR*4)
            CALL MDLONG ('IRPB', K(26), MR)
            CALL MDLONG ('IPBF', K(27), MP)
            CALL MDLONG ('NPPF', K(28), MP)
            CALL MDLONG ('IFPB', K(29), MP)
            CALL MDLONG ('LISTPB', K(30), MP*2)
            CALL MDLONG ('ILBF', K(31), ML)
            CALL MDLONG ('NLPF', K(32), ML)
            CALL MDLONG ('IFLB', K(33), ML)
            CALL MDLONG ('LISTLB', K(34), ML*2)
            CALL MDLONG ('ISBF', K(35), ML)
            CALL MDLONG ('NSPF', K(36), ML)
            CALL MDLONG ('IFSB', K(37), ML)
            CALL MDLONG ('LISTSB', K(38), ML*2)
            CALL MDLONG ('ATTRIB', K(39), MA*(MR + MS))
            CALL MDLONG ('LINKP', K(40), MP*2)
            CALL MDLONG ('LINKL', K(41), ML*2)
            CALL MDLONG ('LINKS', K(42), MS*2)
            CALL MDLONG ('LINKB', K(43), MS*2)
            CALL MDLONG ('LINKR', K(44), MR*2)
            CALL MDLONG ('LINKM', K(45), (MS + MR)*2)
            CALL MDLONG ('LINKSC', K(46), MR*2)
            CALL MDLONG ('LINKPB', K(47), MP*2)
            CALL MDLONG ('LINKLB', K(48), ML*2)
            CALL MDLONG ('LINKSB', K(49), ML*2)
            CALL MDLONG ('REXTRM', K(50), MR*4)
            CALL MDLONG ('IHOLDP', K(51), MP*2)
            CALL MDLONG ('IHOLDL', K(52), ML*2)
            CALL MDLONG ('IHOLDS', K(53), MS*2)
            CALL MDLONG ('IHOLDB', K(54), MS*2)
            CALL MDLONG ('IHOLDR', K(55), MR*2)
            CALL MDLONG ('IHOLDM', K(56), (MS + MR)*2)
            CALL MDLONG ('IHOLD1', K(57), MP*2)
            CALL MDLONG ('IHOLD2', K(58), ML*2)
            CALL MDLONG ('IHOLD3', K(59), ML*2)
            CALL MDLONG ('IWTPBF', K(60), MP*3)
            CALL MDLONG ('IWTLBF', K(61), ML*3)
            CALL MDLONG ('IWTSBF', K(62), ML*3)
            CALL MDLONG ('RSIZE', K(63), MR)
            CALL MDLONG ('IFHOLE', K(64), MR)
            CALL MDLONG ('NHPR', K(65), MR)
            CALL MDLONG ('IHLIST', K(66), MR)
            CALL MDLONG ('IRGFLG', K(67), MR)
            CALL MDSTAT (NERR, MUSED)
            IF (NERR .GT. 0) THEN
               CALL MDEROR (6)
               STOP' '
            END IF

C  RESORT THE LINK ARRAYS

            CALL LTNEW (ML, IA(K(51)))
            CALL LTADD (ML, MLOLD, N(1), IA(K(40)), IA(K(51)))
            CALL LTNEW (ML, IA(K(51)))
            CALL LTNEW (ML, IA(K(52)))
            CALL LTADD (ML, MLOLD, N(2), IA(K(41)), IA(K(52)))
            CALL LTNEW (ML, IA(K(52)))
            CALL LTNEW (MS, IA(K(53)))
            CALL LTADD (MS, MSOLD, N(3), IA(K(42)), IA(K(53)))
            CALL LTNEW (MS, IA(K(53)))
            CALL LTNEW (MS, IA(K(54)))
            CALL LTADD (MS, MSOLD, N(5), IA(K(43)), IA(K(54)))
            CALL LTNEW (MS, IA(K(54)))
            CALL LTNEW (MR, IA(K(55)))
            CALL LTADD (MR, MROLD, N(7), IA(K(44)), IA(K(55)))
            CALL LTNEW (MS + MR, IA(K(56)))
            CALL LTADD (MS + MR, MSOLD + MROLD, MSOLD + MROLD,
     &         IA(K(45)), IA(K(56)))
            CALL LTNEW (MS + MR, IA(K(56)))
            CALL LTNEW (MR, IA(K(55)))
            CALL LTADD (MR, MROLD, N(8), IA(K(46)), IA(K(55)))
            CALL LTNEW (MR, IA(K(55)))
            CALL LTNEW (MP, IA(K(57)))
            CALL LTADD (MP, MPOLD, N(11), IA(K(47)), IA(K(57)))
            CALL LTNEW (MP, IA(K(57)))
            CALL LTNEW (ML, IA(K(58)))
            CALL LTADD (ML, MLOLD, N(13), IA(K(48)), IA(K(58)))
            CALL LTNEW (ML, IA(K(58)))
            CALL LTNEW (ML, IA(K(59)))
            CALL LTADD (ML, MLOLD, N(15), IA(K(49)), IA(K(59)))
            CALL MESAGE ('KEYIN OPTION CAN NOW BE CONTINUED')
            GO TO 160
         END IF

C  LIST OPTION - LISTS FASTQ DATA

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'L') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'l')) THEN
         ICOM = ICOM + 1
         CALL LIST (MP, ML, MS, MR, MSC, MCOM, ICOM, JCOM, CIN, RIN,
     &      IIN, KIN, N, IA(K(1)), A(K(2)), IA(K(3)), IA(K(4)),
     &      IA(K(5)), IA(K(6)), A(K(7)), IA(K(8)), IA(K(9)), IA(K(10)),
     &      IA(K(11)), IA(K(12)), IA(K(13)), IA(K(14)), IA(K(15)),
     &      IA(K(16)), IA(K(17)), IA(K(18)), IA(K(19)), IA(K(20)),
     &      IA(K(21)), IA(K(22)), IA(K(23)), IA(K(24)), IA(K(25)),
     &      IA(K(26)), IA(K(27)), IA(K(28)), IA(K(29)), IA(K(30)),
     &      IA(K(31)), IA(K(32)), IA(K(33)), IA(K(34)), IA(K(35)),
     &      IA(K(36)), IA(K(37)), IA(K(38)), IA(K(40)), IA(K(41)),
     &      IA(K(42)), IA(K(43)), IA(K(44)), IA(K(46)), IA(K(47)),
     &      IA(K(48)), IA(K(49)), IA(K(60)), IA(K(61)), IA(K(62)),
     &      A(K(63)), IA(K(64)), IA(K(65)), IA(K(66)), IA(K(67)), ISCHM,
     &      SCHEME, NUMBER, DEFSCH, DEFSIZ, TITLE, OPTIM, THREE, EIGHT,
     &      NINE, VAXVMS, WROTE, TIME1, VERSN, BATCH)

C  READ OPTION - READS FASTQ DATA

      ELSE IF (((CIN(ICOM)(1:1) .EQ. 'R') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'r')) .AND.
     &   (CIN(ICOM)(2:2).NE.'P') .AND. (CIN(ICOM)(2:2).NE.'p') .AND.
     &   (CIN(ICOM)(2:2).NE.'W') .AND. (CIN(ICOM)(2:2).NE.'w')) THEN
         ICOM = ICOM + 1
         IF ((N(1) .GT. 0) .OR. (N(2) .GT. 0)) THEN

C  CHECK TO SEE IF A FASTQ DATA MERGING IS DESIRED

            CALL INTRUP ('MERGE FILE WITH EXISTING DATA', MERGE, MCOM,
     &         ICOM, JCOM, CIN, IIN, RIN, KIN)
            IF (MERGE) THEN
               CALL LTNEW (MP, IA(K(51)))
               CALL LTNEW (ML, IA(K(52)))
               CALL LTNEW (MS, IA(K(53)))
               CALL LTNEW (MS, IA(K(54)))
               CALL LTNEW (MR, IA(K(55)))
               CALL LTNEW (MS + MR, IA(K(56)))
               CALL LTNEW (MP, IA(K(57)))
               CALL LTNEW (ML, IA(K(58)))
               CALL LTNEW (ML, IA(K(59)))
            ELSE
               IF (.NOT.WROTE) THEN
                  CALL MESAGE('CONTINUATION WILL OVERWRITE OLD DATA')
                  CALL INTRUP('DO YOU WISH TO CONTINUE', IANS, MCOM,
     &               ICOM, JCOM, CIN, IIN, RIN, KIN)
                  IF (.NOT.IANS) GO TO 120
               END IF
               DO 170 I = 1, 29
                  N(I) = 0
  170          CONTINUE
               CALL LTNEW (MP, IA(K(40)))
               CALL LTNEW (ML, IA(K(41)))
               CALL LTNEW (MS, IA(K(42)))
               CALL LTNEW (MS, IA(K(43)))
               CALL LTNEW (MR, IA(K(44)))
               CALL LTNEW (MS + MR, IA(K(45)))
               CALL LTNEW (MR, IA(K(46)))
               CALL LTNEW (MP, IA(K(47)))
               CALL LTNEW (ML, IA(K(48)))
               CALL LTNEW (ML, IA(K(49)))
            END IF
         END IF

         IUNIT = 1
         ITRY = 0
  180    CONTINUE
         IF (ITRY .GT. 0) THEN
            CALL STRLNG (FNAME, LEN)
            IF (FNAME(1:6).NE.'FOR001') WRITE (*, 10000) FNAME(1:LEN)
         END IF
         ITRY = ITRY + 1
         IF (((ITRY .LE. 3) .AND. (.NOT.BATCH)) .OR.
     &      ((BATCH) .AND. (ITRY .LE. 1)) .OR.
     &      ((START) .AND. (.NOT.BATCH))) THEN
            IUNIT = 1
            IF (BATCH) THEN
               IDUMP = 6
               MERGE = .FALSE.
               CIN(1) = 'MESH'
               CIN(2) = 'PROC'
               CIN(3) = 'WRITE'
               CIN(4) = 'EXIT'
               ICOM = 1
               JCOM = 4
               CALL EXNAME (IUNIT, FNAME, LN)
               OPEN (UNIT = IUNIT, FILE = FNAME(1:LN), STATUS = 'OLD',
     &            ERR = 180)
            ELSE IF (START) THEN
               START = .FALSE.
               CALL EXNAME (IUNIT, FNAME, LN)
               IF ((.NOT. VAXVMS) .AND. (FNAME .EQ. 'tty')) GO TO 120
               OPEN (UNIT = IUNIT, FILE = FNAME(1:LN), STATUS = 'OLD',
     &            ERR = 120)
            ELSE
               IF (ICOM .LE. JCOM) THEN
                  FNAME = CIN(ICOM)
                  ICOM = ICOM + 1
               ELSE
                  CALL INQSTR ('INPUT FILE: ', FNAME)
               END IF
               OPEN (UNIT = IUNIT, FILE = FNAME, STATUS = 'OLD',
     &            ERR = 180)
               IDUMP = 0
            END IF
  190       CONTINUE
            REWIND IUNIT
            CALL RDFSQ (MP, ML, MS, MR, MSNAP, MSC, MA, IUNIT, IDUMP, N,
     &         IA(K(1)), A(K(2)), IA(K(3)), IA(K(4)), IA(K(5)),
     &         IA(K(6)), A(K(7)), IA(K(8)), IA(K(9)), IA(K(10)),
     &         IA(K(11)), IA(K(12)), IA(K(13)), IA(K(14)), IA(K(15)),
     &         IA(K(16)), IA(K(17)), IA(K(18)), IA(K(19)), IA(K(20)),
     &         IA(K(21)), IA(K(22)), IA(K(23)), IA(K(24)), IA(K(25)),
     &         IA(K(26)), IA(K(27)), IA(K(28)), IA(K(29)), IA(K(30)),
     &         IA(K(31)), IA(K(32)), IA(K(33)), IA(K(34)), IA(K(35)),
     &         IA(K(36)), IA(K(37)), IA(K(38)), A(K(39)), IA(K(40)),
     &         IA(K(41)), IA(K(42)), IA(K(43)), IA(K(44)), IA(K(45)),
     &         IA(K(46)), IA(K(47)), IA(K(48)), IA(K(49)), IA(K(51)),
     &         IA(K(52)), IA(K(53)), IA(K(54)), IA(K(55)), IA(K(56)),
     &         IA(K(57)), IA(K(58)), IA(K(59)), IA(K(60)), IA(K(61)),
     &         IA(K(62)), A(K(63)), IA(K(64)), IA(K(65)), IA(K(66)),
     &         IA(K(67)), ISCHM, SCHEME, NUMBER, DEFSCH, DEFSIZ, TITLE,
     &         OPTIM, MERGE, THREE, EIGHT, NINE, SNAP, SNAPDX, NSNAP,
     &         RATIO, NOROOM, EXODUSII)
            IF (NOROOM) THEN
               MPOLD = MP
               MLOLD = ML
               MSOLD = MS
               MROLD = MR
               MP = NINT(DBLE(MP)*RATIO)
               ML = NINT(DBLE(ML)*RATIO)
               MS = NINT(DBLE(MS)*RATIO)
               MR = NINT(DBLE(MR)*RATIO)
               CALL MDLONG ('IPOINT', K(1), MP)
               CALL MDLONG ('COOR', K(2), MP*2)
               CALL MDLONG ('IPBOUN', K(3), MP)
               CALL MDLONG ('ILINE', K(4), ML)
               CALL MDLONG ('LTYPE', K(5), ML)
               CALL MDLONG ('NINT', K(6), ML)
               CALL MDLONG ('FACTOR', K(7), ML)
               CALL MDLONG ('LCON', K(8), ML*3)
               CALL MDLONG ('ILBOUN', K(9), ML)
               CALL MDLONG ('ISBOUN', K(10), ML)
               CALL MDLONG ('ISIDE', K(11), MS)
               CALL MDLONG ('NLPS', K(12), MS)
               CALL MDLONG ('IFLINE', K(13), MS)
               CALL MDLONG ('ILLIST', K(14), MS*3)
               CALL MDLONG ('IBARST', K(15), MS)
               CALL MDLONG ('JMAT', K(16), MS)
               CALL MDLONG ('JCENT', K(17), MS)
               CALL MDLONG ('NLPB', K(18), MS)
               CALL MDLONG ('JFLINE', K(19), MS)
               CALL MDLONG ('JLLIST', K(20), MS*3)
               CALL MDLONG ('IREGN', K(21), MR)
               CALL MDLONG ('IMAT', K(22), MR)
               CALL MDLONG ('NSPR', K(23), MR)
               CALL MDLONG ('IFSIDE', K(24), MR)
               CALL MDLONG ('ISLIST', K(25), MR*4)
               CALL MDLONG ('IRPB', K(26), MR)
               CALL MDLONG ('IPBF', K(27), MP)
               CALL MDLONG ('NPPF', K(28), MP)
               CALL MDLONG ('IFPB', K(29), MP)
               CALL MDLONG ('LISTPB', K(30), MP*2)
               CALL MDLONG ('ILBF', K(31), ML)
               CALL MDLONG ('NLPF', K(32), ML)
               CALL MDLONG ('IFLB', K(33), ML)
               CALL MDLONG ('LISTLB', K(34), ML*2)
               CALL MDLONG ('ISBF', K(35), ML)
               CALL MDLONG ('NSPF', K(36), ML)
               CALL MDLONG ('IFSB', K(37), ML)
               CALL MDLONG ('LISTSB', K(38), ML*2)
               CALL MDLONG ('ATTRIB', K(39), MA*(MR + MS))
               CALL MDLONG ('LINKP', K(40), MP*2)
               CALL MDLONG ('LINKL', K(41), ML*2)
               CALL MDLONG ('LINKS', K(42), MS*2)
               CALL MDLONG ('LINKB', K(43), MS*2)
               CALL MDLONG ('LINKR', K(44), MR*2)
               CALL MDLONG ('LINKM', K(45), (MS + MR)*2)
               CALL MDLONG ('LINKSC', K(46), MR*2)
               CALL MDLONG ('LINKPB', K(47), MP*2)
               CALL MDLONG ('LINKLB', K(48), ML*2)
               CALL MDLONG ('LINKSB', K(49), ML*2)
               CALL MDLONG ('REXTRM', K(50), MR*4)
               CALL MDLONG ('IHOLDP', K(51), MP*2)
               CALL MDLONG ('IHOLDL', K(52), ML*2)
               CALL MDLONG ('IHOLDS', K(53), MS*2)
               CALL MDLONG ('IHOLDB', K(54), MS*2)
               CALL MDLONG ('IHOLDR', K(55), MR*2)
               CALL MDLONG ('IHOLDM', K(56), (MS + MR)*2)
               CALL MDLONG ('IHOLD1', K(57), MP*2)
               CALL MDLONG ('IHOLD2', K(58), ML*2)
               CALL MDLONG ('IHOLD3', K(59), ML*2)
               CALL MDLONG ('IWTPBF', K(60), MP*3)
               CALL MDLONG ('IWTLBF', K(61), ML*3)
               CALL MDLONG ('IWTSBF', K(62), ML*3)
               CALL MDLONG ('RSIZE', K(63), MR)
               CALL MDLONG ('IFHOLE', K(64), MR)
               CALL MDLONG ('NHPR', K(65), MR)
               CALL MDLONG ('IHLIST', K(66), MR)
               CALL MDLONG ('IRGFLG', K(67), MR)
               CALL MDSTAT (NERR, MUSED)
               IF (NERR .GT. 0) THEN
                  CALL MDEROR (6)
                  STOP' '
               END IF

C  RESORT THE LINK ARRAYS

               CALL LTNEW (ML, IA(K(51)))
               CALL LTADD (ML, MLOLD, N(1), IA(K(40)), IA(K(51)))
               CALL LTNEW (ML, IA(K(51)))
               CALL LTNEW (ML, IA(K(52)))
               CALL LTADD (ML, MLOLD, N(2), IA(K(41)), IA(K(52)))
               CALL LTNEW (ML, IA(K(52)))
               CALL LTNEW (MS, IA(K(53)))
               CALL LTADD (MS, MSOLD, N(3), IA(K(42)), IA(K(53)))
               CALL LTNEW (MS, IA(K(53)))
               CALL LTNEW (MS, IA(K(54)))
               CALL LTADD (MS, MSOLD, N(5), IA(K(43)), IA(K(54)))
               CALL LTNEW (MS, IA(K(54)))
               CALL LTNEW (MR, IA(K(55)))
               CALL LTADD (MR, MROLD, N(7), IA(K(44)), IA(K(55)))
               CALL LTNEW (MS + MR, IA(K(56)))
               CALL LTADD (MS + MR, MSOLD + MROLD, MSOLD + MROLD,
     &            IA(K(45)), IA(K(56)))
               CALL LTNEW (MS + MR, IA(K(56)))
               CALL LTNEW (MR, IA(K(55)))
               CALL LTADD (MR, MROLD, N(8), IA(K(46)), IA(K(55)))
               CALL LTNEW (MR, IA(K(55)))
               CALL LTNEW (MP, IA(K(57)))
               CALL LTADD (MP, MPOLD, N(11), IA(K(47)), IA(K(57)))
               CALL LTNEW (MP, IA(K(57)))
               CALL LTNEW (ML, IA(K(58)))
               CALL LTADD (ML, MLOLD, N(13), IA(K(48)), IA(K(58)))
               CALL LTNEW (ML, IA(K(58)))
               CALL LTNEW (ML, IA(K(59)))
               CALL LTADD (ML, MLOLD, N(15), IA(K(49)), IA(K(59)))
               CALL MESAGE('FILE WILL NOW BE READ AGAIN AS NEW INPUT')
               GO TO 190
            END IF
            TBZOOM = .FALSE.
            REWIND IUNIT
            CLOSE (IUNIT)
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
         GO TO 120

C  WRITE OPTION - WRITES A FASTQ DATA FILE

      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'W') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'w') .OR.
     &   (CIN(ICOM)(1:2) .EQ. 'BW') .OR.
     &   (CIN(ICOM)(1:2) .EQ. 'bw') .OR.
     &   (CIN(ICOM)(1:2) .EQ. 'RW') .OR.
     &   (CIN(ICOM)(1:2) .EQ. 'rw')) THEN
         IF ((CIN(ICOM)(1:2) .EQ. 'RW') .OR.
     &      (CIN(ICOM)(1:2) .EQ. 'rw')) THEN
            REGWRT = .TRUE.
            BARWRT = .FALSE.
         ELSEIF ((CIN(ICOM)(1:2) .EQ. 'BW') .OR.
     &      (CIN(ICOM)(1:2) .EQ. 'bw')) THEN
            REGWRT = .FALSE.
            BARWRT = .TRUE.
         ELSE
            REGWRT = .FALSE.
            BARWRT = .FALSE.
         ENDIF
         ICOM = ICOM + 1
         IUNIT = 1
  200    CONTINUE
         IF (ICOM .LE. JCOM) THEN
            FNAME = CIN(ICOM)
            ICOM = ICOM + 1
         ELSE
            CALL INQSTR ('FASTQ DATA FILE: ', FNAME)
         END IF
         OPEN (UNIT = IUNIT, FILE = FNAME, STATUS = 'NEW', ERR = 200)
         CALL WRFSQ (IUNIT, MP, ML, MS, MR, MSNAP, MSC, MCOM, ICOM,
     &      JCOM, CIN, RIN, IIN, KIN, N, IA(K(1)), A(K(2)), IA(K(3)),
     &      IA(K(4)), IA(K(5)), IA(K(6)), A(K(7)), IA(K(8)), IA(K(9)),
     &      IA(K(10)), IA(K(11)), IA(K(12)), IA(K(13)), IA(K(14)),
     &      IA(K(15)), IA(K(16)), IA(K(17)), IA(K(18)), IA(K(19)),
     &      IA(K(20)), IA(K(21)), IA(K(22)), IA(K(23)), IA(K(24)),
     &      IA(K(25)), IA(K(26)), IA(K(27)), IA(K(28)), IA(K(29)),
     &      IA(K(30)), IA(K(31)), IA(K(32)), IA(K(33)), IA(K(34)),
     &      IA(K(35)), IA(K(36)), IA(K(37)), IA(K(38)), IA(K(40)),
     &      IA(K(41)), IA(K(42)), IA(K(43)), IA(K(44)), IA(K(46)),
     &      IA(K(47)), IA(K(48)), IA(K(49)), IA(K(60)), IA(K(61)),
     &      IA(K(62)), A(K(63)), IA(K(64)), IA(K(65)), IA(K(66)),
     &      IA(K(67)), ISCHM, SCHEME, NUMBER, DEFSCH, DEFSIZ, TITLE,
     &      OPTIM, THREE, EIGHT, NINE, SNAP, SNAPDX, NSNAP, REGWRT,
     &      BARWRT)
         WROTE = .TRUE.
         CLOSE (IUNIT)

C  GET THE APPROPRIATE HELP MESAGE

      ELSE
         ICOM = ICOM + 1
         CALL HELP_FQ (1)
      END IF
      GO TO 120

10000 FORMAT (' ', 'ERROR OPENING FILE: ', A)
      END
