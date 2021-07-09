C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE COMMAND (A, IA, TITLE, TIME, ITMSEL, MAT, DISP,
     *   CRD, LINK, DENSTY, WAVE, ISEVOK,
     *   NAMEGL, NAMENV, NAMEEL,
     *   NQAREC, QAREC, NINFO, INFREC, DBNAME)

C        READ AND INTERPRET ALL INPUT DATA

      include 'exodusII.inc'
      CHARACTER*80 TITLE, COMMENT
      INTEGER IA(*)
      DIMENSION A(*), TIME(*), MAT(6,*), DISP(NUMNP,*), CRD(NUMNP,*),
     *   LINK(*), DENSTY(*), WAVE(*)
      INTEGER ITMSEL(*), ISEVOK(*)
      CHARACTER*(MXSTLN) NAMEGL(*), NAMENV(*), NAMEEL(*)
      CHARACTER*(MXLNLN) INFREC(*)
      CHARACTER*(MXSTLN) QAREC(4,*)
      CHARACTER*(*) DBNAME

      include 'nu_cvty.blk'
      include 'nu_ptim.blk'
      include 'nu_mass.blk'
      include 'nu_logs.blk'
      include 'nu_numg.blk'
      include 'nu_varcnt.blk'
      include 'nu_cav.blk'
      include 'nu_nset.blk'
      include 'nu_io.blk'

      DIMENSION TRANGE(3), CENT(3)
      PARAMETER (MXNAM = 20)
      DIMENSION KV(MXNAM), RV(MXNAM), IVAL(MXNAM)
      CHARACTER*32  CV(MXNAM), NAME, CTMP
      CHARACTER*16 LABEL(32), TYPE, PROMPT
      CHARACTER*8 GMTHD, SORTYP, LISTYP, OPT
      LOGICAL FIRST, EXOSAV, ALLTIM, ISON, SORUP, ISHELP, LTMP,
     *   LTMP2, LTMP3, LTMP4, CENTER
      LOGICAL FFNUMB, FFMATC, HELP, MATSTR, FFEXST

      CHARACTER*8 CMDTBL(37), SORTBL(12), LISTBL(20)
      SAVE CMDTBL, SORTBL, LISTBL
C      --CMDTBL - the valid commands table
C      --SORTBL - the valid sort options table
C      --LISTBL - the valid list options table

C   --Command table follows.  Remember to change the dimensioned size when
C   --changing the table.

      DATA CMDTBL /
     *   'AXISYMME', 'CAVITY  ', 'COMMENT ', 'DELTIME ', 'DENSITY ',
     *   'EXIT    ', 'END     ', 'EXODUS  ', 'GAP     ', 'HELP    ',
     *   'LIMITS  ', 'LIST    ', 'LOCATE  ', 'OVERLAP ', 'PLANE   ',
     *   'PLANAR  ', 'PROPERTI', 'SORT    ', 'TIMES   ', 'TMAX    ',
     *   'TMIN    ', 'MASS    ', 'PRINT   ', 'ECHO    ', 'SELECT  ',
     *   'WAVESPEE', 'TIMESTEP', 'ALLTIMES', 'NINTV   ', 'ZINTV   ',
     *   'SUM     ', 'AVERAGE ', 'CONDITIO', 'ESUM    ', 'EAVERAGE',
     *   'QUIT    ', '        '/

      DATA SORTBL /
     *   'X       ', 'Y       ', 'Z       ', 'T       ', 'DISTANCE',
     *   'RADIAL  ', 'PARAMETR', 'ANGLE   ', 'THETA   ', 'PHI     ',
     *   'NONE    ', '        '/

      DATA LISTBL /
     *   'SSETS   ', 'SIDESETS', 'NSETS   ', 'NODESETS', 'VARS    ',
     *   'VARIABLE', 'BLOCKS  ', 'MATERIAL', 'TIMES   ', 'STEPS   ',
     *   'COMMANDS', 'SORT    ', 'STIMES  ', 'NAMES   ', 'SELECTED',
     *   'INFORMAT', 'QA      ', 'VOLUME  ', 'NODALVOL', '        '/

      DATA FIRST /.TRUE./, PROMPT /' NUMBERS> '/
      DATA GMTHD /'DISTANCE'/, SORTYP /'NONE'/, SORUP /.TRUE./

      AXI = .TRUE.
      EXOSAV = EXODUS
      IF (FIRST) THEN
         FIRST = .FALSE.
         CALL HEADER (NDIM, TITLE, NUMEL, NUMNP, AXI)
      END IF

      DO 10 IBLK = 1, NELBLK
         MAT(5,IBLK) = 1
 10   CONTINUE
      if (nstep .gt. 0) then
         CALL INILOG (NSTEP, .TRUE., ITMSEL)
         TMIN  = TIME(1)
         TMAX  = TIME(NSTEP)
      else
         tmin = 0.0
         tmax = 0.0
      end if
      STMIN = TMIN
      STMAX = TMAX
      LSTSEL= NSTEP
      IOMIN = 6
      IOMAX = 7
      NQUAD = 1

C     ... GET SOME SCRATCH SPACE

      CALL MDRSRV ('SCRTCH', ISCR, NDIM*NUMNP)
      CALL MDRSRV ('SCRTC2', ISCR2, NDIM*NUMNP)
      CALL MDRSRV ('SORTMP', ISMP, MAX(NUMEL, NUMNP) )
      CALL MDRSRV ('NODSEL', INDSEL, NUMNP)
      CALL MDRSRV ('ELMSEL', IELSEL, NUMEL)
      CALL MDRSRV ('BLKSCR', IBLSC, NELBLK+1)
      call MDRSRV ('NBLSEL', INSEL, NELBLK+1)

      CALL MDSTAT (NERRS, NUSED)
      IF (NERRS .GT. 0) THEN
         CALL MEMERR
         STOP
      END IF
      CALL INILOG (NUMNP, .TRUE., IA(INDSEL))
      CALL INILOG (NUMEL, .TRUE., IA(IELSEL))
      NSELND = NUMNP
      NSELEL = NUMEL

   20 CONTINUE
      PRINT *, ' '
      CALL FREFLD (0, 0, PROMPT(:LENSTR(PROMPT)+1), MXNAM, IOS, NF, KV,
     *   CV, IVAL, RV)
      IF (IOS .NE. 0) CV(1) = 'EXIT'
      KV(MIN(MXNAM,NF)+1) = -999
      NAME = CV(1)
      CALL ABRSTR (NAME, CV(1), CMDTBL)
      IFLD = 2
      IF (NAME .EQ. ' ') NAME = CV(1)

C     ...SCAN FOR INPUT CARD MATCH

      IF (NAME.EQ.'        ') THEN
         GO TO 20

C-----------------------------------------------------------------------
      ELSE IF (NAME .EQ. 'ECHO') THEN
         CALL FFONOF (IFLD, KV, CV, ISON, *20)
         IF (ISON) THEN
            IOMIN = ITERM
            IF (IOMAX .EQ. IHARD) THEN
               CALL PRTERR ('CMDSPEC',
     *            'Results output to both the terminal and list file')
            ELSE
               CALL PRTERR ('CMDSPEC',
     *            'Results output to terminal only')
            END IF
         ELSE
            IOMIN = IHARD
            IOMAX = IHARD
            CALL PRTERR ('CMDSPEC',
     *         'Results output to list file only')
         END IF

C-----------------------------------------------------------------------
      ELSE IF (NAME .EQ. 'PRINT') THEN
         CALL FFONOF (IFLD, KV, CV, ISON, *20)
         IF (ISON) THEN
            IOMAX = IHARD
            IF (IOMIN .EQ. ITERM) THEN
               CALL PRTERR ('CMDSPEC',
     *            'Results output to both the terminal and list file')
            ELSE
               CALL PRTERR ('CMDSPEC',
     *            'Results output to list file only')
            END IF
         ELSE
            IOMIN = ITERM
            IOMAX = ITERM
            CALL PRTERR ('CMDSPEC', 'Results output to terminal only')
         END IF

C-----------------------------------------------------------------------
      ELSE IF (NAME .EQ. 'COMMENT') THEN
         IFLD = 2
         IF (MATSTR(CV(2), 'PAGE', 1)) THEN
            WRITE (IHARD, 40)
            IFLD = 3
         END IF
         ICOM = MIN(1, IVAL(IFLD))
         DO 30 IC = 1, ICOM
            CALL GETINP (0, 0, 'COMMENT> ', COMMENT, IOSTAT)
            WRITE (IHARD,50) COMMENT(:LENSTR(COMMENT))
   30    CONTINUE
   40    FORMAT ('1')
   50    FORMAT (1X,A)

C-----------------------------------------------------------------------
      ELSE IF (NAME .EQ. 'EXODUS') THEN
         CALL CKEXOD (EXOSAV, *20)
         CALL FFONOF (IFLD, KV, CV, EXODUS, *20)
         IF (EXODUS) THEN
            CALL PRTERR ('CMDSPEC',
     *         'Calculations performed for all selected time steps')
         ELSE
            CALL PRTERR ('CMDSPEC',
     *         'Calculations performed for undeformed geometry only')
         END IF

C-----------------------------------------------------------------------
      ELSE IF (NAME .EQ. 'HELP') THEN
         ISHELP = HELP ('NUMBERS', 'COMMANDS', CV(2) )
         IF (.NOT. ISHELP) CALL SHOCMD ('COMMANDS', CMDTBL)

C-----------------------------------------------------------------------
      ELSE IF (NAME .EQ. 'AXISYMME') THEN
         AXI = .TRUE.
         CALL PRTERR ('CMDSPEC', 'Axisymmetric Body Geometry')

C-----------------------------------------------------------------------
      ELSE IF (NAME .EQ. 'PLANAR' .OR. NAME .EQ. 'PLANE') THEN
         AXI = .FALSE.
         CALL PRTERR ('CMDSPEC', 'Planar Body Geometry')

C-----------------------------------------------------------------------
      ELSE IF (NAME .EQ.'PROPERTI' .OR. NAME .EQ. 'MASS') THEN

C ... SET NUMBER OF QUADRATURE POINTS

         CALL FFINTG (IFLD, KV, IVAL,
     *      'number of quadrature points', 1, NQUAD, *20)
         IF (NQUAD .NE. 1 .AND. NQUAD .NE. 2**NDIM) THEN
            IF (NDIM .EQ. 2) CALL PRTERR ('CMDERR',
     *         'quadrature order must be 1 or 4')
            IF (NDIM .EQ. 3) CALL PRTERR ('CMDERR',
     *         'quadrature order must be 1 or 8')
            GO TO 20
         END IF

         CALL FFREAL (IFLD, KV, RV,
     *      'common material density', 0.0, CDENS, *20)
         IF (CDENS .GT. 0.) THEN
            CALL INIREA (NELBLK, CDENS, DENSTY)
            CALL INISTR (NELBLK, ' ', LABEL)
         ELSE IF (CDENS .LT. 0.) THEN
            CALL PRTERR ('CMDERR',
     *         'density must be greater than zero')
            GO TO 20
         END IF

         CALL MASSPR (A, TIME, ITMSEL, DENSTY, MAT,
     *      DISP, NQUAD, LABEL)

C-----------------------------------------------------------------------
      ELSE IF (NAME .EQ. 'TIMESTEP') THEN

         CALL FFREAL (IFLD, KV, RV,
     *      'common material wavespeed', 0.0, CWAVE, *20)
         IF (CWAVE .GT. 0.) THEN
            CALL INIREA (NELBLK, CWAVE, WAVE)
            CALL INISTR (NELBLK, ' ', LABEL)
         ELSE IF (CWAVE .LT. 0.) THEN
            CALL PRTERR ('CMDERR',
     *         'wavespeed must be greater than zero')
            GO TO 20
         END IF

         IF (WAVE(1) .EQ. 0.) THEN
            CALL GETWAV (MAT, WAVE, NELBLK, LABEL)
         END IF

         CALL FFREAL (IFLD, KV, RV,
     *      'critical damping fraction', 0.06, EPSLON, *20)
         CALL ESTIME (CRD, WAVE, LINK, MAT, LABEL, NDIM, 2**NDIM,
     *      NELBLK, A(ISCR), A(ISCR2), EPSLON, NUMNP )

C-----------------------------------------------------------------------
      ELSE IF (NAME .EQ. 'DENSITY') THEN
         CALL GETDEN (MAT, DENSTY, NELBLK, LABEL)

C-----------------------------------------------------------------------
      ELSE IF (NAME .EQ. 'WAVESPEE') THEN
         CALL GETWAV (MAT, WAVE, NELBLK, LABEL)

C-----------------------------------------------------------------------
      ELSE IF (NAME .EQ. 'CAVITY') THEN
         NUMCAV = 0
   60    CONTINUE
         IF (FFNUMB (IFLD, KV)) THEN
           NUMCAV = MIN(NUMCAV + 1, MAXCAV)
           CALL FFINTG (IFLD, KV, IVAL,
     *       'cavity boundary flag', 0, ICAV(NUMCAV), *20)
           GO TO 60
         END IF

         CENT(1) = 0.0
         CENT(2) = 0.0
         CENT(3) = 0.0
         CENTER = .FALSE.
C ... For 3d, if no center specified use 0,0,0 unless user enters 'centroid'
         IF (NDIM .eq. 3) then
           center = .true.
         end if
         IF (FFMATC (IFLD, KV, CV, 'CENTROID', 5)) THEN
           CENTER = .FALSE.
         ELSE IF (FFMATC (IFLD, KV, CV, 'CENTER', 1)) THEN
           CALL FFREAL (IFLD, KV, RV, 'X coordinate of center',
     *       0.0, CENT(1), *20)
           CALL FFREAL (IFLD, KV, RV, 'Y coordinate of center',
     *       0.0, CENT(2), *20)
           IF (NDIM .EQ. 3) THEN
             CALL FFREAL (IFLD, KV, RV, 'Z coordinate of center',
     *         0.0, CENT(3), *20)
           END IF
           CENTER = .TRUE.
         ELSE
         END IF
         CALL CAVITY (A, CRD, A(IBC1), A(IBC2), A(IBC3), A(IBC4),
     *      A(IBC5),A(IBC6), A(IBC7), A(IBC8), DISP, NUMNP, NDIM,
     *      NUMESS, TIME, ITMSEL, TITLE, CENT, CENTER)

C-----------------------------------------------------------------------
      ELSE IF (NAME .EQ. 'GAP') THEN
         CALL FFINTG (IFLD, KV, IVAL,
     *      'gap master surface', 0, IMAS, *20)
         CALL FFINTG (IFLD, KV, IVAL,
     *      'gap slave surface', 0, ISLV, *20)
         CALL FFREAL (IFLD, KV, RV,
     *      'maximum gap matching distance', 0.0, DMAX, *20)
         CALL FFCHAR (IFLD, KV, CV, GMTHD, GMTHD)
         IF (MATSTR(GMTHD, 'DISTANCE', 1)) THEN
            GMTHD = 'DISTANCE'
         ELSE IF (MATSTR(GMTHD, 'NORMAL', 1)) THEN
            GMTHD = 'NORMAL'
         ELSE
            CALL PRTERR ('CMDERR', '"' // GMTHD(:LENSTR(GMTHD))
     &         // '" is an invalid matching option')
            GO TO 20
         END IF
         CALL GAPINI (A, CRD, A(IBC1), A(IBC2), A(IBC3), A(IBC4),
     *      A(IBC5),A(IBC6), A(IBC7), A(IBC8), DISP, NUMNP, NDIM,
     *      NUMESS, TIME, ITMSEL, TITLE, IMAS, ISLV,
     *      DMAX,GMTHD)

C-----------------------------------------------------------------------
      ELSE IF (NAME .EQ. 'OVERLAP') THEN
         CALL FFINTG (IFLD, KV, IVAL,
     *      'overlap master surface', 0, IMAS, *20)
         CALL FFINTG (IFLD, KV, IVAL,
     *      'overlap slave surface', 0, ISLV, *20)
         CALL OVRLAP(A, CRD, A(IBC1), A(IBC2), A(IBC3), A(IBC4),
     *      A(IBC5),A(IBC6), A(IBC7), A(IBC8), DISP, NUMNP, NDIM,
     *      NUMESS, TIME, ITMSEL, TITLE, IMAS, ISLV, NUMEL)

C-----------------------------------------------------------------------
      ELSE IF (NAME .EQ. 'SELECT') THEN
         IF (MATSTR(CV(2), 'ALL', 3)) THEN
            DO 70 IBLK = 1, NELBLK
               MAT(5,IBLK) = 1
   70       CONTINUE
            CALL INILOG(NUMNP, .TRUE., IA(INDSEL))
            CALL INILOG(NUMEL, .TRUE., IA(IELSEL))
            CALL PRTERR ('CMDSPEC',
     *         'All Blocks, Elements, and Nodes Selected')
         ELSE IF (MATSTR(CV(2), 'BLOCKS', 1)) THEN
            IF (MATSTR(CV(3), 'ALL', 1) .OR. NF .EQ. 2) THEN
               DO 80 IBLK = 1, NELBLK
                  MAT(5,IBLK) = 1
   80          CONTINUE
               CALL PRTERR ('CMDSPEC',
     *            'All Blocks, Elements, and Nodes Selected')
            ELSE
               DO 90 IBLK = 1, NELBLK
                  MAT(5,IBLK) = 0
   90          CONTINUE
               DO 100 IBLK = 3, NF
                  IF (IVAL(IBLK).GT.0 .AND. IVAL(IBLK).LE.NELBLK) THEN
                     MAT(5,IVAL(IBLK)) = 1
                  ELSE
                     COMMENT = 'Invalid block number: '
                     CALL FFADDI (IVAL(IBLK), COMMENT)
                     CALL PRTERR ('CMDERR', COMMENT(:LENSTR(COMMENT)))
                  END IF
  100          CONTINUE
            END IF
            CALL SELNOD (MAT, LINK, A(INDSEL), NUMNP, 2**NDIM, NELBLK,
     *         NSELND)
            CALL SELELM (MAT, A(IELSEL), NUMEL, NELBLK, NSELEL)
            CALL SHWBLK (NELBLK, MAT, NSELND, NSELEL)
         ELSE IF (MATSTR(CV(2), 'MATERIAL', 1)) THEN
            IF (MATSTR(CV(3), 'ALL', 1)) THEN
               DO 110 IBLK = 1, NELBLK
                  MAT(5,IBLK) = 1
  110          CONTINUE
            ELSE
               DO 120 IBLK = 1, NELBLK
                  MAT(5,IBLK) = 0
  120          CONTINUE
               DO 140 IFLD = 3, NF
                  IMAT = 0
                  DO 130 IBLK = 1, NELBLK
                     IF (IVAL(IFLD) .EQ. MAT(1,IBLK)) IMAT = IBLK
  130             CONTINUE
                  IF (IMAT .NE. 0) THEN
                     MAT(5,IMAT) = 1
                  ELSE
                     COMMENT = 'Invalid Material Number: '
                     CALL FFADDI (IVAL(IFLD), COMMENT)
                     CALL PRTERR ('CMDERR', COMMENT(:LENSTR(COMMENT)))
                  END IF
  140          CONTINUE
            END IF
            CALL SELNOD (MAT, LINK, A(INDSEL), NUMNP, 2**NDIM, NELBLK,
     *         NSELND)
            CALL SELELM (MAT, A(IELSEL), NUMEL, NELBLK, NSELEL)
            CALL SHWBLK (NELBLK, MAT, NSELND, NSELEL)
         ELSE IF (MATSTR(CV(2), 'SIDESETS', 1) .OR.
     *      MATSTR(CV(2), 'SSETS', 1)) THEN
            CALL SELSSN (A(INDSEL), NUMNP, NF-2, IVAL(3),
     *         A(IBC1), A(IBC3), A(IBC5), A(IBC7), NUMESS,
     *         NUMSEL)
         ELSE IF (MATSTR(CV(2), 'NODESETS', 1) .OR.
     *      MATSTR(CV(2), 'NSETS', 1)) THEN
            CALL SELSSN (A(INDSEL), NUMNP, NF-2, IVAL(3),
     *         A(INS1), A(INS2), A(INS3), A(INS4), NUMNPS,
     *         NUMSEL)
         ELSE IF (MATSTR(CV(2), 'BOX', 1)) THEN
            CALL SELBOX (A(IR), NUMNP, NDIM, RV(3), A(INDSEL), 'Nodes')
            CALL SELBOX (A(IECEN), NUMEL, NDIM, RV(3), A(IELSEL),
     *         'Elements')
         ELSE
            CALL PRTERR ('CMDERR', '"' // CV(2)(:LENSTR(CV(2)))
     &         // '" is an invalid or nonunique SELECT option')
         END IF

C-----------------------------------------------------------------------
      ELSE IF (NAME .EQ. 'LIST') THEN
         CALL ABRSTR (LISTYP, CV(2), LISTBL)
         IF (CV(2) .EQ. ' ') THEN
            CALL SHOCMD ('Valid LIST options', LISTBL)
         ELSE IF (LISTYP .EQ. ' ') THEN
            CALL PRTERR ('CMDERR', '"' // CV(2)(:LENSTR(CV(2)))
     &         // '" is an invalid or nonunique LIST option')
            CALL SHOCMD ('Valid LIST options', LISTBL)
         ELSE IF (LISTYP .EQ. 'SSETS' .OR. LISTYP .EQ. 'SIDESETS') THEN
            IF (NUMESS .EQ. 0) THEN
               CALL PRTERR ('CMDSPEC', 'No side sets to list')
            ELSE
               CALL SHOWFL ('S', NUMESS, A(IBC1), A(IBC2), A(IBC3))
            END IF
         ELSE IF (LISTYP .EQ. 'NSETS' .OR. LISTYP .EQ. 'NODESETS') THEN
            IF (NUMNPS .EQ. 0) THEN
               CALL PRTERR ('CMDSPEC', 'No node sets to list')
            ELSE
               CALL SHOWFL ('N', NUMNPS, A(INS1), A(INS2), A(1))
            END IF
         ELSE IF (LISTYP .EQ.'VARS' .OR. LISTYP .EQ.'VARIABLE') THEN
            CALL DBPINI ('*', NDB, TITLE, NDIM, NUMNP, NUMEL, NELBLK,
     *         NUMNPS, LNPSNL, LNPSDF, NUMESS, LESSEL, LESSNL, LESSDF,
     *         NVARGL, NVARNP, NVAREL, DBNAME)
         ELSE IF (LISTYP .EQ.'BLOCKS' .OR. LISTYP .EQ.'MATERIAL') THEN
            CALL SHWBLK (NELBLK, MAT, NSELND, NSELEL)
         ELSE IF (LISTYP .EQ. 'TIMES') THEN
            CALL CKEXOD (EXOSAV, *20)
            CALL DBPTIM ('*', NSTEP, TIME)
         ELSE IF (LISTYP .EQ. 'STEPS') THEN
            CALL CKEXOD (EXOSAV, *20)
            CALL DBPTIM ('NM', NSTEP, TIME)
         ELSE IF (LISTYP .EQ. 'STIMES' .OR.
     *      LISTYP .EQ. 'SELECTED' .AND. MATSTR(CV(3),'TIMES',1)) THEN
            CALL CKEXOD (EXOSAV, *20)
            CALL STIMES ('*', .FALSE., .TRUE., NSTEP, TIME, ITMSEL)
         ELSE IF (LISTYP .EQ. 'COMMANDS') THEN
            CALL SHOCMD ('COMMANDS', CMDTBL)
         ELSE IF (LISTYP .EQ. 'SORT') THEN
            CALL SHOCMD ('Valid SORT fields', SORTBL)
         ELSE IF (LISTYP .EQ. 'NAMES') THEN
            CALL CKEXOD (EXOSAV, *20)
            CALL DBPNAM ('*', NVARGL, NVARNP, NVAREL,
     *        NAMEGL, NAMENV, NAMEEL)
         ELSE IF (LISTYP .EQ. 'SELECTED') THEN
            IF (MATSTR(CV(4),'RANGE',1)) THEN
               OPT = 'R'
            ELSE
               OPT = ' '
            END IF
            IF (MATSTR(CV(3),'NODES',1)) THEN
               CALL LISSEL (OPT//'L', 'Nodes', IOMIN, IOMAX,
     *            IDUM, A(INDSEL), NUMNP)
            ELSE IF (MATSTR(CV(3),'ELEMENTS',1)) THEN
               CALL LISSEL (OPT//'L', 'Elements', IOMIN, IOMAX,
     *            IDUM, A(IELSEL), NUMEL)
            ELSE
               CALL PRTERR ('CMDERR', '"' // CV(3)(:LENSTR(CV(3)))
     &            // '" is an invalid LIST SELECTED option')
            END IF
         ELSE IF (LISTYP .EQ. 'QA' .OR. LISTYP .EQ. 'INFORMAT') THEN
            IF ((NQAREC .GT. 0) .OR. (NINFO .GT. 0)) THEN
               CALL DBPQA ('*', NQAREC, QAREC, NINFO, INFREC)
            END IF

         ELSE IF (LISTYP .EQ. 'VOLUME') THEN
           CALL PRVOL (ndim, CRD, link, numnp, numel, 8,
     &       a(ismp), IHARD)
           CALL PRTERR ('CMDSPEC',
     *       'Element Volumes were written to the list file')
         ELSE IF (LISTYP .EQ. 'NODALVOL') THEN
           CALL PRNVOL (ndim, CRD, link, numnp, numel, 8,
     &       a(ismp), IHARD)
           CALL PRTERR ('CMDSPEC',
     *       'Nodal Volumes were written to the list file')
         ELSE
            CALL PRTERR ('CMDERR', '"' // CV(2)(:LENSTR(CV(2)))
     &         // '" is an invalid or nonunique LIST option')
         END IF

C-----------------------------------------------------------------------
      ELSE IF (NAME .EQ. 'TMIN') THEN
         CALL CKEXOD (EXOSAV, *20)
         CALL FFREAL (IFLD, KV, RV,
     *      'minimum selected time', TMIN, STMIN, *20)
         CALL SELTIM (TIME, ITMSEL)

C-----------------------------------------------------------------------
      ELSE IF (NAME .EQ. 'TMAX') THEN
         CALL CKEXOD (EXOSAV, *20)
         CALL FFREAL (IFLD, KV, RV,
     *      'maximum selected time', TMAX, STMAX, *20)
         CALL SELTIM (TIME, ITMSEL)

C-----------------------------------------------------------------------
      ELSE IF (NAME .EQ.'TIMES') THEN
         CALL CKEXOD (EXOSAV, *20)
         CALL FFRRNG (IFLD, KV, CV, RV, 'time selection',
     *      TIME(NSTEP), TRANGE, *20)
         IF (TRANGE(3) .GT. 0.0) THEN
            STMIN = TRANGE(1)
            STMAX = TRANGE(2)
            STDEL = TRANGE(3)
         ELSE
            STMIN = TRANGE(2)
            STMAX = TRANGE(1)
            STDEL = TRANGE(3)
         END IF
         CALL SELTIM (TIME, ITMSEL)

C-----------------------------------------------------------------------
      ELSE IF (NAME .EQ. 'DELTIME') THEN
         CALL CKEXOD (EXOSAV, *20)
         CALL FFREAL (IFLD, KV, RV,
     *      'selected time interval', 0.0, STDEL, *20)
         CALL SELTIM (TIME, ITMSEL)

C-----------------------------------------------------------------------
      ELSE IF (NAME .EQ. 'ALLTIMES') THEN
         CALL CKEXOD (EXOSAV, *20)
         STMIN = TMIN
         STMAX = TMAX
         STDEL = 0.0
         CALL SELTIM (TIME, ITMSEL)
         EXODUS = .TRUE.

C-----------------------------------------------------------------------
      ELSE IF (NAME .EQ. 'NINTV') THEN
         CALL CKEXOD (EXOSAV, *20)
         CALL FFINTG (IFLD, KV, IVAL,
     *      'time step interval', 10, NINTV, *20)
         CALL SELINV (TIME, ITMSEL, NINTV)

C-----------------------------------------------------------------------
      ELSE IF (NAME .EQ. 'ZINTV') THEN
         CALL CKEXOD (EXOSAV, *20)
         CALL FFINTG (IFLD, KV, IVAL,
     *      'time step interval', 10, NINTV, *20)
         CALL SELINV (TIME, ITMSEL, -NINTV)

C-----------------------------------------------------------------------
      ELSE IF (NAME .EQ. 'EXIT' .OR. NAME .EQ. 'END'
     &        .or. name .eq. 'QUIT') THEN
         GO TO 190

C ... LOCATE NODES|ELEMENTS WITHIN toler OF LINE|PLANE|POINT

C --- FORMAT: locate nodes line  x1,y1,[z1] x2,y2,[z2], toler1, toler2 type
C             locate nodes plane x1,y1,[z1] i2,j2,[k2], toler1, toler2
C             locate nodes point x1,y1,[z1] toler1, toler2

C             locate elements line  x1,y1,[z1] x2,y2,[z2], toler1, toler2 type
C             locate elements plane x1,y1,[z1] i2,j2,[k2], toler1, toler2
C             locate elements point x1,y1,[z1] toler1, toler2

C             x1, y1, z1, x2, y2, z2 = Coordinate locations
C             i2, j2, k2 = Normal Vector to plane
C             If TOLER2 .EQ. 0, then TOLER1 = Maximum Distance for locate
C             If TOLER2 .NE. 0, then TOLER1 = Minimum Distance for locate,
C                                    TOLER2 = Maximum Distance for locate.
C             If TYPE .EQ. BOUNDED,   then only search within line
C             If TYPE .EQ. UNBOUNDED, then search along projection of line

C-----------------------------------------------------------------------
      ELSE IF (NAME .EQ. 'LOCATE') THEN

         TYPE = 'UNBOUNDE'
         DO 150 I=NF,4,-1
            IF (KV(I) .EQ. 0) THEN
               TYPE = CV(I)(:16)
               GO TO 160
            END IF
  150    CONTINUE
  160    CONTINUE
         IF (MATSTR(CV(2),'NODES',1)) THEN
            NUMLEN = NUMNP
            IPTR   = IR
            ISEL   = INDSEL
         ELSE IF (MATSTR(CV(2),'ELEMENTS',1)) THEN
            NUMLEN = NUMEL
            IPTR   = IECEN
            ISEL   = IELSEL
         ELSE
            CALL PRTERR ('CMDERR', '"' // CV(2)(:LENSTR(CV(2)))
     &         // '" is an invalid option')
            GO TO 170
         END IF

         CALL LOCTOL (CV(3), NDIM, RV, KV)

         IF (MATSTR(CV(3), 'LINE', 1)) THEN
            IF (NDIM .EQ. 3) THEN
               CALL LINE3 (A(IPTR), NUMLEN, A(ISCR), A(ISCR+NUMNP),
     *            NDIM, RV(4), RV(7), RV(10), CV(2), TYPE, SORTYP,
     *            A(ISMP), SORUP, IDUM, '*', A(ISEL))
            ELSE
               CALL LINE2 (A(IPTR), NUMLEN, A(ISCR), A(ISCR+NUMNP),
     *            NDIM, RV(4), RV(6), RV(8), CV(2), TYPE, SORTYP,
     *            A(ISMP), SORUP, IDUM, '*', A(ISEL))
            END IF
         ELSE IF (MATSTR(CV(3), 'PLANE', 2)) THEN
            IF (NDIM .EQ. 3) THEN
               CALL PLANE3 (A(IPTR), NUMLEN, A(ISCR), A(ISCR+NUMNP),
     *            NDIM, RV(4), RV(7), RV(10), CV(2), SORTYP,
     *            A(ISMP), SORUP, IDUM, '*', A(ISEL))
            ELSE
               CALL PRTERR ('CMDERR',
     *            'use locate node/element line for 2D')
            END IF
         ELSE IF (MATSTR(CV(3), 'POINT', 2)) THEN
            IF (NDIM .EQ. 3) THEN
               CALL POINT3 (A(IPTR), NUMLEN, A(ISCR),
     *            NDIM, RV(4), RV(7), CV(2), SORTYP,
     *            A(ISMP), A(ISCR2), SORUP, IDUM, '*', A(ISEL))
            ELSE
               CALL POINT2 (A(IPTR), NUMLEN, A(ISCR),
     *            NDIM, RV(4), RV(6), CV(2), SORTYP,
     *            A(ISMP), A(ISCR2), SORUP, IDUM, '*', A(ISEL))
            END IF
         ELSE
            CALL PRTERR ('CMDERR', '"' // CV(3)(:LENSTR(CV(3)))
     &         // '" is an invalid LOCATE option')
            GO TO 170
         END IF

  170    CONTINUE

C-----------------------------------------------------------------------
      ELSE IF (NAME .EQ. 'SORT') THEN
         IF (CV(2) .EQ. ' ')  CV(2) = 'NONE'
         IF (CV(3) .EQ. ' ')  CV(3) = 'UP'
         CALL ABRSTR (SORTYP, CV(2), SORTBL)
         IF (SORTYP .EQ. ' ') THEN
            CALL PRTERR ('CMDERR', '"' // CV(2)(:LENSTR(CV(2)))
     &         // '" is an invalid or nonunique SORT option')
            SORTYP = 'NONE'
         END IF
         IF (MATSTR(CV(3), 'UP', 1) .OR.
     *      MATSTR(CV(3), 'ASCENDIN', 1) .OR.
     *      MATSTR(CV(3), 'INCREASI', 1) ) THEN
            SORUP = .TRUE.
         ELSE IF (MATSTR(CV(3), 'DOWN', 1) .OR.
     *      MATSTR(CV(3), 'DESCENDI', 1) .OR.
     *      MATSTR(CV(3), 'DECREASI', 1) ) THEN
            SORUP = .FALSE.
         ELSE
            CALL PRTERR ('CMDERR', '"' // CV(3)(:LENSTR(CV(3)))
     &         // '" is an invalid SORT order')
            SORUP = .TRUE.
         END IF
         IF (SORUP) THEN
            CALL PRTERR ('CMDSPEC',
     *         'Sorting on field ' // SORTYP(:LENSTR(SORTYP)) //
     *         ' in ascending order.')
         ELSE
            CALL PRTERR ('CMDSPEC',
     *         'Sorting on field ' // SORTYP(:LENSTR(SORTYP)) //
     *         ' in descending order.')
         END IF

C-----------------------------------------------------------------------
      ELSE IF (NAME .EQ. 'LIMITS') THEN
         CALL MDRSRV ('XYZMIN', IXYZMN, NELBLK*NDIM)
         CALL MDRSRV ('XYZMAX', IXYZMX, NELBLK*NDIM)
         CALL MDSTAT (NERRS, NUSED)
         IF (NERRS .GT. 0) THEN
            CALL MEMERR
            STOP
         END IF
         IF (MATSTR(CV(2), 'ALLTIMES', 1) .AND. EXOSAV .AND. ISDIS) THEN
            ALLTIM = .TRUE.
         ELSE
            ALLTIM = .FALSE.
         END IF
         IF (ALLTIM) THEN
            CALL LIMITS (A(IXYZMN), A(IXYZMX), CRD, LINK, MAT, NDIM,
     *         NELBLK, 2**NDIM, ALLTIM, TIME, ITMSEL,
     *         DISP, NUMNP)
         ELSE
            CALL LIMITS (A(IXYZMN), A(IXYZMX), CRD, LINK, MAT, NDIM,
     *         NELBLK, 2**NDIM, ALLTIM, TIME, ITMSEL,
     *         CRD, NUMNP)
         END IF
         CALL MDDEL ('XYZMIN')
         CALL MDDEL ('XYZMAX')
         CALL MDSTAT (NERRS, NUSED)
         IF (NERRS .GT. 0) THEN
            CALL MEMERR
            STOP 'MEMORY'
         END IF
C-----------------------------------------------------------------------
      ELSE IF (NAME .EQ. 'SUM' .OR. NAME .EQ. 'AVERAGE') THEN
         LTMP2 = .FALSE.
         LTMP3 = .FALSE.
         CALL CKEXOD (EXOSAV, *20)
         IF (NAME .EQ. 'SUM') THEN
            LTMP = .FALSE.
         ELSE
            LTMP = .TRUE.
         END IF
         CALL FFCHAR (IFLD, KV, CV, ' ', NAME)
         IND = LOCSTR (NAME, NVARNP, NAMENV)
         IF (IND .EQ. 0) THEN
            CALL PRTERR ('ERROR', '"' // NAME(:LENSTR(NAME))
     &         // '" is not a valid nodal variable name')
         ELSE
            IPTR = INDSEL
            LTMP2 = .FALSE.
            LTMP3 = .FALSE.
  180       CONTINUE
            CALL FFCHAR (IFLD, KV, CV, 'RADIAL', CTMP)
            IF (MATSTR(CTMP, 'LINEAR', 1)) THEN
               LTMP2 = .FALSE.
            ELSE IF (MATSTR(CTMP, 'RADIAL', 1)) THEN
               LTMP2 = .TRUE.
            ELSE IF (MATSTR(CTMP, 'ABSOLUTE', 1)) THEN
               LTMP3 = .TRUE.
            ELSE
               CALL PRTERR ('ERROR', '"' // CTMP(:LENSTR(CTMP))
     &            // '" is not a valid option')
               GO TO 20
            END IF

            IF (FFEXST (IFLD, KV)) GO TO 180

            CALL SUMNOD (CRD, DISP, A(ISCR), NDIM, NUMNP, IND,
     *         A(IPTR), NAMENV(IND), TIME, ITMSEL,
     $         LTMP2, LTMP, LTMP3)
         END IF
C-----------------------------------------------------------------------
      ELSE IF (NAME .EQ. 'ESUM' .OR. NAME .EQ. 'EAVERAGE') THEN
         LTMP2 = .FALSE.
         LTMP3 = .FALSE.
         CALL CKEXOD (EXOSAV, *20)
         IF (NAME .EQ. 'ESUM') THEN
            LTMP = .FALSE.
         ELSE
            LTMP = .TRUE.
         END IF
         CALL FFCHAR (IFLD, KV, CV, ' ', NAME)
         IND = LOCSTR (NAME, NVAREL, NAMEEL)
         IF (IND .EQ. 0) THEN
            CALL PRTERR ('ERROR', '"' // NAME(:LENSTR(NAME))
     &         // '" is not a valid element variable name')
         ELSE
            IPTR = IELSEL
            LTMP4 = .FALSE.
            LTMP2 = .TRUE.
            LTMP3 = .FALSE.
 182        CONTINUE
            CALL FFCHAR (IFLD, KV, CV, 'DENSITY', CTMP)
C ... Determine whether element variable is per unit volume or is total
            IF (MATSTR(CTMP, 'DENSITY', 1)) THEN
               LTMP4 = .TRUE.
            ELSE IF (MATSTR(CTMP, 'TOTAL', 1)) THEN
               LTMP2 = .FALSE.
            ELSE IF (MATSTR(CTMP, 'ABSOLUTE', 1)) THEN
               LTMP3 = .TRUE.
            ELSE IF (MATSTR(CTMP, 'AVERAGE', 1)) THEN
               LTMP3 = .TRUE.
            ELSE
               CALL PRTERR ('ERROR', '"' // CTMP(:LENSTR(CTMP))
     &            // '" is not a valid option')
               GO TO 20
            END IF

            IF (FFEXST (IFLD, KV)) GO TO 182

            CALL SUMELM (CRD, DISP, A(ISCR), MAT, NDIM, NUMNP, IND,
     *         A(IPTR), NAMEEL(IND), TIME, ITMSEL,
     $         LTMP, LTMP2, LTMP3, LTMP4, numel, link, 2**ndim,
     &           nelblk, a(ismp), ISEVOK, A(IBLSC), A(INSEL))
         END IF
C ----------------------------------------
      ELSE IF (NAME .EQ. 'CONDITIO' ) THEN
         LTMP = .FALSE.
         CALL FFCHAR (IFLD, KV, CV, 'NODEBUG', CTMP)
         IF (MATSTR(CTMP, 'DEBUG', 1)) LTMP = .TRUE.
         CALL MDRSRV ('SUMRY',  ISUMR, 16*NELBLK)
         CALL MDRSRV ('ISUMR',  IISUM,  8*NELBLK)
         CALL MDRSRV ('ASPECT', IASPEC, NUMEL)
         CALL MDRSRV ('SKEW',   ISKEW,  NUMEL)
         CALL MDRSRV ('TAPER',  ITAPER, NUMEL)
         CALL MDRSRV ('AREA',   IAREA, NUMEL)
         IF (NDIM .EQ. 3) THEN
            CALL MDRSRV ('SKEWX', ISKX, NUMEL)
            CALL MDRSRV ('SKEWY', ISKY, NUMEL)
            CALL MDRSRV ('SKEWZ', ISKZ, NUMEL)
            CALL MDRSRV ('TAPRX', ITPX, NUMEL)
            CALL MDRSRV ('TAPRY', ITPY, NUMEL)
            CALL MDRSRV ('TAPRZ', ITPZ, NUMEL)
            CALL MDRSRV ('JACOB', IJAC, NUMEL)
         END IF
         CALL MDSTAT (NERRS, NUSED)
         IF (NERRS .GT. 0) THEN
            CALL MEMERR
            STOP 'MEMORY'
         END IF
         NNODES = 2**NDIM
         IF (NDIM .EQ. 2) THEN
            CALL CON2D(CRD, NDIM, NUMNP, LINK, NNODES, NUMEL, MAT,
     *         NELBLK, A(IELSEL), A(IASPEC), A(ISKEW), A(ITAPER),
     *         A(IAREA), A(ISUMR), A(IISUM), LTMP)
         ELSE
            CALL CON3D(CRD, NDIM, NUMNP, LINK, NNODES, NUMEL, MAT,
     *         NELBLK, A(IELSEL), A(IASPEC), A(ISKEW), A(ITAPER),
     *         A(IAREA), A(ISUMR), A(IISUM), A(ISKX), A(ISKY), A(ISKZ),
     *         A(ITPX), A(ITPY), A(ITPZ), A(IJAC), LTMP)
         END IF
         CALL MDDEL ('SUMRY')
         CALL MDDEL ('ISUMR')
         CALL MDDEL ('ASPECT')
         CALL MDDEL ('SKEW')
         CALL MDDEL ('TAPER')
         CALL MDDEL ('AREA')
         IF (NDIM .EQ. 3) THEN
            CALL MDDEL ('SKEWX')
            CALL MDDEL ('SKEWY')
            CALL MDDEL ('SKEWZ')
            CALL MDDEL ('TAPRX')
            CALL MDDEL ('TAPRY')
            CALL MDDEL ('TAPRZ')
            CALL MDDEL ('JACOB')
         END IF
         CALL MDSTAT (NERRS, NUSED)
         IF (NERRS .GT. 0) THEN
            CALL MEMERR
            STOP 'MEMORY'
         END IF
      ELSE
         CALL PRTERR ('CMDERR', '"' // NAME(:LENSTR(NAME))
     &      // '" is an invalid or nonunique command')
      END IF
      GO TO 20
  190 CONTINUE
      RETURN

      END
