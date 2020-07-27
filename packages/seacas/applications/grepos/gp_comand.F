C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C     -*- Mode: fortran -*-
C=======================================================================
      SUBROUTINE COMAND (NDBIN, EXECUT,
     &  IDELB, NUMELB, NUMLNK, NUMATR,
     &  IDNPS, NNNPS, NDNPS, IXNNPS, IXDNPS,
     $  LTNNPS, FACNPS,
     &  IDESS, NEESS, NNESS, IXEESS, IXNESS,
     $  LTEESS, LTSESS, FACESS,
     &  XN, YN, ZN, XEXPL, YEXPL, ZEXPL, MODBLK,
     *  ISATRB, ATRSCL, IEQUIV, MAPNOD, MAPEL,
     $  IELBST, INPSST, IESSST,
     *  NQAREC, QAREC, NINFO, INFREC, BLKTYP,
     *  ebname, nsname, ssname, atname,
     *  namegv, nvargl, namenv, nvarnp, nameev, nvarel,
     *  namemv, nvarns, namesv, nvarss, INOD2EL,
     *  SWPSS, SMOOTH, USRSUB, CENTRD,
     &  NSTEPS, TIMES, ITIMST, A, IA, *)
#if defined(__INTEL_COMPILER) && (__INTEL_COMPILER == 1300)
cDEC$ OPTIMIZE:2
#endif
C=======================================================================

C     --*** COMAND *** (GREPOS) Input and process commands
C     --   Written by Amy Gilkey - revised 05/09/88
C     --   Modified by Greg Sjaardema - 02/06/89
C     --
C     --COMAND inputs and executes an user command.
C     --
C     --Parameters:
C     --   IDNPS - IN - the IDs of existing node sets
C     --   IDESS - IN - the IDs of existing side sets
C     --   IBPARM - IN - the block parameters (defined by the block type)
C     --   IDELB - IN - the ids for each block
C     --   NUMELB - IN - the number of elements for each block
C     --   NUMLNK - IN - the number of nodes per element for each block
C     --   XN, YN, ZN - IN - the nodal coordinates
C     --   XEXPL, YEXPL, ZEXPL - IN - the offsets if EXPLODE option
C     --   * - return statement for QUIT
C     --
C     --Common Variables:
C     --   Sets NDIM, NUMNP, NUMEL, NELBLK,
C     --      NUMNPS, LNPSNL, NUMESS, LESSEL, LESSNL of /DBNUMS/
C     --   Sets DOTRAN, NNREPL, NEREPL, DIM3, NRTRAN, D3TRAN, ZGRAD,
C     --      CENTER, NUMCOL of /PARAMS/
C     --   Sets XOFFS, YOFFS, ZOFFS of /XYZOFF/
C     --   Sets ROT3D, ROTMAT, ROTCEN of /XYZROT/

      PARAMETER (MAXFLD = 64)

      include 'gp_params.blk'
      include 'gp_namlen.blk'
      include 'gp_dbnums.blk'
      include 'gp_dbtitl.blk'
      include 'gp_xyzoff.blk'
      include 'gp_xyzrot.blk'
      include 'gp_xyzmir.blk'
      include 'gp_xyzero.blk'
      include 'gp_xyzscl.blk'
      include 'gp_xyzwrp.blk'
      include 'gp_nsset.blk'
      include 'gp_splxyz.blk'
      include 'gp_smooth.blk'
      include 'gp_snap.blk'
      include 'gp_deform.blk'
      include 'gp_combine.blk'
      include 'gp_attrot.blk'

      INTEGER IDELB(*)
      INTEGER NUMELB(*)
      INTEGER NUMLNK(*)
      INTEGER NUMATR(*)
      INTEGER IELBST(*), INPSST(*), IESSST(*), ITIMST(*)
      INTEGER IEQUIV(*), MAPNOD(*), MAPEL(*)

      logical glonod, absolute

      INTEGER IDNPS(*), NNNPS(*), NDNPS(*), IXNNPS(*),
     $  IXDNPS(*), LTNNPS(*)
      REAL    FACNPS(*)

      INTEGER IDESS(*), NEESS(*), NNESS(*), IXEESS(*), IXNESS(*)
      INTEGER LTEESS(*), LTSESS(*)
      REAL    FACESS(*)
      REAL    TIMES(*)

      REAL XN(*), YN(*), ZN(*)
      REAL XEXPL(*), YEXPL(*), ZEXPL(*)
      REAL ATRSCL(2,*)
      REAL A(*)
      INTEGER IA(*)

      LOGICAL FFEXST, MATSTR, HELP

      CHARACTER*(128) STRING
      CHARACTER*(mxlnln) INFREC(*)
      CHARACTER*(mxstln) QAREC(4,*)
      CHARACTER*(mxstln) BLKTYP(*)
      character*(maxnam) ebname(*), nsname(*), ssname(*),
     *  namenv(*), nameev(*), namemv(*), namegv(*), namesv(*), ATNAME(*)

      INTEGER INOD2EL(*)
      CHARACTER*(256) CTEMP
      CHARACTER*8  STRA
      CHARACTER*1  ITYPE
      CHARACTER*(128)  WORD, VERB, LISTYP, WORD2, NEWNM, OLDNM
      INTEGER      INTYP(MAXFLD+1)
      CHARACTER*(128)  CFIELD(MAXFLD)
      INTEGER      IFIELD(MAXFLD)
      REAL         RFIELD(MAXFLD)
      LOGICAL      ISHELP, DELOK, SMOOTH
      LOGICAL      SWPSS, USRSUB, EXECUT, CENTRD
      LOGICAL      DOALLA, DOALLB, ISATRB, ATRWRN

      INTEGER      IRNG(3)

      REAL RNUM(9)
      CHARACTER*20 RSTR(9)

      LOGICAL      FIRST
      SAVE FIRST

      CHARACTER*(mxstln) MYNAME
      CHARACTER*(mxstln) CMDTBL(41), LISTBL(15)
      SAVE CMDTBL, LISTBL
C     --CMDTBL - the valid commands table

C     --Command table follows.  Remember to change the dimensioned size when
C     --changing the table.
      DATA CMDTBL /
     1  'OFFSET  ',  'REVOLVE ', 'REVCEN  ', 'MIRROR  ', 'NAME    ',
     2  'LIST    ',  'SHOW    ', 'HELP    ', 'ZERO    ', 'RENAME  ',
     3  'END     ',  'EXIT    ', 'QUIT    ', 'SCALE   ',
     4  'CHANGE  ',  'DELETE  ', 'EXPLODE ', 'SHIFT   ',
     5  'INCREMENT', 'LIMITS  ', 'ADJUST  ', 'RANDOMIZE',
     6  'SMOOTH  ',  'SWAP    ', 'USERSUBROUTINE', 'SNAP    ',
     7  'MOVE    ',  'WARP    ', 'EXECUTE ', 'COMBINE ', 'CENTROIDS',
     8  'UNDELETE',  'KEEP    ', 'EQUIVALENCE', 'DEFORM  ',
     9  'ELEMENTIZE','ROTATE  ', 'MERGE   ', 'TMIN    ', 'TMAX     ',
     $  '        ' /

      DATA LISTBL /
     *  'SSETS   ', 'SIDESETS', 'NSETS   ', 'NODESETS', 'VARS    ',
     *  'VARIABLE', 'BLOCKS  ', 'MATERIAL', 'COMMANDS', 'INFORMATION',
     *  'QA      ', 'TIMES   ', 'STEPS   ', 'NAMES   ', '        '/

      DATA MYNAME /'GREPOS'/
      DATA FIRST  /.TRUE./
C     --Initialize

      EXECUT = .FALSE.

      LNAM = LENSTR(MYNAME)
      MBLK = MAX(NELBLK, NUMNPS)
      CALL INIREA ( MBLK, 0.0, XEXPL)
      CALL INIREA ( MBLK, 0.0, YEXPL)
      IF (NDIM .EQ. 3) THEN
        CALL INIREA ( MBLK, 0.0, ZEXPL)
      END IF

      CALL INIINT(NVARNP, 0, INOD2EL)

C...Initialize attribute rotation
      REVATT = .FALSE.
      ROTALL = .FALSE.
      ROTTYP = .FALSE.
      ATTIND = -1
      ATTBLK = -1
      ROTBLK = "NONE"
      CALL INIREA (3*3, 0.0, ROTATT)
      DO I = 1, 3
        ROTATT(I,I) = 1.0
      END DO

C...Initialize attribute scaling
      ATRWRN = .FALSE.
      CALL SCLATR(NELBLK, IDELB, NUMATR, 0, .TRUE., 0, .TRUE., 1.0,
     *  ATRSCL, .FALSE., .TRUE.)

      if (FIRST) then
C... Initialization code that should only be executed once

C... Node equivalencing
        equiv = .FALSE.
        eqtoler = -1.0
        do i=1, numnp
          iequiv(i) = i
        end do

C...  Initialize Material, NodeSet, SideSet status
C     --      0 = same
C     --      - = delete
C     --      n = combine with entity n
        do 10 i=1, nelblk
          ielbst(i) = 0
 10     continue
        do 20 i=1, numnps
          inpsst(i) = 0
 20     continue
        do 30 i=1, numess
          iessst(i) = 0
 30     continue
        do 40 i=1, nsteps
          itimst(i) = 0
 40     continue
        first = .false.
      end if
C...Initialize transformation variables
      ADJTYP = 0
      XOFFS = 0.0
      YOFFS = 0.0
      ZOFFS = 0.0
      XMIRR = 1.0
      YMIRR = 1.0
      ZMIRR = 1.0
      XZERO = 0.0
      YZERO = 0.0
      ZZERO = 0.0
      XSCAL = 1.0
      YSCAL = 1.0
      ZSCAL = 1.0
      XRAND = 0.0
      YRAND = 0.0
      ZRAND = 0.0

      IWARP = 0
      MODBLK = 0
      IDEFST = 0

      SPLOFF = .FALSE.
      SWPSS  = .FALSE.
      SMOOTH = .FALSE.
      USRSUB = .FALSE.
      CENTRD = .FALSE.

      ROT3D  = .FALSE.
      CALL INIREA (3*3, 0.0, ROTMAT)
      DO 100 I = 1, 3
        ROTMAT(I,I) = 1.0
 100  CONTINUE

      CALL MINMAX (NUMNP, XN, XMIN, XMAX)
      ROTCEN(1) = XMIN
      CALL MINMAX (NUMNP, YN, YMIN, YMAX)
      ROTCEN(2) = YMIN
      IF (NDIM .EQ. 3) THEN
        CALL MINMAX (NUMNP, ZN, ZMIN, ZMAX)
        ROTCEN(3) = ZMIN
      ELSE
        ROTCEN(3) = 0.0
      END IF

C ... Initialize SNAP variables
      numsnp = 0

 110  CONTINUE

C     --Read command line

      WRITE (*, *)
      CALL FREFLD (0, 0, 'GREPOS> ', MAXFLD,
     &  IOSTAT, NUMFLD, INTYP, CFIELD, IFIELD, RFIELD)
      IF (IOSTAT .LT. 0) GOTO 430
      IF (NUMFLD .EQ. 0) GOTO 110
      INTYP(MIN(MAXFLD,NUMFLD)+1) = -999

      IFLD = 1
      CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
      CALL ABRSTR (VERB, WORD, CMDTBL)
      IF (VERB .EQ. ' ') VERB = WORD

C     --Perform command

      IF (VERB .EQ. '?') THEN
        CALL SHOCMD ('COMMANDS', CMDTBL)
        VERB = ' '
C=======================================================================
      ELSE IF (VERB .EQ. 'OFFSET' .OR. VERB .EQ. 'SHIFT') THEN

C     ... Originally offset just asked for NDIM values.  It was changed
C     to go by axis type (OFFSET Y 1.0). We do maintain compatibility
C     and check for both types of input.

        RMULT = 0.0
        IF ((INTYP(IFLD) .EQ. 0) .OR. ((CFIELD(IFLD) .GE. 'A')
     &    .AND. (CFIELD(IFLD) .LE. 'Z'))) THEN
 120      CONTINUE
          IF (FFEXST (IFLD, INTYP)) THEN
            CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
            IF (MATSTR (WORD, 'SPLINE', 1)) THEN
              SPLOFF = .TRUE.
              CALL GETSPL(A)
            ELSE IF (MATSTR (WORD, 'RESET', 1)) THEN
              XOFFS = 0.0
              YOFFS = 0.0
              ZOFFS = 0.0
              RMULT = 0.0
            ELSE IF (MATSTR (WORD, 'ADD', 2)) THEN
C     ... Set for cumulative offsets/shifts
              RMULT = 1.0
            ELSE IF (MATSTR (WORD, 'ALL', 2)) THEN
              CALL FFREAL (IFLD, INTYP, RFIELD,
     &          'model offset', 0.0, TOFFS, *130)
              XOFFS = RMULT * XOFFS + TOFFS
              YOFFS = RMULT * YOFFS + TOFFS
              ZOFFS = RMULT * ZOFFS + TOFFS
            ELSE IF (WORD .EQ. 'X') THEN
              CALL FFREAL (IFLD, INTYP, RFIELD,
     &          'X coordinate offset', 0.0, TOFFS, *130)
              XOFFS = RMULT * XOFFS + TOFFS
            ELSE IF (WORD .EQ. 'Y') THEN
              CALL FFREAL (IFLD, INTYP, RFIELD,
     &          'Y coordinate offset', 0.0, TOFFS, *130)
              YOFFS = RMULT * YOFFS + TOFFS
            ELSE IF (WORD .EQ. 'Z') THEN
              IF (NDIM .EQ. 3) THEN
                CALL FFREAL (IFLD, INTYP, RFIELD,
     &            'Z coordinate offset', 0.0, TOFFS, *130)
                ZOFFS = RMULT * ZOFFS + TOFFS
              ELSE
                CALL PRTERR ('CMDERR',
     *            'Z allowed for 3D database only')
              END IF
            ELSE
              IF (NDIM .EQ. 3) CALL PRTERR ('CMDERR',
     &          'Expected "X", "Y", "Z", "ALL", "ADD", or "RESET"')
              IF (NDIM .EQ. 2) CALL PRTERR ('CMDERR',
     &          'Expected "X", "Y", "ALL", "ADD", or "RESET"')
              GOTO 130
            END IF
            GOTO 120
          END IF
        ELSE
          CALL FFREAL (IFLD, INTYP, RFIELD,
     &      'X coordinate offset', XOFFS, XOFFS, *130)
          CALL FFREAL (IFLD, INTYP, RFIELD,
     &      'Y coordinate offset', YOFFS, YOFFS, *130)
          IF (NDIM .EQ. 3) CALL FFREAL (IFLD, INTYP, RFIELD,
     &      'Z coordinate offset', ZOFFS, ZOFFS, *130)
        END IF
 130    CONTINUE

C=======================================================================
      ELSE IF (VERB .EQ. 'ADJUST') THEN
C...  ADJUST {CENTER|MINIMUM|MAXIMUM} {X|Y|Z} adjustment
C...  ADJUST {MINIMUM|MAXIMUM} {X|Y|Z} adjustment

C...  Calculate offset such that the value in the generated mesh
C     (ignoring any other offsets, scales, or rotations)
C     will equal the input value.

C...EXAMPLE:adjust maximum x 10
C     will set the maximum x in the generated mesh equal to 10.

 140    CONTINUE
        IF (FFEXST (IFLD, INTYP)) THEN
          CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
          IF (MATSTR(WORD, 'MAXIMUM', 3)) THEN
            ADJTYP = 1
          ELSE IF (MATSTR(WORD, 'MINIMUM', 3)) THEN
            ADJTYP = 2
          ELSE IF (MATSTR(WORD, 'CENTER', 3)) THEN
            ADJTYP = 3
          ELSE
            CALL PRTERR ('CMDERR',
     &        'Expected "MINIMUM", "MAXIMUM", or "CENTER"')
            GO TO 150
          END IF
        END IF
        IF (FFEXST (IFLD, INTYP)) THEN
          CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
          IF (WORD .EQ. 'X') THEN
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &        'X adjustment', 0.0, TOFFS, *150)
            IF (ADJTYP .EQ. 1) THEN
              XOFFS = TOFFS - XMAX
            ELSE IF (ADJTYP .EQ. 2) THEN
              XOFFS = TOFFS - XMIN
            ELSE IF (ADJTYP .EQ. 3) THEN
              XOFFS = TOFFS - (XMAX+XMIN)/2.0
            END IF
          ELSE IF (WORD .EQ. 'Y') THEN
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &        'Y adjustment', 0.0, TOFFS, *150)
            IF (ADJTYP .EQ. 1) THEN
              YOFFS = TOFFS - YMAX
            ELSE IF (ADJTYP .EQ. 2) THEN
              YOFFS = TOFFS - YMIN
            ELSE IF (ADJTYP .EQ. 3) THEN
              YOFFS = TOFFS - (YMAX+YMIN)/2.0
            END IF
          ELSE IF (WORD .EQ. 'Z') THEN
            IF (NDIM .EQ. 3) THEN
              CALL FFREAL (IFLD, INTYP, RFIELD,
     &          'Z adjustment', 0.0, TOFFS, *150)
              IF (ADJTYP .EQ. 1) THEN
                ZOFFS = TOFFS - ZMAX
              ELSE IF (ADJTYP .EQ. 2) THEN
                ZOFFS = TOFFS - ZMIN
              ELSE IF (ADJTYP .EQ. 3) THEN
                ZOFFS = TOFFS - (ZMAX+ZMIN)/2.0
              END IF
            ELSE
              CALL PRTERR ('CMDERR',
     *          'Z allowed for 3D database only')
            END IF
          ELSE
            IF (NDIM .EQ. 3) CALL PRTERR ('CMDERR',
     &        'Expected "X", "Y", or "Z"')
            IF (NDIM .EQ. 2) CALL PRTERR ('CMDERR',
     &        'Expected "X" or "Y"')
            GOTO 150
          END IF
          GOTO 140
        END IF
 150    CONTINUE
        VERB = 'OFFSET'

C=======================================================================
      ELSE IF (VERB .EQ. 'LIMITS') THEN
        call limits('Input Mesh Limits:', ndim,
     &    xmin, xmax, ymin, ymax, zmin, zmax)
        VERB = ' '
C=======================================================================
      ELSE IF (VERB .EQ. 'EXPLODE') THEN
        DO 170 IBLK = 1, NELBLK
          CALL INTSTR (1, 0, IDELB(IBLK), WORD, LW)
          WRITE (CTEMP, 10010) 'Offset for Material ', WORD(:LW),'> '
10010     FORMAT (1X, 3A)
          CALL FREFLD (0, 0, CTEMP(:LENSTR(CTEMP)+1), MAXFLD,
     &      IOSTAT, NUMFLD, INTYP, CFIELD, IFIELD, RFIELD)
          IF (IOSTAT .LT. 0) GOTO 430
          INTYP(MIN(MAXFLD,NUMFLD)+1) = -999
          IFLD = 1
          CALL FFREAL (IFLD, INTYP, RFIELD,
     &      CTEMP, 0.0, XEXPL(IBLK), *160)
          CALL FFREAL (IFLD, INTYP, RFIELD,
     &      CTEMP, 0.0, YEXPL(IBLK), *160)
          IF (NDIM .EQ. 3) THEN
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &        CTEMP, 0.0, ZEXPL(IBLK), *160)
          END IF
 160      CONTINUE
 170    CONTINUE
        MODBLK = 1
        VERB = ' '

C=======================================================================
      ELSE IF (VERB .EQ. 'ZERO') THEN
 180    CONTINUE
        IF (FFEXST (IFLD, INTYP)) THEN
          CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
          IF (MATSTR (WORD, 'RESET', 1)) THEN
            XZERO = 0.0
            YZERO = 0.0
            ZZERO = 0.0
          ELSE IF (WORD .EQ. 'X') THEN
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &        'minimum X coordinate', 0.0, XZERO, *190)
          ELSE IF (WORD .EQ. 'Y') THEN
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &        'minimum Y coordinate', 0.0, YZERO, *190)
          ELSE IF (WORD .EQ. 'Z') THEN
            IF (NDIM .EQ. 3) THEN
              CALL FFREAL (IFLD, INTYP, RFIELD,
     &          'minimum Z coordinate', 0.0, ZZERO, *190)
            ELSE
              CALL PRTERR ('CMDERR',
     *          'Z allowed for 3D database only')
            END IF
          ELSE
            CALL PRTERR ('CMDERR',
     &        'Expected "X", "Y", "Z" or "RESET"')
            GOTO 190
          END IF
          GOTO 180
        END IF

 190    CONTINUE

C=======================================================================
      ELSE IF (VERB .EQ. 'SCALE') THEN
 200    CONTINUE
        IF (FFEXST (IFLD, INTYP)) THEN
          CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)

          IF (.NOT. MATSTR(WORD, 'ATTRIBUTE', 2) .AND.
     *      ISATRB .AND. .NOT. ATRWRN) THEN
            CALL PRTERR ('WARNING',
     *        'Some element blocks have attributes which ' //
     *        'may also need scaling.')
            ATRWRN = .TRUE.
          END IF

          IF (MATSTR (WORD, 'RESET', 1)) THEN
            XSCAL = 1.0
            YSCAL = 1.0
            ZSCAL = 1.0

          ELSE IF (MATSTR (WORD, 'ALL', 2)) THEN
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &        'model scale factor', 1.0, TSCAL, *230)
            XSCAL = XSCAL * ABS(TSCAL)
            YSCAL = YSCAL * ABS(TSCAL)
            ZSCAL = ZSCAL * ABS(TSCAL)

          ELSE IF (WORD .EQ. 'X') THEN
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &        'X scale factor', 1.0, TSCAL, *230)
            XSCAL = XSCAL * ABS(TSCAL)

          ELSE IF (WORD .EQ. 'Y') THEN
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &        'Y scale factor', 1.0, TSCAL, *230)
            YSCAL = YSCAL * ABS(TSCAL)

          ELSE IF (WORD .EQ. 'Z') THEN
            IF (NDIM .EQ. 3) THEN
              CALL FFREAL (IFLD, INTYP, RFIELD,
     &          'Z scale factor', 1.0, TSCAL, *230)
              ZSCAL = ZSCAL * ABS(TSCAL)
            ELSE
              CALL PRTERR ('CMDERR',
     *          'Z allowed for 3D database only')
            END IF

          ELSE IF (MATSTR(WORD, 'ATTRIBUTE', 2)) THEN
C ... SCALE ATTRIBUTE {num|ALL} BLOCK {id|ALL} scale
C ... SCALE ATTRIBUTE RESET
C ... See if next field is integer (id) or word ('all')
            VERB = ' '
            STRING = 'SCALE ATTRIBUTE {id|ALL} BLOCK {id|ALL} {scale}'
            IF (.NOT. ISATRB) THEN
              CALL PRTERR ('CMDERR',
     *          'There are no attributes on this database')
              GO TO 230
            END IF

            ATRWRN = .TRUE.
            DOALLA = .FALSE.
            DOALLB = .FALSE.

            IF ((INTYP(IFLD) .EQ. 0) .OR. ((CFIELD(IFLD) .GE. 'A')
     &        .AND. (CFIELD(IFLD) .LE. 'Z'))) THEN
C ... Field is character
              CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
              IF (MATSTR (WORD, 'ALL', 2)) THEN
                DOALLA = .TRUE.
              ELSE IF (MATSTR (WORD, 'RESET', 2)) THEN
                CALL SCLATR(NELBLK, IDELB, NUMATR, 0, .TRUE., 0,
     *            .TRUE., 1.0, ATRSCL, .TRUE., .TRUE.)
                GO TO 230
              ELSE
                CALL PRTERR('CMDERR',
     *            'Expected keyword "ALL"')
                CALL PRTERR('CMDREQ', STRING(:LENSTR(STRING)))
                GO TO 230
              END IF
            ELSE
C ... Field is Attribute number
              CALL FFINTG (IFLD, INTYP, IFIELD,
     $          'attribute ID', 0, IDATR, *420)
            END IF

C ... Next word should be 'BLOCK'
            CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
            IF (.NOT.(MATSTR (WORD, 'BLOCK', 1))) THEN
              CALL PRTERR('CMDERR',
     *          'Expected keyword "BLOCK"')
              CALL PRTERR('CMDREQ', STRING(:LENSTR(STRING)))
              GO TO 230
            END IF
C ... Now find block id or 'ALL'
            IF ((INTYP(IFLD) .EQ. 0) .OR. ((CFIELD(IFLD) .GE. 'A')
     &        .AND. (CFIELD(IFLD) .LE. 'Z'))) THEN
C ... Field is character
              CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
              IF (MATSTR (WORD, 'ALL', 2)) THEN
                DOALLB = .TRUE.
              ELSE
                CALL PRTERR('CMDERR',
     *            'Expected keyword "ALL"')
                CALL PRTERR('CMDREQ', STRING(:LENSTR(STRING)))
                GO TO 230
              END IF
            ELSE
C ... Field is Attribute number
              CALL FFINTG (IFLD, INTYP, IFIELD,
     $          'block ID', 0, IDBLK, *420)
            END IF

            CALL FFREAL (IFLD, INTYP, RFIELD,
     *        'attribute scale factor', 1.0, ASCALE, *230)

C ... Setup the scale and echo to the user
            CALL SCLATR(NELBLK, IDELB, NUMATR, IDATR, DOALLA,
     *        IDBLK, DOALLB, ASCALE, ATRSCL, .TRUE., .TRUE.)

          ELSE IF (MATSTR(WORD, 'BLOCK', 1)) THEN
            DO 220 IBLK = 1, NELBLK
              CALL INTSTR (1, 0, IDELB(IBLK), WORD, LW)
              WRITE (CTEMP, 10020)
     $          'Scale factor for Material ', WORD(:LW),'> '
10020         FORMAT (1X, 3A)
              CALL FREFLD (0, 0, CTEMP(:LENSTR(CTEMP)+1), MAXFLD,
     &          IOSTAT, NUMFLD, INTYP, CFIELD, IFIELD, RFIELD)
              IF (IOSTAT .LT. 0) GOTO 430
              INTYP(MIN(MAXFLD,NUMFLD)+1) = -999
              IFLD = 1
              CALL FFREAL (IFLD, INTYP, RFIELD,
     &          CTEMP, 0.0, XEXPL(IBLK), *210)
              CALL FFREAL (IFLD, INTYP, RFIELD,
     &          CTEMP, 0.0, YEXPL(IBLK), *210)
              IF (NDIM .EQ. 3) THEN
                CALL FFREAL (IFLD, INTYP, RFIELD,
     &            CTEMP, 0.0, ZEXPL(IBLK), *210)
              END IF
 210          CONTINUE
 220        CONTINUE
            MODBLK = 2
            VERB = ' '
          ELSE
            IF (NDIM .EQ. 3) CALL PRTERR ('CMDERR',
     &        'Expected "X", "Y", "Z", "ALL", "ATTRIBUTE", or "RESET"')
            IF (NDIM .EQ. 2) CALL PRTERR ('CMDERR',
     &        'Expected "X", "Y", "ALL", "ATTRIBUTE", or "RESET"')
            GOTO 230
          END IF
          GOTO 200
        END IF

 230    CONTINUE

C=======================================================================
      ELSE IF (VERB .EQ. 'DEFORM') THEN
        IF (FFEXST (IFLD, INTYP)) THEN
          CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
          IF (MATSTR (WORD, 'RESET', 1)) THEN
            IDEFST = 0
          ELSE IF (WORD .EQ. 'STEP') THEN
            CALL FFINTG (IFLD, INTYP, IFIELD,
     &        'deform step', 0, IDEFST, *235)
          ELSE
            CALL PRTERR ('CMDERR',
     &        'Expected "STEP", or "RESET"')
            GOTO 235
          END IF
        END IF

 235    CONTINUE

C=======================================================================
      ELSE IF (VERB .EQ. 'RANDOMIZE') THEN
 240    CONTINUE
        IF (FFEXST (IFLD, INTYP)) THEN
          CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
          IF (MATSTR (WORD, 'RESET', 1)) THEN
            XRAND = 0.0
            YRAND = 0.0
            ZRAND = 0.0
          ELSE IF (MATSTR (WORD, 'ALL', 1)) THEN
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &        'model random factor', 0.0, TRAND, *270)
            XRAND = ABS(TRAND)
            YRAND = ABS(TRAND)
            ZRAND = ABS(TRAND)
          ELSE IF (WORD .EQ. 'X') THEN
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &        'X random factor', 0.0, TRAND, *270)
            XRAND = ABS(TRAND)
          ELSE IF (WORD .EQ. 'Y') THEN
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &        'Y random factor', 0.0, TRAND, *270)
            YRAND = ABS(TRAND)
          ELSE IF (WORD .EQ. 'Z') THEN
            IF (NDIM .EQ. 3) THEN
              CALL FFREAL (IFLD, INTYP, RFIELD,
     &          'Z random factor', 0.0, TRAND, *270)
              ZRAND = ABS(TRAND)
            ELSE
              CALL PRTERR ('CMDERR',
     *          'Z allowed for 3D database only')
            END IF
          ELSE IF (MATSTR(WORD, 'BLOCK', 1)) THEN
            DO 260 IBLK = 1, NELBLK
              CALL INTSTR (1, 0, IDELB(IBLK), WORD, LW)
              WRITE (CTEMP, 10030)
     $          'Random factor for Material ', WORD(:LW),'> '
10030         FORMAT (1X, 3A)
              CALL FREFLD (0, 0, CTEMP(:LENSTR(CTEMP)+1), MAXFLD,
     &          IOSTAT, NUMFLD, INTYP, CFIELD, IFIELD, RFIELD)
              IF (IOSTAT .LT. 0) GOTO 430
              INTYP(MIN(MAXFLD,NUMFLD)+1) = -999
              IFLD = 1
              CALL FFREAL (IFLD, INTYP, RFIELD,
     &          CTEMP, 0.0, XEXPL(IBLK), *250)
              CALL FFREAL (IFLD, INTYP, RFIELD,
     &          CTEMP, 0.0, YEXPL(IBLK), *250)
              IF (NDIM .EQ. 3) THEN
                CALL FFREAL (IFLD, INTYP, RFIELD,
     &            CTEMP, 0.0, ZEXPL(IBLK), *250)
              END IF
 250          CONTINUE
 260        CONTINUE
            MODBLK = 3
            VERB = ' '
          ELSE IF (MATSTR(WORD,'NODESETS',3) .OR.
     $        MATSTR(WORD,'NSETS',   3) ) THEN
            DO 265 IBLK = 1, NUMNPS
              CALL INTSTR (1, 0, IDNPS(IBLK), WORD, LW)
              WRITE (CTEMP, 10030)
     $          'Random factors (x y z) for Nodeset ', WORD(:LW),'> '
              CALL FREFLD (0, 0, CTEMP(:LENSTR(CTEMP)+1), MAXFLD,
     &          IOSTAT, NUMFLD, INTYP, CFIELD, IFIELD, RFIELD)
              IF (IOSTAT .LT. 0) GOTO 430
              INTYP(MIN(MAXFLD,NUMFLD)+1) = -999
              IFLD = 1
              CALL FFREAL (IFLD, INTYP, RFIELD,
     &          CTEMP, 0.0, XEXPL(IBLK), *255)
              CALL FFREAL (IFLD, INTYP, RFIELD,
     &          CTEMP, 0.0, YEXPL(IBLK), *255)
              IF (NDIM .EQ. 3) THEN
                CALL FFREAL (IFLD, INTYP, RFIELD,
     &            CTEMP, 0.0, ZEXPL(IBLK), *255)
              END IF
 255          CONTINUE
 265        CONTINUE
            MODBLK = 4
            VERB = ' '
          ELSE
            IF (NDIM .EQ. 3) CALL PRTERR ('CMDERR',
     &        'Expected "X", "Y", "Z", "ALL", or "RESET"')
            IF (NDIM .EQ. 2) CALL PRTERR ('CMDERR',
     &        'Expected "X", "Y", "ALL", or "RESET"')
            GOTO 270
          END IF
          GOTO 240
        END IF

 270    CONTINUE

C=======================================================================
      ELSE IF (VERB .EQ. 'MIRROR') THEN

 280    CONTINUE
        IF (FFEXST (IFLD, INTYP)) THEN
          CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
          IF (MATSTR (WORD, 'RESET', 1)) THEN
            XMIRR = 1.
            YMIRR = 1.
            ZMIRR = 1.
          ELSE IF (WORD .EQ. 'X') THEN
            XMIRR = -1.
          ELSE IF (WORD .EQ. 'Y') THEN
            YMIRR = -1.
          ELSE IF (WORD .EQ. 'Z') THEN
            IF (NDIM .EQ. 3) THEN
              ZMIRR = -1.
            ELSE
              CALL PRTERR ('CMDERR',
     *          'Z not allowed for 2D database')
            END IF
          ELSE
            IF (NDIM .EQ. 3) CALL PRTERR ('CMDERR',
     &        'Expected "X", "Y", "Z" or "RESET"')
            IF (NDIM .EQ. 2) CALL PRTERR ('CMDERR',
     &        'Expected "X", "Y", or "RESET"')
            GOTO 290
          END IF
          GOTO 280
        END IF

 290    CONTINUE

C=======================================================================
      ELSE IF (VERB .EQ. 'REVOLVE' .OR. VERB .EQ. 'ROTATE') THEN
        DEGANG = ATAN2(0.0, -1.0) / 180.0

 300    CONTINUE
        IF (FFEXST (IFLD, INTYP)) THEN
          CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
          IF (MATSTR (WORD, 'RESET', 1)) THEN
            ROT3D = .FALSE.
            CALL INIREA (3*3, 0.0, ROTMAT)
            DO 310 I = 1, 3
              ROTMAT(I,I) = 1.0
 310        CONTINUE
          ELSE IF (NDIM .EQ. 3 .AND.
     *        ((WORD .EQ. 'X') .OR. (WORD .EQ. 'Y')
     &        .OR. (WORD .EQ. 'Z')) ) THEN
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &        'angle of rotation', 0.0, DEG, *320)
            ROT3D = .TRUE.
            CALL ROTXYZ (WORD, DEG * DEGANG, ROTMAT)
          ELSE IF (NDIM .EQ. 2 .AND. (WORD .EQ. 'Z') ) THEN
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &        'angle of rotation', 0.0, DEG, *320)
            ROT3D = .TRUE.
            CALL ROTXYZ (WORD, DEG * DEGANG, ROTMAT)

          ELSE IF (MATSTR(WORD, 'ATTRIBUTE', 2)) THEN
C ... REVOLVE ATTRIBUTE {which} BLOCK {id|ALL}  {x|y|z} {angle}
C ... REVOLVE ATTRIBUTE {which} TYPE {blk_type}
C ... REVOLVE ATTRIBUTE RESET
C ... Next field is integer (which)
            VERB = ' '
            STRING = 'REVOLVE/ROTATE ATTRIBUTE {index} BLOCK {id|ALL}'
            IF (.NOT. ISATRB) THEN
              CALL PRTERR ('CMDERR',
     *          'There are no attributes on this database')
              GO TO 320
            END IF

            ROTALL = .FALSE.
            IF ((INTYP(IFLD) .EQ. 0) .OR. ((CFIELD(IFLD) .GE. 'A')
     &        .AND. (CFIELD(IFLD) .LE. 'Z'))) THEN
C ... Field is character
              CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
              IF (MATSTR (WORD, 'HELP', 2)) THEN
                CALL PRTERR('CMDREQ',
     *'REVOLVE ATTRIBUTE {index} BLOCK {id|ALL} {X|Y|Z} {angle}')
                CALL PRTERR('CMDREQ',
     *'REVOLVE ATTRIBUTE {index} TYPE {type} {X|Y|Z} {angle}')
                CALL PRTERR('CMDREQ',
     *'REVOLVE ATTRIBUTE RESET')
              ELSE IF (MATSTR (WORD, 'RESET', 2)) THEN
                REVATT = .FALSE.
                CALL INIREA (3*3, 0.0, ROTATT)
                DO I = 1, 3
                  ROTATT(I,I) = 1.0
                END DO
              ELSE
                CALL PRTERR('CMDERR',
     *            'Expected keyword "RESET"')
                CALL PRTERR('CMDREQ', STRING(:LENSTR(STRING)))
                GO TO 320
              END IF
            ELSE
C ... Field is Attribute index -- first of ndim attributes involed in revolve. 1-based
              CALL FFINTG (IFLD, INTYP, IFIELD,
     $          'attribute index', 0, ATTIND, *420)
            END IF

C ... Next word should be 'BLOCK'
            CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
            IF ((MATSTR (WORD, 'BLOCK', 1))) THEN
C ... Now find block id or 'ALL'
              ROTTYP = .FALSE.
              IF ((INTYP(IFLD) .EQ. 0) .OR. ((CFIELD(IFLD) .GE. 'A')
     &          .AND. (CFIELD(IFLD) .LE. 'Z'))) THEN
C ... Field is character
                CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
                IF (MATSTR (WORD, 'ALL', 2)) THEN
                  ROTALL = .TRUE.
                ELSE
                  CALL PRTERR('CMDERR',
     *              'Expected keyword "ALL"')
                  CALL PRTERR('CMDREQ', STRING(:LENSTR(STRING)))
                  GO TO 320
                END IF
              ELSE
C ... Field is block id
                ROTALL = .FALSE.
                CALL FFINTG (IFLD, INTYP, IFIELD,
     $            'block ID', 0, ATTBLK, *420)
              END IF

            ELSE IF ((MATSTR (WORD, 'TYPE', 1))) THEN
              ROTTYP = .TRUE.
C ... Now find block type
              CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', ROTBLK)
            ELSE
              CALL PRTERR('CMDERR',
     *          'Expected keyword "BLOCK" or "TYPE"')
              CALL PRTERR('CMDREQ', STRING(:LENSTR(STRING)))
              GO TO 320
            END IF
 339        CONTINUE
            IF (FFEXST (IFLD, INTYP)) THEN
              CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
              IF (NDIM .EQ. 3 .AND. ((WORD .EQ. 'X') .OR.
     *          (WORD .EQ. 'Y') .OR. (WORD .EQ. 'Z')) ) THEN
                CALL FFREAL (IFLD, INTYP, RFIELD,
     &            'angle of rotation', 0.0, DEG, *320)
                REVATT = .TRUE.
                CALL ROTXYZ (WORD, DEG * DEGANG, ROTATT)
              ELSE IF (NDIM .EQ. 2 .AND. (WORD .EQ. 'Z') ) THEN
                CALL FFREAL (IFLD, INTYP, RFIELD,
     &            'angle of rotation', 0.0, DEG, *320)
                REVATT = .TRUE.
                CALL ROTXYZ (WORD, DEG * DEGANG, ROTATT)
              END IF
              GO TO 339
            END IF

            IF (REVATT) THEN
              WRITE (*, 9940) 'Rotation matrix for selected attributes:'
              DO I = 1, 3
                IX = (I-1) * 3
                DO J = 1, 3
                  RNUM(IX+J) = ROTATT(I,J)
                END DO
              END DO
              CALL NUMSTR (9, 4, RNUM, RSTR, LR)
              DO I = 1, 3
                IX = (I-1) * 3
                WRITE (*, 9940) ('   ', RSTR(IX+J)(:LR), J=1,3)
 9940           FORMAT (1X, 20A)
              END DO
            ELSE
              WRITE (*, 9940) 'No rotation defined for generated mesh'
            END IF
C ... End of "revolve attribute" parsing

          ELSE
            IF (NDIM .EQ. 3) CALL PRTERR ('CMDERR',
     &        'Expected "X", "Y", "Z", "ATTRIBUTE" or "RESET"')
            IF (NDIM .EQ. 2) CALL PRTERR ('CMDERR',
     &        'Expected "Z", "ATTRIBUTE" or "RESET"')
            GOTO 320
          END IF
          GOTO 300
        END IF

 320    CONTINUE

C=======================================================================
      ELSE IF (VERB .EQ. 'REVCEN') THEN
        CALL MINMAX (NUMNP, XN, XMIN, XMAX)
        CALL MINMAX (NUMNP, YN, YMIN, YMAX)
        IF (NDIM .EQ. 3) THEN
          CALL MINMAX (NUMNP, ZN, ZMIN, ZMAX)
        ELSE
          ZMIN = 0.0
        END IF

        CALL FFREAL (IFLD, INTYP, RFIELD,
     &    'X revolution center', XMIN, ROTCEN(1), *420)
        CALL FFREAL (IFLD, INTYP, RFIELD,
     &    'Y revolution center', YMIN, ROTCEN(2), *420)
        IF (NDIM .EQ. 3) CALL FFREAL (IFLD, INTYP, RFIELD,
     &    'Z revolution center', ZMIN, ROTCEN(3), *420)
        IF (NDIM .EQ. 2) ROTCEN(3) = 0.0

C=======================================================================
      ELSE IF (VERB .EQ. 'LIST') THEN
        CALL ABRSTR (LISTYP, CFIELD(2), LISTBL)
        IF (CFIELD(2) .EQ. ' ') THEN
          CALL SHOCMD ('Valid LIST options', LISTBL)
        ELSE IF (LISTYP .EQ. ' ') THEN
          CALL PRTERR ('CMDERR', '"' // CFIELD(2)(:LENSTR(CFIELD(2)))
     &      // '" is an invalid or nonunique LIST option')
          CALL SHOCMD ('Valid LIST options', LISTBL)
        ELSE IF (LISTYP .EQ. 'SSETS' .OR. LISTYP .EQ. 'SIDESETS') THEN
          IF (NUMESS .EQ. 0) THEN
            CALL PRTERR ('CMDSPEC', 'No side sets to list')
          ELSE
            CALL SHOWFL ('S',NUMESS,IDESS,IA(KNESS),IA(KNDSS),SSNAME)
          END IF
        ELSE IF (LISTYP .EQ. 'NSETS' .OR. LISTYP .EQ. 'NODESETS') THEN
          IF (NUMNPS .EQ. 0) THEN
            CALL PRTERR ('CMDSPEC', 'No node sets to list')
          ELSE
            CALL SHOWFL ('N', NUMNPS, IDNPS, IA(KNNNS), IA(1), NSNAME)
          END IF
        ELSE IF (LISTYP .EQ.'VARS' .OR. LISTYP .EQ.'VARIABLE') THEN
          CALL DBPINI ('TISV', NDBIN, TITLE, NDIM, NUMNP, NUMEL,
     &      NELBLK, NUMNPS, LNPSNL, LNPSNL, NUMESS, LESSEL, LESSDF,
     *      NVARGL, NVARNP, NVAREL, NVARNS, NVARSS, ' ')
        ELSE IF (LISTYP .EQ.'BLOCKS' .OR. LISTYP .EQ.'MATERIAL') THEN
          CALL DBPELB ('N', NELBLK, IDELB, NUMELB, NUMLNK, NUMATR,
     &      BLKTYP, EBNAME, ATNAME, IDUM, IDUM, IDUM, IDUM)
        ELSE IF (LISTYP .EQ. 'COMMANDS') THEN
          CALL SHOCMD ('COMMANDS', CMDTBL)
        ELSE IF (LISTYP .EQ. 'QA' .OR. LISTYP .EQ. 'INFORMAT') THEN
          IF ((NQAREC .GT. 0) .OR. (NINFO .GT. 0)) THEN
            CALL DBPQA ('*', NQAREC, QAREC, NINFO, INFREC)
          END IF
        ELSE IF (LISTYP .EQ. 'TIMES' .OR. LISTYP .EQ. 'STEPS') THEN
          if (nsteps .eq. 0) then
            CALL PRTERR ('CMDSPEC', 'No timesteps on database')
          else
            do 324 itim = 1, nsteps
              if (itimst(itim) .eq. 0) then
                STRA = '  Active'
              else
                STRA = 'Inactive'
              end if
              write (*, 10060) itim, times(itim), stra
10060         FORMAT( 'Step ', i6, ', Time = ', 1pe12.5,1x, A)
 324        continue
          end if
        ELSE IF (LISTYP .EQ. 'NAMES') THEN
          WRITE (*, 10070) 'Global:   ', (NAMEGV(I), I=1,NVARGL)
          WRITE (*, 10070) 'Nodal:    ', (NAMENV(I), I=1,NVARNP)
          WRITE (*, 10070) 'Element:  ', (NAMEEV(I), I=1,NVAREL)
          WRITE (*, 10070) 'Nodeset:  ', (NAMEMV(I), I=1,NVARNS)
          WRITE (*, 10070) 'Sideset:  ', (NAMESV(I), I=1,NVARSS)
10070     FORMAT (4X, A, :, 2 (2X, A), :, /, (14X, 2 (2X, A)))
        ELSE
          CALL PRTERR ('CMDERR', '"' // CFIELD(2)(:LENSTR(CFIELD(2)))
     &      // '" is an invalid or nonunique LIST option')
        END IF
        VERB = ' '

      ELSE IF ((VERB .EQ. 'SHOW') .OR. (VERB .EQ. 'LIST')) THEN
        CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
        CALL ABRSTR (VERB, WORD, CMDTBL)
        CALL SHOW (VERB, WORD)
        VERB = ' '

C=======================================================================
C NAME BLOCK|SIDESET|SSET|NODESET|NSET id newname
C NAME ATTRIBUTE attindex BLOCK id newname

      ELSE IF (VERB .EQ. 'NAME') THEN
        CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)

        itype = ' '
        IF (MATSTR(WORD, 'BLOCK', 2)) THEN
          itype = 'B'
        ELSE IF (MATSTR(WORD,'SIDESETS',3) .OR.
     $      MATSTR(WORD,'SSETS',   3) ) THEN
          itype = 'S'
        ELSE IF (MATSTR(WORD,'NODESETS',3) .OR.
     $      MATSTR(WORD,'NSETS',   3) ) THEN
          itype = 'N'
        ELSE IF (MATSTR(WORD,'ATTRIBUTE',3)) THEN
          itype = 'A'
        END IF

        CALL FFINTG (IFLD, INTYP, IFIELD,
     $    'entity ID', 0, IDOLD, *420)
        if (itype .eq. 'A') then
C.. Skip over 'block' and get block id
          CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
          CALL FFINTG (IFLD, INTYP, IFIELD,
     $      'block ID', 0, IDBLK, *420)
        end if
        CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
        if (itype .eq. 'B') then
          CALL NEWNAM ('B', IDELB, EBNAME, NELBLK, IDOLD, WORD)
        else if (itype .eq. 'N') then
          CALL NEWNAM ('N', IDNPS, NSNAME, NUMNPS, IDOLD, WORD)
        else if (itype .eq. 'S') then
          CALL NEWNAM ('S', IDESS, SSNAME, NUMESS, IDOLD, WORD)
        else if (itype .eq. 'A') then
          CALL NEWATT (IDELB, IDBLK,  NELBLK, IDOLD, NUMATR, ATNAME,
     *      WORD)
        end if
        verb = ' '
C=======================================================================
      ELSE IF (VERB .EQ. 'RENAME') THEN
        write (*,*) 'In rename'
        CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
        CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD2)
        if (.NOT. MATSTR (WORD2, 'VARIABLE', 2)) THEN
          call PRTERR('CMDERR',
     *      'SYNTAX: RENAME [global|nodal|node|element|nodeset|nset|'//
     *      'sideset|sset] VARIABLE [oldname] [newname]')
          go to 399
        end if

        CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', OLDNM)
        CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', NEWNM)

        IF (MATSTR(WORD, 'GLOBAL', 2)) THEN
          CALL RENAM(WORD, NAMEGV, NVARGL, OLDNM, NEWNM)
        ELSE IF (MATSTR(WORD,'NODESET',5) .OR.
     $      MATSTR(WORD,'NSET',   3) ) THEN
          CALL RENAM(WORD, NAMEMV, NVARNS, OLDNM, NEWNM)
        ELSE IF (MATSTR(WORD, 'NODAL', 4) .OR.
     *      MATSTR(WORD,'NODE', 4)) THEN
          CALL RENAM(WORD, NAMENV, NVARNP, OLDNM, NEWNM)
        ELSE IF (MATSTR(WORD, 'ELEMENT', 2)) THEN
          CALL RENAM(WORD, NAMEEV, NVAREL, OLDNM, NEWNM)
        ELSE IF (MATSTR(WORD,'SIDESET',3) .OR.
     $      MATSTR(WORD,'SSET',   3) ) THEN
          CALL RENAM(WORD, NAMESV, NVARSS, OLDNM, NEWNM)
        END IF

 399    CONTINUE
        verb = ' '
C=======================================================================
      ELSE IF (VERB .EQ. 'CHANGE') THEN
        CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)

        IF (MATSTR(WORD, 'TYPE', 2)) THEN
          CALL FFINTG (IFLD, INTYP, IFIELD,
     $      'block ID', 0, IDOLD, *420)
          CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
          CALL NEWNAM ('b', IDELB, BLKTYP, NELBLK, IDOLD, WORD)
        ELSE IF (MATSTR(WORD,'SIDESETS',3) .OR.
     $      MATSTR(WORD,'SSETS',   3) ) THEN
C========================================================================
C ... CHANGE SIDESET {id|ALL} FACTOR {value}
C ... See if next field is integer (id) or word ('all')
          VERB = ' '
          STRING = 'CHANGE SIDESET {id|ALL} FACTOR {value}'
          IF (NUMESS .EQ. 0) THEN
            CALL PRTERR ('CMDERR',
     *        'There are no sidesets on this database')
            GO TO 420
          END IF

          DOALLA = .FALSE.

          IF ((INTYP(IFLD) .EQ. 0) .OR. ((CFIELD(IFLD) .GE. 'A')
     &      .AND. (CFIELD(IFLD) .LE. 'Z'))) THEN
C ... Field is character
            CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
            IF (MATSTR (WORD, 'ALL', 2)) THEN
              DOALLA = .TRUE.
            ELSE
              CALL PRTERR('CMDERR',
     *          'Expected keyword "ALL"')
              CALL PRTERR('CMDREQ', STRING(:LENSTR(STRING)))
              GO TO 420
            END IF
          ELSE
C ... Field is sideset id
            CALL FFINTG (IFLD, INTYP, IFIELD,
     $        'sideset ID', 0, IDSET, *420)
          END IF

C ... Next word should be 'FACTOR' or an id
          IF ((INTYP(IFLD) .EQ. 0) .OR. ((CFIELD(IFLD) .GE. 'A')
     &      .AND. (CFIELD(IFLD) .LE. 'Z'))) THEN
C     ... Field is character
            CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
            IF (.NOT.(MATSTR (WORD, 'FACTOR', 1))) THEN
              CALL PRTERR('CMDERR',
     *          'Expected keyword "FACTOR"')
              CALL PRTERR('CMDREQ', STRING(:LENSTR(STRING)))
              GO TO 420
            END IF

            CALL FFREAL (IFLD, INTYP, RFIELD,
     *        'distribution factor value', 1.0, AVALUE, *420)

C     ... Setup the scale and echo to the user
            if (DOALLA) then
              CALL INIREA (LESSNL, AVALUE, FACESS)
            ELSE
              iss = locint(idset, numess, idess)
              IF (iss .EQ. 0) THEN
                WRITE (STRING, 10040) 'Sideset', idset
10040           FORMAT (A,I5,' does not exist')
                CALL SQZSTR (STRING, LSTR)
                CALL PRTERR ('ERROR', STRING(:LSTR))
                goto 420
              END IF
              call inirea(nness(iss), avalue, facess(ixness(iss)))
            ENDIF
          ELSE
            CALL FFINTG (IFLD, INTYP, IFIELD,'new ID', 0, IDNEW, *420)
            CALL NEWID ('S', IDESS, NUMESS, IDNEW, IDSET)
          END IF

        ELSE IF (MATSTR(WORD,'NODESETS',3) .OR.
     $      MATSTR(WORD,'NSETS',   3) ) THEN
C========================================================================
C ... CHANGE NODESET {id|ALL} FACTOR {value}
C ... See if next field is integer (id) or word ('all')
          VERB = ' '
          STRING = 'CHANGE NODESET {id|ALL} FACTOR {value}'
          IF (NUMNPS .EQ. 0) THEN
            CALL PRTERR ('CMDERR',
     *        'There are no nodesets on this database')
            GO TO 420
          END IF

          DOALLA = .FALSE.

          IF ((INTYP(IFLD) .EQ. 0) .OR. ((CFIELD(IFLD) .GE. 'A')
     &      .AND. (CFIELD(IFLD) .LE. 'Z'))) THEN
C ... Field is character
            CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
            IF (MATSTR (WORD, 'ALL', 2)) THEN
              DOALLA = .TRUE.
            ELSE
              CALL PRTERR('CMDERR',
     *          'Expected keyword "ALL"')
              CALL PRTERR('CMDREQ', STRING(:LENSTR(STRING)))
              GO TO 420
            END IF
          ELSE
C ... Field is nodeset id
            CALL FFINTG (IFLD, INTYP, IFIELD,
     $        'nodeset ID', 0, IDSET, *420)
          END IF

C ... Next word should be 'FACTOR' or an id
          IF ((INTYP(IFLD) .EQ. 0) .OR. ((CFIELD(IFLD) .GE. 'A')
     &      .AND. (CFIELD(IFLD) .LE. 'Z'))) THEN
C     ... Field is character
            CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
            IF (.NOT.(MATSTR (WORD, 'FACTOR', 1))) THEN
              CALL PRTERR('CMDERR',
     *          'Expected keyword "FACTOR"')
              CALL PRTERR('CMDREQ', STRING(:LENSTR(STRING)))
              GO TO 420
            END IF

            CALL FFREAL (IFLD, INTYP, RFIELD,
     *        'distribution factor value', 1.0, AVALUE, *420)

C     ... Setup the scale and echo to the user
            if (DOALLA) then
              CALL INIREA (LNPSNL, AVALUE, FACNPS)
            ELSE
              iss = locint(idset, numnps, idnps)
              IF (iss .EQ. 0) THEN
                WRITE (STRING, 10040) 'Nodeset', idset
                CALL SQZSTR (STRING, LSTR)
                CALL PRTERR ('ERROR', STRING(:LSTR))
                goto 420
              END IF
              call inirea(ndnps(iss), avalue, facnps(ixdnps(iss)))
            ENDIF
          ELSE
C     ... Changing id. Get the new id.
            CALL FFINTG (IFLD, INTYP, IFIELD, 'new ID', 0, IDNEW, *420)
            CALL NEWID ('N', IDNPS, NUMNPS, IDNEW, IDSET)
          END IF

        ELSE IF (MATSTR(WORD, 'ATTRIBUTE', 3)) THEN
C========================================================================
C ... CHANGE ATTRIBUTE {num|ALL} BLOCK {id|ALL} {value}
C ... CHANGE ATTRIBUTE RESET
C ... See if next field is integer (id) or word ('all')
          VERB = ' '
          STRING = 'CHANGE ATTRIBUTE {id|ALL} BLOCK {id|ALL} {value}'
          IF (.NOT. ISATRB) THEN
            CALL PRTERR ('CMDERR',
     *        'There are no attributes on this database')
            GO TO 420
          END IF

          ATRWRN = .TRUE.
          DOALLA = .FALSE.
          DOALLB = .FALSE.

          IF ((INTYP(IFLD) .EQ. 0) .OR. ((CFIELD(IFLD) .GE. 'A')
     &      .AND. (CFIELD(IFLD) .LE. 'Z'))) THEN
C ... Field is character
            CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
            IF (MATSTR (WORD, 'ALL', 2)) THEN
              DOALLA = .TRUE.
            ELSE IF (MATSTR (WORD, 'RESET', 2)) THEN
              CALL SCLATR(NELBLK, IDELB, NUMATR, 0, .TRUE., 0,
     *          .TRUE., 0.0, ATRSCL, .TRUE., .FALSE.)
              GO TO 420
            ELSE
              CALL PRTERR('CMDERR',
     *          'Expected keyword "ALL"')
              CALL PRTERR('CMDREQ', STRING(:LENSTR(STRING)))
              GO TO 420
            END IF
          ELSE
C ... Field is Attribute number
            CALL FFINTG (IFLD, INTYP, IFIELD,
     $        'attribute ID', 0, IDATR, *420)
          END IF

C ... Next word should be 'BLOCK'
          CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
          IF (.NOT.(MATSTR (WORD, 'BLOCK', 1))) THEN
            CALL PRTERR('CMDERR',
     *        'Expected keyword "BLOCK"')
            CALL PRTERR('CMDREQ', STRING(:LENSTR(STRING)))
            GO TO 420
          END IF
C ... Now find block id or 'ALL'
          IF ((INTYP(IFLD) .EQ. 0) .OR. ((CFIELD(IFLD) .GE. 'A')
     &      .AND. (CFIELD(IFLD) .LE. 'Z'))) THEN
C ... Field is character
            CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
            IF (MATSTR (WORD, 'ALL', 2)) THEN
              DOALLB = .TRUE.
            ELSE
              CALL PRTERR('CMDERR',
     *          'Expected keyword "ALL"')
              CALL PRTERR('CMDREQ', STRING(:LENSTR(STRING)))
              GO TO 420
            END IF
          ELSE
C ... Field is Attribute number
            CALL FFINTG (IFLD, INTYP, IFIELD,
     $        'block ID', 0, IDBLK, *420)
          END IF

          CALL FFREAL (IFLD, INTYP, RFIELD,
     *      'attribute value', 1.0, AVALUE, *420)

C ... Setup the scale and echo to the user
          CALL SCLATR(NELBLK, IDELB, NUMATR, IDATR, DOALLA,
     *      IDBLK, DOALLB, AVALUE, ATRSCL, .TRUE., .FALSE.)

        ELSE
          CALL FFINTG (IFLD, INTYP, IFIELD,
     *      'old ID', 0, IDOLD, *420)
          CALL FFINTG (IFLD, INTYP, IFIELD,
     *      'new ID', 0, IDNEW, *420)

          IF (MATSTR(WORD,'MATERIALS',1) .OR.
     &      MATSTR(WORD, 'BLOCKS',1)) THEN
            CALL NEWID ('M', IDELB, NELBLK, IDNEW, IDOLD)
          ELSE IF (MATSTR(WORD,'NODESETS',1) .OR.
     $        MATSTR(WORD,'NSETS',   1) ) THEN
            CALL NEWID ('N', IDNPS, NUMNPS, IDNEW, IDOLD)
          ELSE IF (MATSTR(WORD,'SIDESETS',1) .OR.
     $        MATSTR(WORD,'SSETS',   1)) THEN
            CALL NEWID ('S', IDESS, NUMESS, IDNEW, IDOLD)
          END IF
        END IF
        VERB = ' '

C=======================================================================
      ELSE IF (VERB .EQ. 'INCREMENT') THEN
        CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
        CALL FFINTG (IFLD, INTYP, IFIELD,
     *    'increment', 0, IDINC, *420)

        IF (MATSTR(WORD,'MATERIALS',1) .OR.
     &    MATSTR(WORD, 'BLOCKS',1)) THEN
          CALL INCID ('M', IDELB, NELBLK, IDINC)
        ELSE IF (MATSTR(WORD,'NODESETS',1) .OR.
     $      MATSTR(WORD,'NSETS',   1)) THEN
          CALL INCID ('N', IDNPS, NUMNPS, IDINC)
        ELSE IF (MATSTR(WORD,'SIDESETS',1) .OR.
     $      MATSTR(WORD,'SSETS',   1)) THEN
          CALL INCID ('S', IDESS, NUMESS, IDINC)
        END IF
        VERB = ' '

C=======================================================================
      else if (verb .eq. 'MERGE') THEN
C     ... Syntax: MERGE {global|local} Node {node_1} into node {node_2}
        CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
        if  (matstr(word, 'GLOBAL', 3)) then
          glonod = .true.
        else if  (matstr(word, 'LOCAL', 3)) then
          glonod = .false.
        else
          CALL PRTERR ('CMDERR',
     *      'Expected "global" or "local" following "merge"')
          go to 420
        end if
c ... skip 'node'
        CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)

C ... Get the two nodes to be equivalenced...
        CALL FFINTG (IFLD, INTYP, IFIELD, 'Node1', 0, node1, *420)
C ... Skip 'node into'
        CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
        CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
        CALL FFINTG (IFLD, INTYP, IFIELD, 'Node2', 0, node2, *420)

        if (glonod) then
          node1 = ng2l('node', node1, numnp, mapnod)
          node2 = ng2l('node', node2, numnp, mapnod)
        end if

        if (node1 .gt. numnp .or. node2 .gt. numnp) then
          CALL PRTERR ('CMDERR', 'Node index exceeds node count')
        else
C ... Grepos has requirement that the highest numbered node be removed
C     in the equiv process.  However, in user specification, they want
C     to specify which node is removed.  If it is the lower-numbered node,
C     simply replace the coordinates of the higher numbered node with the
C     coordinates of the lower-numbered node.
          if (node1 .lt. node2) then
            if (IEQUIV(node2) .ne. node2) then
                call prterr ('CMDERR',
     *              'Node already merged, cannot specify twice')
            else
                IEQUIV(node2) = -node1
                xn(node1) = xn(node2)
                if (ndim .ge. 2) yn(node1) = yn(node2)
                if (ndim .eq. 3) zn(node1) = zn(node2)
            endif
          else
            if (IEQUIV(node1) .ne. node1) then
                call prterr ('CMDERR',
     *              'Node already merged, cannot specify twice')
            else
                IEQUIV(node1) = -node2
            end if
          end if
          equiv = .TRUE.
          EQTOLER = -1.0
        end if
        verb = ' '

      ELSE IF (VERB .EQ. 'EQUIVALENCE') THEN
C     ... Syntax: EQUIVALENCE tolerance
C     ... Syntax: EQUIVALENCE RESET

        IF ((INTYP(IFLD) .EQ. 0) .OR. ((CFIELD(IFLD) .GE. 'A')
     &    .AND. (CFIELD(IFLD) .LE. 'Z'))) THEN
C     ... Field is character
          CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
          IF (MATSTR (WORD, 'RESET', 3)) THEN
            EQUIV = .FALSE.
            EQTOLER = 0.0
            do i=1, numnp
              iequiv(i) = i
            end do
          ELSE
            CALL PRTERR ('CMDERR', 'Expected "RESET"')
          END IF
          verb = ' '
        ELSE
          EQUIV = .TRUE.
          CALL FFREAL (IFLD, INTYP, RFIELD,
     &      'Tolerance', 0.0, EQTOLER, *420)
          IF (EQTOLER .LT. 0.0) THEN
            CALL PRTERR ('CMDERR',
     &        "Tolerance must be greater than 0.0.")
            EQUIV = .FALSE.
          END IF
        END IF

C=======================================================================
      ELSE IF (VERB .EQ. 'COMBINE') THEN
C ... Syntax: COMBINE entity   RESET (Not Implemented yet)
C ... Syntax: COMBINE BLOCK    id_final WITH BLOCK    id1, id2, ..., idn
C ... Syntax: COMBINE MATERIAL id_final WITH MATERIAL id1, id2, ..., idn
C ... Syntax: COMBINE NODESET  id_final WITH NODESET  id1, id2, ..., idn
C ... Syntax: COMBINE SIDESET  id_final WITH SIDESET  id1, id2, ..., idn
C     Sets ids of id1, id2, ..., idn to -idfinal
        CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
        IF (FFEXST (IFLD, INTYP)) THEN
          CALL FFINTG (IFLD, INTYP, IFIELD,
     *      'ID', 0, IDFINL, *420)
        ELSE
          CALL PRTERR ('CMDERR',
     $      "Missing Keyword 'MATERIALS/NODESETS/SIDESETS'")
          GO TO 420
        END IF

C ... This should be the keyword 'WITH'
        CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD2)
        if (.NOT. MATSTR(WORD2, 'WITH', 1)) THEN
          CALL PRTERR('CMDERR',
     $      'Expected keyword "WITH", found "' //
     $      WORD2(:lenstr(word2)) // '".')
          go to 420
        end if

C ... This should be the entity type. Must match WORD.
C     Matching assumed, not checked.
        CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD2)
        if (.NOT. MATSTR(WORD2, WORD, 1)) THEN
          CALL PRTERR('CMDERR',
     $      'Expected keyword "'//word(:lenstr(word))//'" found "' //
     $      WORD2(:lenstr(word2)) // '".')
          go to 420
        end if

 325    CONTINUE
        IF (FFEXST (IFLD, INTYP)) THEN
          CALL FFINTG (IFLD, INTYP, IFIELD,
     *      'ID', 0, ID, *420)
          IF (MATSTR(WORD,'MATERIALS',1) .OR.
     &      MATSTR(WORD, 'BLOCKS',1)) THEN
            CALL COMID ('Material', IDELB, IELBST, NELBLK, IDFINL, ID)
          ELSE IF (MATSTR(WORD,'NODESETS',1) .OR.
     &        MATSTR(WORD,'NSETS',   1) ) THEN
            CALL COMID ('Nodeset', IDNPS, INPSST, NUMNPS, IDFINL, ID)
          ELSE IF (MATSTR(WORD,'SIDESETS',1) .OR.
     &        MATSTR(WORD,'SSETS',   1) ) THEN
            CALL COMID ('Sideset', IDESS, IESSST, NUMESS, IDFINL, ID)
          ELSE
            write (string, 326) word(:lenstr(word))
 326        format('"',A,
     $        '" is an invalid or nonunique COMBINE option')
            call sqzstr(string, lstr)
            CALL PRTERR ('CMDERR', string(:lstr))
          END IF
          GOTO 325
        END IF
        VERB = ' '
C========================================================================
      ELSE IF (VERB .EQ. 'DELETE') THEN
        CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
        IF (MATSTR(WORD,'QAINFO',3)) THEN
          NQAREC = 0
          NINFO  = 1
          INFREC(NINFO) =
     $      'All previous QA and INFO records deleted by GREPOS'
          CALL PRTERR ('CMDSPEC',
     $      'All QA and INFO records deleted.')
        ELSE IF (MATSTR(WORD,'QA',2)) THEN
          NQAREC = 0
          NINFO  = NINFO + 1
          INFREC(NINFO) =
     $      'All previous QA records deleted by GREPOS'
          CALL PRTERR ('CMDSPEC',
     $      'All QA records deleted.')
        ELSE IF (MATSTR(WORD,'INFORMATION',2)) THEN
          NINFO  = 1
          INFREC(NINFO) =
     $      'All previous INFO records deleted by GREPOS'
          CALL PRTERR ('CMDSPEC',
     $      'All INFO records deleted.')
        ELSE
          IF (FFEXST (IFLD, INTYP)) THEN
            IF ((INTYP(IFLD) .EQ. 0) .OR. ((CFIELD(IFLD) .GE. 'A')
     &        .AND. (CFIELD(IFLD) .LE. 'Z'))) THEN
C ... Field is character
              CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD2)
              IF (MATSTR (WORD2, 'MAP', 3)) THEN
                IF (MATSTR (WORD, 'NODE', 2)) THEN
                  do i=1, numnp
                    mapnod(i) = i
                  end do
                ELSE IF (MATSTR (WORD, 'ELEMENT', 2)) THEN
                  do i=1, numel
                    mapel(i) = i
                  end do
                ELSE
                  CALL PRTERR ('CMDERR', '"' // WORD(:LENSTR(WORD))
     &              // '" is an invalid DELETE MAP -- NODE or ELEMENT')
                  GO TO 333
                END IF
              ELSE IF (MATSTR (WORD2, 'ALL', 2)) THEN
                IF (MATSTR(WORD,'MATERIALS',1) .OR.
     &            MATSTR(WORD, 'BLOCKS',1)) THEN
                  do i=1, nelblk
                    ielbst(i) = -idelb(i)
                  end do
                  CALL PRTERR ('CMDSPEC','All material blocks deleted.')
                ELSE IF (MATSTR(WORD,'NODESETS',1) .OR.
     &              MATSTR(WORD,'NSETS',   1) ) THEN
                  do i=1, numnps
                    inpsst(i) = -idnps(i)
                  end do
                  CALL PRTERR ('CMDSPEC','All nodesets deleted.')
                ELSE IF (MATSTR(WORD,'SIDESETS',1) .OR.
     &              MATSTR(WORD,'SSETS',   1) ) THEN
                  do i=1, numess
                    iessst(i) = -idess(i)
                  end do
                  CALL PRTERR ('CMDSPEC','All sidesets deleted.')
                ELSE
                  CALL PRTERR ('CMDERR', '"' // WORD(:LENSTR(WORD))
     &              // '" is an invalid or nonunique DELETE option')
                  GO TO 333
                END IF

              ELSE IF (MATSTR (WORD2, 'NAMES', 2)) THEN
                IF (MATSTR(WORD,'MATERIALS',1) .OR.
     &            MATSTR(WORD, 'BLOCKS',1)) THEN
                  do i=1, nelblk
#if defined(__XLF__)
                    ebname(i) = " "
#else
                    ebname(i) = ""
#endif
                  end do
                  CALL PRTERR ('CMDSPEC',
     *              'All material block names deleted.')
                ELSE IF (MATSTR(WORD,'NODESETS',1) .OR.
     &              MATSTR(WORD,'NSETS',   1) ) THEN
                  do i=1, numnps
#if defined(__XLF__)
                    nsname(i) = " "
#else
                    nsname(i) = ""
#endif
                  end do
                  CALL PRTERR ('CMDSPEC','All nodeset names deleted.')
                ELSE IF (MATSTR(WORD,'SIDESETS',1) .OR.
     &              MATSTR(WORD,'SSETS',   1) ) THEN
                  do i=1, numess
#if defined(__XLF__)
                    ssname(i) = " "
#else
                    ssname(i) = ""
#endif
                  end do
                  CALL PRTERR ('CMDSPEC','All sideset names deleted.')
                ELSE
                  CALL PRTERR ('CMDERR', '"' // WORD(:LENSTR(WORD))
     &              // '" is an invalid DELETE NAME option')
                  GO TO 333
                END IF

              ELSE IF (MATSTR (WORD2, 'RESET', 2)) THEN
                IF (MATSTR(WORD,'MATERIALS',1) .OR.
     &            MATSTR(WORD, 'BLOCKS',1)) THEN
                  do i=1, nelblk
                    ielbst(i) = 0
                  end do
                  CALL PRTERR ('CMDSPEC',
     *              'All material blocks undeleted.')
                ELSE IF (MATSTR(WORD,'NODESETS',1) .OR.
     &              MATSTR(WORD,'NSETS',   1) ) THEN
                  do i=1, numnps
                    inpsst(i) = 0
                  end do
                  CALL PRTERR ('CMDSPEC','All nodesets undeleted.')
                ELSE IF (MATSTR(WORD,'SIDESETS',1) .OR.
     &              MATSTR(WORD,'SSETS',   1) ) THEN
                  do i=1, numess
                    iessst(i) = 0
                  end do
                  CALL PRTERR ('CMDSPEC','All sidesets undeleted.')
                ELSE
                  CALL PRTERR ('CMDERR', '"' // WORD(:LENSTR(WORD))
     &              // '" is an invalid or nonunique DELETE option')
                  GO TO 333
                END IF

              ELSE
                CALL PRTERR('CMDERR',
     *            'Expected keyword "ALL" or "RESET"')
                CALL PRTERR('CMDREQ', STRING(:LENSTR(STRING)))
                GO TO 333
              END IF
            ENDIF
          END IF

C ... Ids only.  Either a list or range (or both)
 330      CONTINUE
          IF (FFEXST (IFLD, INTYP)) THEN

C         --Scan numeric range

            CALL FFVRNG (IFLD, INTYP, CFIELD, IFIELD,
     &        WORD, -1, IRNG, *333)

C         --Store the range selected

            DO 332 I = IRNG(1), IRNG(2), IRNG(3)
              IF (MATSTR(WORD,'MATERIALS',1) .OR.
     &          MATSTR(WORD, 'BLOCKS',1)) THEN
                CALL DELID ('M', IDELB, IELBST, NELBLK, I, .TRUE.)
              ELSE IF (MATSTR(WORD,'NODESETS',1) .OR.
     &            MATSTR(WORD,'NSETS',   1) ) THEN
                CALL DELID ('N', IDNPS, INPSST, NUMNPS, I, .TRUE.)
              ELSE IF (MATSTR(WORD,'SIDESETS',2) .OR.
     &            MATSTR(WORD,'SSETS',   1) ) THEN
                CALL DELID ('S', IDESS, IESSST, NUMESS, I, .TRUE.)
              ELSE IF (MATSTR(WORD,'STEPS',2) .OR.
     &            MATSTR(WORD,'TIMES',   1) ) THEN
                if (i .le. nsteps) then
                  ITIMST(I) = 1
                else
                  CALL PRTERR ('CMDERR',
     *              'Step number is out of range')
                end if
              ELSE
                CALL PRTERR ('CMDERR', '"' // WORD(:LENSTR(WORD))
     &            // '" is an invalid or nonunique DELETE option')
              END IF
 332        CONTINUE

            GOTO 330
          END IF
        END IF
 333    CONTINUE
        VERB = ' '

C========================================================================
      ELSE IF (VERB .EQ. 'UNDELETE' .OR. VERB .EQ. 'KEEP') THEN
        CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
        IF ((INTYP(IFLD) .EQ. 0) .OR. ((CFIELD(IFLD) .GE. 'A')
     &    .AND. (CFIELD(IFLD) .LE. 'Z'))) THEN
C ... Field is character
          CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD2)
          IF (MATSTR (WORD2, 'ALL', 2)) THEN
            IF (MATSTR(WORD,'MATERIALS',1) .OR.
     &        MATSTR(WORD, 'BLOCKS',1)) THEN
              do i=1, nelblk
                ielbst(i) = 0
              end do
              CALL PRTERR ('CMDSPEC','All material blocks undeleted.')
            ELSE IF (MATSTR(WORD,'NODESETS',1) .OR.
     &          MATSTR(WORD,'NSETS',   1) ) THEN
              do i=1, numnps
                inpsst(i) = 0
              end do
              CALL PRTERR ('CMDSPEC','All nodesets undeleted.')
            ELSE IF (MATSTR(WORD,'SIDESETS',1) .OR.
     &          MATSTR(WORD,'SSETS',   1) ) THEN
              do i=1, numess
                iessst(i) = 0
              end do
              CALL PRTERR ('CMDSPEC','All sidesets undeleted.')
            ELSE
              CALL PRTERR ('CMDERR', '"' // WORD(:LENSTR(WORD))
     &          // '" is an invalid or nonunique DELETE option')
              GO TO 336
            END IF

          ELSE
            CALL PRTERR('CMDERR',
     *        'Expected keyword "ALL"')
            CALL PRTERR('CMDREQ', STRING(:LENSTR(STRING)))
            GO TO 336
          END IF
        ELSE
 335      CONTINUE
          IF (FFEXST (IFLD, INTYP)) THEN

C         --Scan numeric range
            CALL FFVRNG (IFLD, INTYP, CFIELD, IFIELD,
     &        WORD, -1, IRNG, *336)

C         --Store the range selected
            DO I = IRNG(1), IRNG(2), IRNG(3)
              IF (MATSTR(WORD,'MATERIALS',1) .OR.
     &          MATSTR(WORD, 'BLOCKS',1)) THEN
                CALL DELID ('M', IDELB, IELBST, NELBLK, I, .FALSE.)
              ELSE IF (MATSTR(WORD,'NODESETS',1) .OR.
     &            MATSTR(WORD,'NSETS',   1) ) THEN
                CALL DELID ('N', IDNPS, INPSST, NUMNPS, I, .FALSE.)
              ELSE IF (MATSTR(WORD,'SIDESETS',2) .OR.
     &            MATSTR(WORD,'SSETS',   1) ) THEN
                CALL DELID ('S', IDESS, IESSST, NUMESS, I, .FALSE.)
              ELSE IF (MATSTR(WORD,'STEPS',2) .OR.
     &            MATSTR(WORD,'TIMES',   1) ) THEN
                if (i .le. nsteps) then
                  ITIMST(I) = 0
                else
                  CALL PRTERR ('CMDERR',
     *              'Step number is out of range')
                end if
              ELSE
                CALL PRTERR ('CMDERR', '"' // WORD(:LENSTR(WORD))
     &            // '" is an invalid or nonunique DELETE option')
              END IF
            END DO

            GOTO 335
          END IF
        END IF
 336    CONTINUE
        VERB = ' '

C========================================================================
      ELSE IF (VERB .EQ. 'SWAP') THEN
        CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
        IF (MATSTR(WORD,'SIDESETS',1) .OR.
     &    MATSTR(WORD,'SSETS',   1) ) THEN
 340      CONTINUE
          IF (FFEXST (IFLD, INTYP)) THEN
            CALL FFINTG (IFLD, INTYP, IFIELD,
     *        'ID', 0, IDNEW, *420)
            CALL SWPID ('S', IDESS, NUMESS, IDNEW, DELOK)
            IF (DELOK) SWPSS = .TRUE.
            GOTO 340
          END IF
        ELSE
          CALL PRTERR ('CMDERR', '"' // WORD(:LENSTR(WORD))
     &      // '" is an invalid or nonunique SWAP option')
        END IF
        VERB = ' '

C========================================================================
      ELSE IF (VERB .EQ. 'TMIN') then
         call ffreal (ifld, intyp, rfield,
     $    'Minimum Time', 0.0, tmin, *345)
         do itim = 1, nsteps
            if (times(itim) .lt. tmin) then
               itimst(itim) = 1
            end if
         end do
 345     continue
        VERB = ' '
C========================================================================
      ELSE IF (VERB .EQ. 'TMAX') then
         call ffreal (ifld, intyp, rfield,
     $    'Maximum Time', 0.0, tmax, *346)
         do itim = 1, nsteps
            if (times(itim) .gt. tmax) then
               itimst(itim) = 1
            end if
         end do
 346     continue
        VERB = ' '
C========================================================================

      ELSE IF (VERB .EQ. 'SMOOTH') then
        SMOOTH = .TRUE.
        call ffreal (ifld, intyp, rfield,
     $    'Tolerance', 1.0e-6, toler, *350)
        call ffintg (ifld, intyp, ifield,
     $    'Iterations', 100, NIT, *350)
        call ffreal (ifld, intyp, rfield,
     $    'Relaxation Factor', 1.0, R0, *350)

 350    continue

C ... SNAP
C========================================================================
      ELSE IF (VERB .EQ. 'SNAP') then
        if (ffexst(ifld, intyp)) then
          IF ((INTYP(IFLD) .EQ. 0) .OR. ((CFIELD(IFLD) .GE. 'A')
     &      .AND. (CFIELD(IFLD) .LE. 'Z'))) THEN
            CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
            IF (MATSTR (WORD, 'RESET', 1)) THEN
              NUMSNP = 0
              go to 360
            ELSE
              call PRTERR('CMDERR', 'Expected "RESET" or sideset id')
              go to 370
            END IF
          END IF

          if (ffexst(ifld, intyp)) then
            numsnp = numsnp + 1
            if (numsnp .gt. mxsnap) then
              call PRTERR('CMDERR', 'Too many snap sets specified')
              go to 370
            end if

            ismtyp(numsnp) = ISNAP
            usnorm(numsnp) = PNORM
            snptol(numsnp) = 0.0
            gap(numsnp)    = 0.0
C ... NOTE: This value is squared in snpnod, don't increase too big.
            delmax(numsnp) = 1.0e15

            call ffintg (ifld, intyp, ifield,
     $        'sideset to snap nodes from', 0, IDSSSL(numsnp), *370)
            imat = locint(idsssl(numsnp), numess, idess)
            IF (IMAT .EQ. 0) THEN
              WRITE (STRING, 10040) 'Sideset', idsssl(numsnp)
              CALL SQZSTR (STRING, LSTR)
              CALL PRTERR ('ERROR', STRING(:LSTR))
              NUMSNP = NUMSNP - 1
              GO TO 370
            END IF

            CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
            IF (.NOT. MATSTR (WORD, 'TO', 1)) THEN
              CALL PRTERR('CMDERR',
     *          'Expected required keyword "TO"')
              NUMSNP = NUMSNP - 1
              GO TO 370
            END IF

            call ffintg (ifld, intyp, ifield,
     $        'sideset to snap nodes to', 0, IDSSMA(numsnp), *370)
            imat = locint(idssma(numsnp), numess, idess)
            IF (IMAT .EQ. 0) THEN
              WRITE (STRING, 10040) 'Sideset', idssma(numsnp)
              CALL SQZSTR (STRING, LSTR)
              CALL PRTERR ('ERROR', STRING(:LSTR))
              NUMSNP = NUMSNP - 1
              GO TO 370
            END IF
          end if

 360      CONTINUE
          IF (FFEXST (IFLD, INTYP)) THEN
            CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
            IF (MATSTR (WORD, 'RESET', 1)) THEN
              numsnp = 0
              go to 360
            else IF (MATSTR(WORD, 'NORMAL', 1)) THEN
              USNORM(NUMSNP) = PNORM
              VECTOR(1,numsnp) = 0.0
              VECTOR(2,numsnp) = 0.0
              VECTOR(3,numsnp) = 0.0
            else IF (MATSTR(WORD, 'EDGE', 1)) THEN
              USNORM(NUMSNP) = PEDGE
              VECTOR(1,numsnp) = 0.0
              VECTOR(2,numsnp) = 0.0
              VECTOR(3,numsnp) = 0.0
            ELSE IF (WORD .EQ. 'X') THEN
              USNORM(NUMSNP) = PVECT
              VECTOR(1,NUMSNP) = 1.0
              VECTOR(2,NUMSNP) = 0.0
              VECTOR(3,NUMSNP) = 0.0
            ELSE IF (WORD .EQ. 'Y') THEN
              USNORM(NUMSNP) = PVECT
              VECTOR(1,NUMSNP) = 0.0
              VECTOR(2,NUMSNP) = 1.0
              VECTOR(3,NUMSNP) = 0.0
            ELSE IF (WORD .EQ. 'Z') THEN
              USNORM(NUMSNP) = PVECT
              VECTOR(1,NUMSNP) = 0.0
              VECTOR(2,NUMSNP) = 0.0
              VECTOR(3,NUMSNP) = 1.0
            ELSE IF (WORD .EQ. 'MINUSX') THEN
              USNORM(NUMSNP) = PVECT
              VECTOR(1,NUMSNP) = -1.0
              VECTOR(2,NUMSNP) =  0.0
              VECTOR(3,NUMSNP) =  0.0
            ELSE IF (WORD .EQ. 'MINUSY') THEN
              USNORM(NUMSNP) = PVECT
              VECTOR(1,NUMSNP) =  0.0
              VECTOR(2,NUMSNP) = -1.0
              VECTOR(3,NUMSNP) =  0.0
            ELSE IF (WORD .EQ. 'MINUSZ') THEN
              USNORM(NUMSNP) = PVECT
              VECTOR(1,NUMSNP) =  0.0
              VECTOR(2,NUMSNP) =  0.0
              VECTOR(3,NUMSNP) = -1.0
            ELSE IF (MATSTR(WORD, 'VECTOR', 1)) THEN
              USNORM(NUMSNP) = PVECT
              CALL FFREAL (IFLD, INTYP, RFIELD,
     &          'vector X coordinate', 0.0, VECTOR(1,NUMSNP), *370)
              CALL FFREAL (IFLD, INTYP, RFIELD,
     &          'vector Y coordinate', 0.0, VECTOR(2,NUMSNP), *370)
              CALL FFREAL (IFLD, INTYP, RFIELD,
     &          'vector Z coordinate', 0.0, VECTOR(3,NUMSNP), *370)
            ELSE IF (MATSTR(WORD, 'RADIAL', 1)) THEN
              USNORM(NUMSNP) = PRAD
              VECTOR(1,NUMSNP) = 0.0
              VECTOR(2,NUMSNP) = 0.0
              VECTOR(3,NUMSNP) = 0.0
            ELSE IF (MATSTR(WORD, 'CENTER', 1)) THEN
              USNORM(NUMSNP) = PRAD
              CALL FFREAL (IFLD, INTYP, RFIELD,
     &          'X center', 0.0, VECTOR(1,NUMSNP), *370)
              CALL FFREAL (IFLD, INTYP, RFIELD,
     &          'Y center', 0.0, VECTOR(2,NUMSNP), *370)
              CALL FFREAL (IFLD, INTYP, RFIELD,
     &          'Z center', 0.0, VECTOR(3,NUMSNP), *370)
            ELSE IF (MATSTR(WORD, 'TOLERANCE', 1)) THEN
              CALL FFREAL (IFLD, INTYP, RFIELD,
     &          'Tolerance', 0.0, SNPTOL(NUMSNP), *370)
            ELSE IF (MATSTR(WORD, 'MAXDELTA', 1)) THEN
              CALL FFREAL (IFLD, INTYP, RFIELD,
     &          'MaxDelta', 1.0E15, DELMAX(NUMSNP), *370)
            ELSE IF (MATSTR(WORD, 'GAP', 1)) THEN
              CALL PRTERR('CMDERR',
     *          '"GAP" not supported for "SNAP" command')
            END IF
            GO TO 360
          END IF
        end if
 370    CONTINUE
        VERB = 'SNAP'

C ... MOVE
C========================================================================
      ELSE IF (VERB .EQ. 'MOVE') then
        if (ffexst(ifld, intyp)) then
          IF ((INTYP(IFLD) .EQ. 0) .OR. ((CFIELD(IFLD) .GE. 'A')
     &      .AND. (CFIELD(IFLD) .LE. 'Z'))) THEN
            CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
            IF (MATSTR (WORD, 'RESET', 1)) THEN
              NUMSNP = 0
              go to 380
            ELSE IF (MATSTR (WORD, 'GLOBAL', 3) .OR.
     *          MATSTR (WORD, 'LOCAL', 3)) THEN
C ...           Syntax here is:
C                   MOVE GLOBAL NODE {id} TO {x} {y} {z} (absolute move)
C                   MOVE LOCAL  NODE {id} TO {x} {y} {z}
C                   MOVE GLOBAL NODE {id} BY {x} {y} {z} (relative move)
C                   MOVE LOCAL  NODE {id} BY {x} {y} {z}
              glonod = .FALSE.
              if (matstr(word, 'GLOBAL', 3)) glonod = .TRUE.

c ... skip 'node'
              CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
              if (.NOT. MATSTR (WORD, 'NODE', 2)) THEN
                call PRTERR('CMDERR',
     *            'SYNTAX: Move {global|local} NODE {id} {TO|BY} x y z')
                go to 390
              end if

C ... Get the id of the node to be moved
              CALL FFINTG (IFLD, INTYP, IFIELD, 'Node', 0, node, *420)
              if (glonod) then
                node = ng2l('node', node, numnp, mapnod)
              end if
              if (node .le. 0 .or. node .gt. numnp) then
                CALL PRTERR ('CMDERR', 'Node index exceeds node count')
                go to 390
              end if

C ... See if absolute or relative
              CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
              absolute = .false.
              if (MATSTR (WORD, 'TO', 2)) THEN
                absolute = .true.
              else if (MATSTR (WORD, 'BY', 2)) THEN
                absolute = .false.
              else
                call PRTERR('CMDERR',
     *            'Expected "TO" or "BY" in move command')
              go to 390

              end if

              x = 0.0
              y = 0.0
              z = 0.0

              CALL FFREAL (IFLD, INTYP, RFIELD,
     &          'X value', 0.0, X, *420)
              if (ndim .gt. 1) then
                CALL FFREAL (IFLD, INTYP, RFIELD,
     &            'Y value', 0.0, Y, *420)
                if (ndim .gt. 2) then
                  CALL FFREAL (IFLD, INTYP, RFIELD,
     &              'Z value', 0.0, Z, *420)
                end if
              end if

C ... Move the node...
              if (absolute) then
                xn(node) = x
                if (ndim .ge. 2) yn(node) = y
                if (ndim .eq. 3) zn(node) = z
              else
                xn(node) = xn(node) + x
                if (ndim .ge. 2) yn(node) = yn(node) + y
                if (ndim .eq. 3) zn(node) = zn(node) + z
              end if
              VERB = ' '
              go to 420
            ELSE
              call PRTERR('CMDERR',
     *          'Expected "RESET", "GLOBAL", "LOCAL", or sideset id')
              go to 390
            END IF
          END IF

          if (ffexst(ifld, intyp)) then
            numsnp = numsnp + 1
            if (numsnp .gt. mxsnap) then
              call PRTERR('CMDERR', 'Too many snap sets specified')
              go to 390
            end if

            ismtyp(numsnp) = IMOVE
            usnorm(numsnp) = PNORM
            snptol(numsnp) = 0.0
            gap(numsnp)    = 0.0
C ... NOTE: This value is squared in snpnod, don't increase too big.
            delmax(numsnp) = 1.0e15

            call ffintg (ifld, intyp, ifield,
     $        'sideset to move nodes from', 0, IDSSSL(numsnp), *390)
            imat = locint(idsssl(numsnp), numess, idess)
            IF (IMAT .EQ. 0) THEN
              WRITE (STRING, 10050) idsssl(numsnp)
10050         FORMAT ('Sideset ',I5,' does not exist')
              CALL SQZSTR (STRING, LSTR)
              CALL PRTERR ('ERROR', STRING(:LSTR))
              NUMSNP = NUMSNP - 1
              GO TO 390
            END IF

            CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
            IF (.NOT. MATSTR (WORD, 'TO', 1)) THEN
              CALL PRTERR('CMDERR',
     *          'Expected required keyword "TO"')
              NUMSNP = NUMSNP - 1
              GO TO 390
            END IF

            call ffintg (ifld, intyp, ifield,
     $        'sideset to move nodes to', 0, IDSSMA(numsnp), *390)
            imat = locint(idssma(numsnp), numess, idess)
            IF (IMAT .EQ. 0) THEN
              WRITE (STRING, 10050) idssma(numsnp)
              CALL SQZSTR (STRING, LSTR)
              CALL PRTERR ('ERROR', STRING(:LSTR))
              NUMSNP = NUMSNP - 1
              GO TO 390
            END IF
          end if

 380      CONTINUE
          IF (FFEXST (IFLD, INTYP)) THEN
            CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
            IF (MATSTR (WORD, 'RESET', 1)) THEN
              numsnp = 0
              go to 380
            else IF (MATSTR(WORD, 'NORMAL', 1)) THEN
              CALL PRTERR('CMDERR',
     *          '"NORMAL" not supported for "MOVE" command')
              go to 390
            else IF (MATSTR(WORD, 'EDGE', 1)) THEN
              CALL PRTERR('CMDERR',
     *          '"EDGE" not supported for "MOVE" command')
              go to 390
            ELSE IF (WORD .EQ. 'X') THEN
              USNORM(NUMSNP) = PVECT
              VECTOR(1,NUMSNP) = 1.0
              VECTOR(2,NUMSNP) = 0.0
              VECTOR(3,NUMSNP) = 0.0
            ELSE IF (WORD .EQ. 'Y') THEN
              USNORM(NUMSNP) = PVECT
              VECTOR(1,NUMSNP) = 0.0
              VECTOR(2,NUMSNP) = 1.0
              VECTOR(3,NUMSNP) = 0.0
            ELSE IF (WORD .EQ. 'Z') THEN
              USNORM(NUMSNP) = PVECT
              VECTOR(1,NUMSNP) = 0.0
              VECTOR(2,NUMSNP) = 0.0
              VECTOR(3,NUMSNP) = 1.0
            ELSE IF (WORD .EQ. 'MINUSX') THEN
              USNORM(NUMSNP) = PVECT
              VECTOR(1,NUMSNP) = -1.0
              VECTOR(2,NUMSNP) =  0.0
              VECTOR(3,NUMSNP) =  0.0
            ELSE IF (WORD .EQ. 'MINUSY') THEN
              USNORM(NUMSNP) = PVECT
              VECTOR(1,NUMSNP) =  0.0
              VECTOR(2,NUMSNP) = -1.0
              VECTOR(3,NUMSNP) =  0.0
            ELSE IF (WORD .EQ. 'MINUSZ') THEN
              USNORM(NUMSNP) = PVECT
              VECTOR(1,NUMSNP) =  0.0
              VECTOR(2,NUMSNP) =  0.0
              VECTOR(3,NUMSNP) = -1.0
            ELSE IF (MATSTR(WORD, 'VECTOR', 1)) THEN
              USNORM(NUMSNP) = PVECT
              CALL FFREAL (IFLD, INTYP, RFIELD,
     &          'vector X coordinate', 0.0, VECTOR(1,NUMSNP), *390)
              CALL FFREAL (IFLD, INTYP, RFIELD,
     &          'vector Y coordinate', 0.0, VECTOR(2,NUMSNP), *390)
              CALL FFREAL (IFLD, INTYP, RFIELD,
     &          'vector Z coordinate', 0.0, VECTOR(3,NUMSNP), *390)
            ELSE IF (MATSTR(WORD, 'RADIAL', 1)) THEN
              CALL PRTERR('CMDERR',
     *          '"RADIAL" not supported for "MOVE" command')
              go to 390
            ELSE IF (MATSTR(WORD, 'CENTER', 1)) THEN
              CALL PRTERR('CMDERR',
     *          '"CENTER" not supported for "MOVE" command')
              go to 390
            ELSE IF (MATSTR(WORD, 'TOLERANCE', 1)) THEN
              CALL FFREAL (IFLD, INTYP, RFIELD,
     &          'Tolerance', 0.0, SNPTOL(NUMSNP), *390)
            ELSE IF (MATSTR(WORD, 'MAXDELTA', 1)) THEN
              CALL FFREAL (IFLD, INTYP, RFIELD,
     &          'MaxDelta', 1.0E15, DELMAX(NUMSNP), *390)
            ELSE IF (MATSTR(WORD, 'GAP', 1)) THEN
              CALL FFREAL (IFLD, INTYP, RFIELD,
     &          'Gap', 0.0, GAP(NUMSNP), *390)
            END IF
            GO TO 380
          END IF
        end if
 390    CONTINUE
        VERB = 'MOVE'

C --- BEGINNING OF WARP
C========================================================================
      ELSE IF (VERB .EQ. 'WARP') THEN
        IF (NDIM .NE. 3) THEN
          CALL PRTERR('CMDERR',
     *      'Warp cannot be specified for 2D databases')
          GO TO 410
        END IF

        CALL FFCHAR (IFLD, INTYP, CFIELD,
     *    'HELP', WORD)

        IF (MATSTR(WORD, 'HELP', 1)) THEN
          CALL PRTERR ('CMDSPEC',
     *      'Syntax: WARP {XAXIS|YAXIS|ZAXIS|ORIGIN} RADIUS {radius}'//
     *      ' NORMAL {X|Y|Z}')
          VERB = ' '
          GO TO 410
        ELSE IF (MATSTR(WORD, 'XAXIS', 1)) THEN
          IWARP = -1
        ELSE IF (MATSTR(WORD, 'YAXIS', 1)) THEN
          IWARP = -2
        ELSE IF (MATSTR(WORD, 'ZAXIS', 1)) THEN
          IWARP = -3
        ELSE IF (MATSTR(WORD, 'ORIGIN', 1)) THEN
          IWARP = 1

        ELSE
          CALL PRTERR ('CMDERR', 'Invalid WARP Option')
          IWARP = 0
          GO TO 410
        END IF

 400    continue
        IF (FFEXST (IFLD, INTYP)) THEN
          CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
          IF (MATSTR (WORD, 'RADIUS', 1)) THEN
            CALL FFREAL (IFLD, INTYP, RFIELD,
     *        'reference distance', 0.0, WRPDIS, *410)
            IF (WRPDIS .LE. 0.0) THEN
              CALL PRTERR ('CMDERR',
     *          'Reference radius must be greater than zero')
              IWARP = 0
              GO TO 410
            END IF
          ELSE IF (MATSTR (WORD, 'NORMAL', 1)) THEN
            NRMWRP = 0
            CALL FFCHAR (IFLD, INTYP, CFIELD, 'X', WORD)
            IF (MATSTR(WORD, 'X', 1)) THEN
              NRMWRP = 1
            ELSE IF (MATSTR(WORD, 'Y', 1)) THEN
              NRMWRP = 2
            ELSE IF (MATSTR(WORD, 'Z', 1)) THEN
              NRMWRP = 3
            ELSE
              CALL PRTERR ('CMDERR', 'Invalid WARP NORMAL Option')
            END IF
          ELSE
            CALL PRTERR ('CMDERR', 'Invalid WARP Option')
            IWARP = 0
            GO TO 410
          END IF
          GO TO 400
        ELSE
          GO TO 410
        END IF

 410    CONTINUE
        IF (IWARP .NE. 0) THEN
          if (wrpdis .le. 0.0) then
            CALL PRTERR ('CMDERR',
     *        'Reference radius must be greater than zero')
            IWARP = 0
          end if
          if (nrmwrp .eq. 0) then
            CALL PRTERR ('CMDERR',
     *        'Normal Vector Axis must be specified')
            IWARP = 0
          end if
        END IF

C... User Subroutine
C========================================================================
      ELSE IF (VERB .EQ. 'USERSUBROUTINE') then
        call ffonof (ifld, intyp, cfield, USRSUB, *420)

C... Convert nodal variable to element variable (average nodal quantities)
C========================================================================
      ELSE IF (VERB .EQ. 'ELEMENTIZE') then
        CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)

        IF (MATSTR(WORD, 'HELP', 1)) THEN
          CALL PRTERR ('CMDSPEC',
     *      'Syntax: ELEMENTIZE {nodal_variable_name}')
          VERB = ' '
          goto 420
        END IF
C ... Search for name in list of nodal variables.
        index = locstr(word, nvarnp, namenv)

        if (index .eq. 0) then
          write (STRING,*) 'Could not find nodal variable ',
     *      word(:lenstr(word)),
     *      ' Enter "list names" to see valid names.'
          CALL SQZSTR (STRING, LSTR)
          CALL PRTERR ('ERROR', STRING(:LSTR))
        else
          inod2el(index) = 1
          write (STRING,*) 'Elementizing nodal variable ',
     *      word(:lenstr(word))
          CALL SQZSTR (STRING, LSTR)
          CALL PRTERR ('CMDSPEC', STRING(:LSTR))
        end if
        VERB = ' '

C... Output element centroids
C========================================================================
      ELSE IF (VERB .EQ. 'CENTROIDS') then
        call ffonof (ifld, intyp, cfield, CENTRD, *420)

C...  Help
C========================================================================
      ELSE IF (VERB .EQ. 'HELP') THEN
        ISHELP = HELP (MYNAME(:LNAM), 'COMMANDS', CFIELD(IFLD))
        IF (.NOT. ISHELP) THEN
          CALL PRTERR('CMDREQ', 'Help is not available at this time')
        END IF
        VERB = ' '

C========================================================================
      ELSE IF ((VERB .EQ. 'END') .OR. (VERB .EQ. 'EXIT')) THEN
        CALL SCNEOF
        GOTO 430

C========================================================================
      ELSE IF (VERB .EQ. 'EXECUTE') THEN
        CALL SHOW(VERB, ' ')
        EXECUT = .TRUE.
        GOTO 430

C========================================================================
      ELSE IF (VERB .EQ. 'QUIT') THEN
        CALL SCNEOF
        RETURN 1

C========================================================================
      ELSE
        CALL PRTERR ('CMDERR', '"' // VERB(:LENSTR(VERB))
     &    // '" is an invalid command')
        VERB = ' '
      END IF

 420  CONTINUE

      IF (VERB .NE. ' ') THEN
        CALL SHOW (VERB, ' ')
      END IF

      GOTO 110

 430  CONTINUE

      RETURN
      END

C=======================================================================
      integer function ng2l (type, idglo, num, map)
C=======================================================================

C   --G2L locates the local node/element corresponding to the global id
C   --
C   --Parameters:
C   --   TYPE - IN - 'node' or 'element'
C   --   NUM  - IN - the number of node/elements
C   --   MAP  - IN - the node/element number map

C ... TYPE = Node or Element
      CHARACTER*(*) TYPE
      INTEGER MAP(*)

      CHARACTER*128 STRA
      LOGICAL FOUND

      ng2l = idglo
      found = .FALSE.
      do i = 1, num
        if (idglo .eq. map(i)) then
          ng2l = i
          return
        end if
      end do

      write (stra, 10000) type, idglo
      call sqzstr(stra, lstra)
      call prterr('WARNING', stra(:lstra))
      RETURN

10000  FORMAT (' No local ',A,' has global id equal to ',I10)
      END
