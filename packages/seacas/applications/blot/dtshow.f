C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DTSHOW (SHOTYP, NAMES, LIDSP)
C=======================================================================

C   --*** DTSHOW *** (DETOUR) Display DETOUR parameter information
C   --   Written by Amy Gilkey - revised 04/08/88
C   --
C   --DTSHOW displays the DETOUR plot parameters.
C   --
C   --The SHOW options with the items they display are:
C   --   PLOT     - display mode and plot variables
C   --   HARDCOPY -
C   --   VIEW     -
C   --   WIREFRAM -
C   --   SOLID    -
C   --   CONTOUR  -
C   --   PAINT    -
C   --   SYMBOL   -
C   --   VECTOR   -
C   --   SIGMAX   -
C   --   SIGMIN   -
C   --   GAUSS    -
C   --   NCNTRS   - minimum and maximum contour interval values, contour
C   --   CRANGE   -    interval and number of intervals
C   --   CMIN     -
C   --   CMAX     -
C   --   CSHIFT   -
C   --   DELCNTR  -
C   --   CINTV    - contour intervals (plus NCNTRS info)
C   --   CLABEL   - number of interior contour letters
C   --   COPEN    - whether painted contour limits are "open"
C   --   CSYMBOLS - number of min/max symbols to disable
C   --   VSCALE   - vector/symbol scale factor
C   --
C   --Parameters:
C   --   SHOTYP - IN - the expanded SHOW option string
C   --   NAMES - IN - the variable names
C   --   LIDSP(0:*)  - IN/OUT - the indices of the selected variables
C   --          whose values will be displayed on the plot legend.
C   --          ABS(LIDSP(0)) = the number of variables in the list.
C   --          SIGN(LIDSP(0)) specifies whether the variables in the
C   --                   list should have their values displayed on
C   --                   the plot legend.  If >0, they should;
C   --                   If <=0, they should not.
C   --          LIDSP(i) identifies the ith variable in the list.
C   --          If LIDSP(i) > 0, LIDSP(i) is the id of a history variable.
C   --          If LIDSP(i) < 0, -LIDSP(i) is the id of a global variable.
C   --          If LIDSP(i) = 0, TIME is to be displayed on the plot legend.
C   --
C   --Common Variables:
C   --   Uses NDIM, NELBLK of /DBNUMS/
C   --   Uses IS3DIM of /D3NUMS/
C   --   Uses VECSCL of /ETCOPT/
C   --   Uses MSHDEF, MSHLIN of /MSHOPT/
C   --   Uses MODDET, MODTYP, IDTVAR, NNDVAR, NEDVAR of /DETOPT/
C   --   Uses CINTOK, NCNTR, NOCMIN, NOCMAX, LABINC, MAXMIN, MAXMAX of /CNTR/

      include 'params.blk'
      PARAMETER (MSHNON=0, MSHBOR=1, MSHDIV=2, MSHSEL=3, MSHALL=4)

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)

      include 'dbnums.blk'
      include 'dbnumgq.blk'
      include 'd3nums.blk'
      include 'etcopt.blk'
      include 'mshopt.blk'
      include 'detopt.blk'
      include 'cntr.blk'

      CHARACTER*(*) SHOTYP
      CHARACTER*(*) NAMES(*)
      INTEGER LIDSP(0:*)

      LOGICAL ISABRT
      INTEGER NUMMOD
      INTEGER MAINVW, NDEFVW, IXVW

      CHARACTER*1024 STRING
      CHARACTER CH
      CHARACTER*(MXNAME) NAM(4)
      CHARACTER*20 RSTR(9)
      REAL RNUM(9)
      CHARACTER*(MXNAME) TNAME(5)

C *** Display mode control and multiple views mode control ***

      IF ((SHOTYP .EQ. 'PLOT') .OR. (SHOTYP .EQ. 'HARDCOPY')
     &   .OR. (SHOTYP .EQ. 'VIEW')
     &   .OR. (SHOTYP .EQ. 'WIREFRAM') .OR. (SHOTYP .EQ. 'SOLID')
     &   .OR. (SHOTYP .EQ. 'NSETS') .OR. (SHOTYP .EQ. 'SSETS')
     &   .OR. (SHOTYP .EQ. 'CONTOUR') .OR. (SHOTYP .EQ. 'PAINT')
     &   .OR. (SHOTYP .EQ. 'EPAINT')
     &   .OR. (SHOTYP .EQ. 'SYMBOL') .OR. (SHOTYP .EQ. 'VECTOR')
     &   .OR. (SHOTYP .EQ. 'SIGMAX') .OR. (SHOTYP .EQ. 'SIGMIN')
     &   .OR. (SHOTYP .EQ. 'GAUSS')) THEN

C      --Get the variable names
         DO 100 I = 1, MAX (NNDVAR, NEDVAR)
            IF (IDTVAR(I) .GT. 0) THEN
               NAM(I) = NAMES(IDTVAR(I))
            ELSE IF ((NUMMOD (MODDET, ' ', 'VECTOR', ' ') .GE. 1)) THEN
               NAM(I) = '0'
            ELSE
               NAM(I) = '---'
            END IF
  100    CONTINUE

C      --Find the "main" view on symmetric views
         MAIN = MAINVW ()

         NVIEWS = NDEFVW (.TRUE.)
         DO 110 IVW = 1, NVIEWS
            IVIEW = IXVW (.TRUE., IVW)
            IF (MSHDEF(IVIEW) .EQ. 'EMPTY') THEN
               STRING = 'Empty view'
            ELSE IF (MODDET(IVIEW) .EQ. 'WIREFRAM') THEN
               WRITE (STRING, 10000) 'Wireframe mesh'
10000           FORMAT (A, ' plot', :, ' of ', 4 (A, :, ' '))
            ELSE IF (MODDET(IVIEW) .EQ. 'SOLID') THEN
               WRITE (STRING, 10000) 'Solid mesh'
            ELSE IF (MODDET(IVIEW) .EQ. 'CONTOUR') THEN
               IF (MODTYP(IVIEW) .EQ. 'LINE') THEN
                  WRITE (STRING, 10000) 'Line contour', NAM(1)
               ELSE IF (MODTYP(IVIEW) .EQ. 'PAINT') THEN
                  WRITE (STRING, 10000) 'Paint contour', NAM(1)
               END IF
            ELSE IF (MODDET(IVIEW) .EQ. 'ELEMCONT') THEN
               IF (MODTYP(IVIEW) .EQ. 'PAINT') THEN
                  WRITE (STRING, 10000) 'Paint element contour', NAM(1)
               END IF
            ELSE IF (MODDET(IVIEW) .EQ. 'VECTOR') THEN
               IF ((MODTYP(IVIEW) .EQ. 'NODE')
     &            .OR. (MODTYP(IVIEW) .EQ. 'ELEMENT')) THEN
                  WRITE (STRING, 10000) 'Vector', (NAM(I), I=1,NDIM)
               ELSE IF (MODTYP(IVIEW) .EQ. 'SIGMAX') THEN
                  WRITE (STRING, 10000)
     &               'Maximum principal stress vector', (NAM(I), I=1,3)
               ELSE IF (MODTYP(IVIEW) .EQ. 'SIGMIN') THEN
                  WRITE (STRING, 10000)
     &               'Minimum principal stress vector', (NAM(I), I=1,3)
               END IF
            ELSE IF (MODDET(IVIEW) .EQ. 'SYMBOL') THEN
               IF (MODTYP(IVIEW) .EQ. 'ANGLE') THEN
                  WRITE (STRING, 10000) 'Angle symbol', NAM(1)
               ELSE IF (MODTYP(IVIEW) .EQ. 'CRACK') THEN
                  WRITE (STRING, 10000) 'Crack symbol', NAM(1)
               ELSE IF (MODTYP(IVIEW) .EQ. 'STATE') THEN
                  WRITE (STRING, 10000) 'State symbol', NAM(1)
               ELSE IF (MODTYP(IVIEW) .EQ. 'SPHERE') THEN
                  STRING = 'Spherical element plot'
               ELSE IF (MODTYP(IVIEW) .EQ. 'FSPHERE') THEN
                  STRING = 'Filled spherical element plot'
               ELSE
                  WRITE (STRING, 10000) 'Symbol', NAM(1)
               END IF
            ELSE IF (MODDET(IVIEW) .EQ. 'GAUSS') THEN
               WRITE (STRING, 10000) 'Gauss symbol', (NAM(I), I=1,4)
            ELSE
               WRITE (STRING, 10000) 'UNKNOWN'
            END IF
            CALL SQZSTR (STRING, LSTR)

            IF (NVIEWS .GT. 1) THEN
               WRITE (CH, '(I1)') IVIEW
               IF (IVIEW .EQ. MAIN) THEN
                  WRITE (*, 10080) 'VIEW ', CH, ':* ', STRING(:LSTR)
               ELSE
                  WRITE (*, 10080) 'VIEW ', CH, ':  ', STRING(:LSTR)
               END IF
            ELSE
               WRITE (*, 10080) STRING(:LSTR)
            END IF
  110    CONTINUE

C *** Contour control ***

      ELSE IF ((SHOTYP .EQ. 'NCNTRS') .OR. (SHOTYP .EQ. 'CRANGE')
     &   .OR. (SHOTYP .EQ. 'CMIN') .OR. (SHOTYP .EQ. 'CMAX')
     &   .OR. (SHOTYP .EQ. 'CSHIFT') .OR. (SHOTYP .EQ. 'DELCNTR')
     &   .OR. (SHOTYP .EQ. 'CINTV')) THEN

         IF ((.NOT. CINTOK) .OR. (SHOTYP .NE. 'CINTV')) THEN
            IF (.NOT. CINTOK) THEN
               RNUM(1) = CMIN
               RNUM(2) = CMAX
               RNUM(3) = DELC
               CALL NUMSTR (3, 4, RNUM, RSTR, LSTR)
               WRITE (STRING, 10010, IOSTAT=IDUM)
     &            (RSTR(I)(:LSTR), I=1,3), NCNTR
10010           FORMAT ('Contour: ', A, ' to ', A,
     &            ' by ', A, ',', I6, ' contour')
            ELSE
               RNUM(1) = CMIN
               RNUM(2) = CMAX
               CALL NUMSTR (2, 4, RNUM, RSTR, LSTR)
               WRITE (STRING, 10020, IOSTAT=IDUM)
     &            (RSTR(I)(:LSTR), I=1,2), NCNTR
10020           FORMAT ('Defined contours: ', A, ' to ', A,
     &            ',', I6, ' contour')
            END IF
            CALL SQZSTR (STRING, LSTR)
            IF (LINCON) THEN
               WRITE (*, 10080) STRING(:LSTR), ' lines'
            ELSE
               WRITE (*, 10080) STRING(:LSTR), ' areas'
            END IF
         END IF

         IF (SHOTYP .EQ. 'CINTV') THEN
            WRITE (*, 10080)
            IF (.NOT. CINTOK) THEN
               WRITE (*, 10080) 'Contours:'
            ELSE
               WRITE (*, 10080) 'Defined Contours:'
            END IF
            RNUM(1) = CMIN
            RNUM(2) = CMAX
            IF (LINCON) THEN
               NC = NCNTR
            ELSE
               NC = NCNTR+1
            END IF
            DO 120 I = 1, NC
               IF (ISABRT ()) RETURN
               RNUM(3) = CNTRI (I)
               CALL NUMSTR (3, 4, RNUM, RSTR, LSTR)
               WRITE (*, 10030, IOSTAT=IDUM) I, RSTR(3)(:LSTR)
10030           FORMAT (1X, I5, ') ', A)
  120       CONTINUE
         END IF

C *** Display options ***

      ELSE IF ((SHOTYP .EQ. 'CLABEL') .or. (shotyp .eq. 'LINES')) THEN
         IF (LABINC .LT. 0) THEN
            WRITE (*, 10080) 'No contour labels'
         ELSE IF (LABINC .EQ. 0) THEN
            WRITE (*, 10080) 'Contours labeled',
     &         ' on mesh and element block boundaries only'
         ELSE
            CALL INTSTR (1, 0, LABINC, STRING, LSTR)
            WRITE (*, 10080) 'Contour label for every ', STRING(:LSTR),
     &         ' interior mesh lines'
         END IF

      ELSE IF ((SHOTYP .EQ. 'COPEN') .or. (shotyp .eq. 'OPENCNTR')) THEN
         IF (NOCMIN) THEN
            WRITE (STRING, 10040) 'minimum'
10040        FORMAT ('Contour ', A, ' open (painted)')
         ELSE
            WRITE (STRING, 10050) 'minimum'
10050        FORMAT ('Contour ', A, ' closed (NOT painted)')
         END IF
         CALL SQZSTR (STRING, LSTR)
         WRITE (*, 10080) STRING(:LSTR)
         IF (NOCMAX) THEN
            WRITE (STRING, 10040) 'maximum'
         ELSE
            WRITE (STRING, 10050) 'maximum'
         END IF
         CALL SQZSTR (STRING, LSTR)
         WRITE (*, 10080) STRING(:LSTR)

      ELSE IF ((SHOTYP .EQ. 'CSYMBOLS')
     &   .or. (shotyp .eq. 'MINMAX')) THEN
         WRITE (STRING, 10060, IOSTAT=IDUM)
     &      'minimum = ', MAXMIN, ', maximum =', MAXMAX
10060     FORMAT (A, I5, A, I5)
         CALL SQZSTR (STRING, LSTR)
         WRITE (*, 10080)
     &      'Maximum number of symbols to display: ', STRING(:LSTR)

      ELSE IF ((SHOTYP .EQ. 'VSCALE') .or. (shotyp .eq. 'VECSCL')) THEN
         CALL NUMSTR1(4, VECSCL, RSTR, LSTR)
         WRITE (*, 10080) 'Vector/symbol scale factor = ',
     &      RSTR(1)(:LSTR)

      ELSE IF (SHOTYP .EQ. 'DISPVAR') THEN
         IF (LIDSP(0) .GT. 0) THEN
            WRITE (*, 10080) 'Values of display variables will',
     &         ' appear on plot legend'
         ELSE
            WRITE (*, 10080) 'Values of display variables will not',
     &         ' appear on plot legend'
         ENDIF
         WRITE (*, 10080)
         WRITE (*, 10080) 'Display variables - '
         NUMDSP = ABS(LIDSP(0))
         NUM = 0
         NLAST = 0
  130    CONTINUE
         IF (NUM .LT. NUMDSP) THEN
            NFIRST = NLAST + 1
            NLAST = NFIRST + 4
            IF (NLAST .GT. NUMDSP) NLAST = NUMDSP
            N = NLAST - NFIRST + 1
            NID = NFIRST
            DO 140 I = 1, N
               IF (LIDSP(NID) .EQ. 0) THEN
                  TNAME(I) = 'TIME'
               ELSE IF (LIDSP(NID) .GT. 0) THEN
                  CALL DBVIX_BL ('H', LIDSP(NID), IDVAR)
                  TNAME(I) = NAMES(IDVAR)
               ELSE IF (LIDSP(NID) .LT. 0) THEN
                  CALL DBVIX_BL ('G', -LIDSP(NID), IDVAR)
                  TNAME(I) = NAMES(IDVAR)
               ENDIF
               NID = NID + 1
  140       CONTINUE
            WRITE (*, 10070)(TNAME(I), I=1,N)
10070        FORMAT (5X, 5(A32,2X))
            NUM = NUM + N
            GO TO 130
         ENDIF

      END IF

      RETURN
10080  FORMAT (1X, 10A)
      END
