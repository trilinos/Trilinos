C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE COMAND (IDNPS, IDESS, IDNSET, IDESET,
     &     BLKTYP, IBPARM, IDELB, NUMELB, NUMLNK, NAMELB, ELATTR,
     &     XN, YN, A, *)
C=======================================================================

C   --*** COMAND *** (GEN3D) Input and process commands
C   --   Written by Amy Gilkey - revised 05/09/88
C   --   Modified by Greg Sjaardema - 02/06/89
C   --       Added Warp Function
C   --       Added Gradient to Rotations (not for center blocks)
C   --  5/11 Added GETINT call to get intervals/distance/gradient
C   --       Added INTERVALS command to set  ""    ""       ""
C   --
C   --COMAND inputs and executes an user command.
C   --
C   --Parameters:
C   --   IDNPS - IN - the IDs of existing node sets
C   --   IDESS - IN - the IDs of existing side sets
C   --   IDNSET - OUT - the IDs of the front and back surface node sets;
C   --      (0) = length
C   --   IDESET - OUT - the IDs of the front and back surface side sets;
C   --      (0) = length
C   --   BLKTYP - IN - the element block type
C   --   IBPARM - IN - the block parameters (defined by the block type)
C   --   IDELB - IN - the ids for each block
C   --   NUMELB - IN - the number of elements for each block
C   --   NUMLNK - IN - the number of nodes per element for each block
C   --   XN, YN - IN - the 2D nodal coordinates
C   --   * - return statement for QUIT
C   --
C   --Common Variables:
C   --   Sets NDIM, NUMNP, NUMEL, NELBLK,
C   --      NUMNPS, LNPSNL, NUMESS, LESSEL, LESSNL of /DBNUMS/
C   --   Sets ITRANT, NNREPL, NEREPL, DIM3, NRTRAN, D3TRAN, ZGRAD,
C   --      CENTER, NUMCOL of /PARAMS/
C   --   Sets XOFFS, YOFFS, ZOFFS of /XYZOFF/
C   --   Sets ROT3D, ROTMAT, ROTCEN of /XYZROT/

      PARAMETER (MAXFLD = 1000)

      INCLUDE 'g3_dbnums.blk'
      INCLUDE 'g3_params.blk'
      INCLUDE 'g3_xyzoff.blk'
      INCLUDE 'g3_xyzrot.blk'
      INCLUDE 'g3_xyzmir.blk'
      INCLUDE 'g3_xyzero.blk'
      INCLUDE 'g3_xyzscl.blk'
      INCLUDE 'g3_twist.blk'
      INCLUDE 'g3_cmdsho.blk'
      INCLUDE 'g3_xxxxx.blk'
      INCLUDE 'g3_splxyz.blk'

      INTEGER IDNPS(*), IDESS(*)
      INTEGER IDNSET(0:MAXSET,2), IDESET(0:MAXSET,2)
      CHARACTER BLKTYP(*)
      CHARACTER*32 NAMELB(*)
      REAL ELATTR(*)
      INTEGER IBPARM(4,*)
      INTEGER IDELB(*)
      INTEGER NUMELB(*)
      INTEGER NUMLNK(*)
      REAL XN(*), YN(*)
      REAL A(*)

      LOGICAL FFEXST, MATSTR, HELP, FFNUMB

      CHARACTER*8 WORD, VERB, MYNAME
      INTEGER     INTYP(MAXFLD+1)
      CHARACTER*8 CFIELD(MAXFLD)
      INTEGER     IFIELD(MAXFLD)
      REAL        RFIELD(MAXFLD)
      LOGICAL ISFRON
      LOGICAL ANYCEN
      LOGICAL ISHELP

      CHARACTER*8 CMDTBL(36)
      SAVE CMDTBL
C      --CMDTBL - the valid commands table

C   --Command table follows.  Remember to change the dimensioned size when
C   --changing the table.
      DATA CMDTBL /
     1   'TRANSLAT', 'ROTATE  ', 'OFFSET  ', 'REVOLVE ', 'REVCEN  ',
     2   'WARP    ', 'BLOCK   ', 'CENTER  ', 'TUNNEL  ', 'SPECIAL ',
     3   'NSETS   ', 'NODESETS', 'SSETS   ', 'SIDESETS', 'MIRROR  ',
     4   'LIST    ', 'SHOW    ', 'HELP    ', 'ZERO    ', 'SCALE   ',
     5   'END     ', 'EXIT    ', 'QUIT    ', 'TWIST   ', 'PROJECT ',
     6   'SUMMARY ', 'INTERVAL', 'ROTCEN  ', 'ROTAXIS ', 'SPLINE  ',
     7   'CHANGE  ', 'TRANSPLI', 'SHIFT   ', 'LIMITS  ', 'ATTRIBUT',
     &   '        ' /

      DATA MYNAME /'GEN3D   '/

C   --Initialize

      LNAM = LENSTR(MYNAME)
      ITRANT = 0
      NEREPL = 1
      NNREPL = NEREPL + 1
      DIM3 = 1.0
      NRTRAN(1) = NEREPL
      D3TRAN(1) = DIM3
      ZGRAD(1) = 1.0
      RGRAD = 1.0
      NDEGR = 0
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
      TWYCEN = 0.0
      TWXCEN = 0.0
      ROT3D = .FALSE.
      ISCENT = .FALSE.
      ROTAX = 0
      CALL INIREA (3*3, 0.0, ROTMAT)
      DO 10 I = 1, 3
         ROTMAT(I,I) = 1.0
   10 CONTINUE

C     -- Zero the elattr array here.  If set by user, then copy these
C     -- attributes into the attrib array in wrelb.  Else, use the
C     -- attributes set in the input genesis file.
      call inirea (nelblk, 0.0, elattr)

      CALL MINMAX (NUMNP, XN, XMIN, XMAX)
      ROTCEN(1) = XMIN
      CALL MINMAX (NUMNP, YN, YMIN, YMAX)
      ROTCEN(2) = YMIN
      ROTCEN(3) = 0.0
      CPOINT = .FALSE.

      IDNSET(0,1) = 0
      IDNSET(0,2) = 0
      IDESET(0,1) = 0
      IDESET(0,2) = 0

      CALL INISTR (NELBLK, ' ', BLKTYP)
      CALL INIINT (4 * NELBLK, 0, IBPARM)
      NUMCOL = 0
      NUMCDM = MAX(1, NUMCOL)
      NUMROW = 0
      ANYCEN = .FALSE.
   20 CONTINUE

C   --Read command line

      WRITE (*, *)
      CALL FREFLD (0, 0, 'GEN3D> ', MAXFLD,
     &   IOSTAT, NUMFLD, INTYP, CFIELD, IFIELD, RFIELD)
      IF (IOSTAT .LT. 0) GOTO 220
      IF (NUMFLD .EQ. 0) GOTO 20
      INTYP(MIN(MAXFLD,NUMFLD)+1) = -999

      IFLD = 1
      CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
      CALL ABRSTR (VERB, WORD, CMDTBL)
      IF (VERB .EQ. ' ') VERB = WORD

C   --Perform command

      IF (VERB .EQ. 'TRANSLAT') THEN
         ITRANT = 1
         CALL GETINT ('translation', IFLD, INTYP, IFIELD, RFIELD,
     *      NBLK, NRTRAN, D3TRAN, ZGRAD, NEREPL, NNREPL, DIM3,
     *      MAXINT, *170)
         NUMCOL = -999
         NUMCDM = MAX(1, NUMCOL)

      ELSE IF (VERB .EQ. 'TRANSPLI') THEN
         ITRANT = 128
         CALL GETINT ('translation', IFLD, INTYP, IFIELD, RFIELD,
     *      NBLK, NRTRAN, D3TRAN, ZGRAD, NEREPL, NNREPL, DIM3,
     *      MAXINT, *170)
         NUMCOL = -999
         NUMCDM = MAX(1, NUMCOL)

         CALL GETSPT (A)

      ELSE IF (VERB .EQ. 'SPLINE') THEN
         ITRANT = 64
         CALL GETINT ('thickness', IFLD, INTYP, IFIELD, RFIELD,
     *      NBLK, NRTRAN, D3TRAN, ZGRAD, NEREPL, NNREPL, DIM3,
     *      MAXINT, *170)
         NUMCOL = -999
         NUMCDM = MAX(1, NUMCOL)

         CALL GETSPL (A)

C --- Beginning of Intervals Command
      ELSE IF (VERB .EQ. 'INTERVAL') THEN
         CALL GETINT ('translation', IFLD, INTYP, IFIELD, RFIELD,
     *      NBLK, NRTRAN, D3TRAN, ZGRAD, NEREPL, NNREPL, DIM3,
     *      MAXINT, *170)

         IF (ITRANT .EQ. 2) THEN
           NDEGR = NINT (DIM3)
           IF (NDEGR .LT. 360) THEN
             NNREPL = NEREPL + 1
           ELSE
             NNREPL = NEREPL
           END IF

           IF (NDEGR .EQ. 360) THEN
             IF (IDNSET(0,2) .GT. 0) THEN
               CALL PRTERR ('CMDWARN',
     &           'Back node sets are deleted'
     &           // ' for 360-degree rotation')
             END IF
             IF (IDESET(0,2) .GT. 0) THEN
               CALL PRTERR ('CMDERR',
     &           'Back side sets are deleted'
     &           // ' for 360-degree rotation')
             END IF
           END IF
         END IF
C --- BEGINNING OF WARP
       ELSE IF (VERB .EQ. 'WARP') THEN
         ITRANT = 4
         CALL FFCHAR (IFLD, INTYP, CFIELD,
     *      'POINT', WORD)
         IF((.NOT. MATSTR(WORD, 'LINE',    1)) .AND.
     *      (.NOT. MATSTR(WORD, 'XAXIS',   1)) .AND.
     *      (.NOT. MATSTR(WORD, 'YAXIS',   1)) .AND.
     *      (.NOT. MATSTR(WORD, 'ELLIPSE', 1)) .AND.
     *      (.NOT. MATSTR(WORD, 'POINT',   1))) THEN
            CALL PRTERR ('CMDERR', 'Invalid WARP Option')
            GOTO 170
         END IF

         CALL GETINT ('translation', IFLD, INTYP, IFIELD, RFIELD,
     *      NBLK, NRTRAN, D3TRAN, ZGRAD, NEREPL, NNREPL, DIM3, 1,
     *      *170)

         IF (MATSTR(WORD, 'LINE', 1)) THEN
            IWARP = 2
C            CALL FFREAL (IFLD, INTYP, RFIELD,
C     *         'line X1 coordinate', 0.0, CWARP(1), *200)
C            CALL FFREAL (IFLD, INTYP, RFIELD,
C     *         'line Y1 coordinate', 0.0, CWARP(2), *200)
C            CALL FFREAL (IFLD, INTYP, RFIELD,
C     *         'line Z1 coordinate', 0.0, CWARP(3), *200)
C            CALL FFREAL (IFLD, INTYP, RFIELD,
C     *         'line X2 coordinate', 0.0, CWARP(4), *200)
C            CALL FFREAL (IFLD, INTYP, RFIELD,
C     *         'line Y2 coordinate', 0.0, CWARP(5), *200)
C            CALL FFREAL (IFLD, INTYP, RFIELD,
C     *         'line Z2 coordinate', 0.0, CWARP(6), *200)

         ELSE IF (MATSTR(WORD, 'XAXIS', 1)) THEN
            IWARP = -1

         ELSE IF (MATSTR(WORD, 'YAXIS', 1)) THEN
            IWARP = -2

         ELSE IF (MATSTR(WORD, 'POINT', 1)) THEN
            IWARP = 1
C            CALL FFREAL (IFLD, INTYP, RFIELD,
C     *         'point X coordinate', 0.0, CWARP(1), *110)
C            CALL FFREAL (IFLD, INTYP, RFIELD,
C     *         'point Y coordinate', 0.0, CWARP(2), *110)
C            CALL FFREAL (IFLD, INTYP, RFIELD,
C     *         'point Z coordinate', 0.0, CWARP(3), *110)
         ELSE IF (MATSTR(WORD, 'ELLIPSE', 1)) THEN
            IWARP = 2
            CALL FFREAL (IFLD, INTYP, RFIELD,
     *         'ellipse major axis', 0.0, HRAD, *170)
         ELSE
            IWARP = 0
            GO TO 170
         END IF
         CALL FFREAL (IFLD, INTYP, RFIELD,
     *      'warping distance', 0.0, DWARP, *170)

         CALL FFCHAR (IFLD, INTYP, CFIELD,
     *      'RADIAL', WORD)

         IF (MATSTR(WORD, 'VERTICAL', 1)) THEN
            VEDGE = .TRUE.
         ELSE IF (MATSTR(WORD, 'RADIAL', 1)) THEN
            VEDGE = .FALSE.
         ELSE
            CALL PRTERR ('CMDERR', 'Invalid edge specification')
            GOTO 170
         END IF
         NUMCOL = -999
         NUMCDM = MAX(1, NUMCOL)

C --- END OF WARP

      ELSE IF (VERB .EQ. 'ROTATE') THEN
         ITRANT = 2
         CALL GETINT ('rotation', IFLD, INTYP, IFIELD, RFIELD,
     *     NBLK, NRTRAN(1), D3TRAN(1), ZGRAD(1), NEREPL, NNREPL,
     *     DIM3, 1, *170)

         CALL FFREAL (IFLD, INTYP, RFIELD,
     &      'center of rotation', 0.0, CENTER, *170)

         NEREPL = NRTRAN(1)
         NDEGR = NINT (D3TRAN(1))
         IF (NDEGR .LT. 360) THEN
            NNREPL = NEREPL + 1
         ELSE
            NNREPL = NEREPL
         END IF
         DIM3 = D3TRAN(1)

         IF (NDEGR .EQ. 360) THEN
            IF (IDNSET(0,2) .GT. 0) THEN
               CALL PRTERR ('CMDWARN',
     &            'Back node sets are deleted'
     &            // ' for 360-degree rotation')
            END IF
            IF (IDESET(0,2) .GT. 0) THEN
               CALL PRTERR ('CMDERR',
     &            'Back side sets are deleted'
     &            // ' for 360-degree rotation')
            END IF
         END IF

         NUMCOL = -999
         NUMCDM = MAX(1, NUMCOL)

         IF (ANYCEN .AND. (.NOT. CPOINT)) THEN
            NUMCOL = -1
            NUMCDM = MAX(1, NUMCOL)

            IF (MOD (NDEGR, 90) .NE. 0) THEN
               CALL PRTERR ('CMDWARN', 'Rotation around the mesh edge'
     &            // ' must be 90, 180, 270 or 360 degrees')
               GOTO 30
            END IF

            N45 = NINT (DIM3 / 45.0)
            IF (MOD (NEREPL, N45) .NE. 0) THEN
               CALL PRTERR ('CMDERR', 'Number of rotations'
     &            // ' around the mesh edge is incorrect')
               WRITE (*, 250) '(multiple of 2 for 90 degrees,'
     &            , ' 4 for 180, 6 for 270, 8 for 360)'
               GOTO 30
            END IF

            NUMCOL = -999
            NUMCDM = MAX(1, NUMCOL)
   30       CONTINUE
         END IF

       ELSE IF (VERB .EQ. 'ROTAXIS ') THEN
         CALL FFCHAR (IFLD, INTYP, CFIELD, 'Y', WORD)
         IF (MATSTR(WORD, 'X', 1)) THEN
           ROTAX = 1
         ELSE IF (MATSTR(WORD, 'Y', 1)) THEN
           ROTAX = 0
         ELSE
           CALL PRTERR('CMDERR',
     *       'Invalid ROTAXIS Option, Expected "X" or "Y"')
           GOTO 170
         END IF

C --- Experimental Rotate Routine -- Sets Center, Intervals set by
C        Intervals command
      ELSE IF (VERB .EQ. 'ROTCEN  ') THEN
         ITRANT = 32

         CALL FFREAL (IFLD, INTYP, RFIELD,
     &      'center of rotation', 0.0, CENTER, *170)

         NEREPL = 0
         DEGR = 0.
         DO 40 IBLK = 1, NBLK
            NEREPL = NEREPL + NRTRAN(IBLK)
            DEGR = DEGR + D3TRAN(IBLK)
   40    CONTINUE

C     ... Give tolerance of 1.0e-5
         deg_min = 360.0 - 1.0e-5
         deg_max = 360.0 + 1.0e-5

         if (degr .gt. deg_min .and. degr .lt. deg_max) then
            degr = 360.0
         end if
         NDEGR = NINT (degr)

         IF (DEGR .GT. 360.0) THEN
            CALL PRTERR ('CMDERR',
     &         'Total number of degrees exceeds 360')
            GOTO 170
         END IF

         IF (NEREPL .EQ. 0) THEN
            CALL PRTERR ('CMDERR',
     $           'INTERVAL command must be specified before ROTCEN')
            GO TO 170
         END IF

         IDEGR = NINT (DEGR/NEREPL)
         IF (IDEGR .GE. 180) THEN
            CALL PRTERR ('CMDERR',
     &         'Single rotation cannot cover 180 degrees')
            GOTO 170
         END IF
         IF (DEGR .LT. 360.0) THEN
            NNREPL = NEREPL + 1
         ELSE
            NNREPL = NEREPL
         END IF

         IF (DEGR .EQ. 360.0) THEN
            IF (IDNSET(0,2) .GT. 0) THEN
               CALL PRTERR ('CMDWARN',
     &            'Back node sets are deleted'
     &            // ' for 360-degree rotation')
            END IF
            IF (IDESET(0,2) .GT. 0) THEN
               CALL PRTERR ('CMDERR',
     &            'Back side sets are deleted'
     &            // ' for 360-degree rotation')
            END IF
         END IF

         NUMCOL = -999
         NUMCDM = MAX(1, NUMCOL)

         IF (ANYCEN .AND. (.NOT. CPOINT)) THEN
            NUMCOL = -1
            NUMCDM = MAX(1, NUMCOL)

            IF (MOD (NDEGR, 90) .NE. 0) THEN
               CALL PRTERR ('CMDWARN', 'Rotation around the mesh edge'
     &            // ' must be 90, 180, 270 or 360 degrees')
               GOTO 50
            END IF

            N45 = NINT (DIM3 / 45.0)
            IF (MOD (NEREPL, N45) .NE. 0) THEN
               CALL PRTERR ('CMDERR', 'Number of rotations'
     &            // ' around the mesh edge is incorrect')
               WRITE (*, 250) '(multiple of 2 for 90 degrees,'
     &            , ' 4 for 180, 6 for 270, 8 for 360)'
               GOTO 50
            END IF

            NUMCOL = -999
            NUMCDM = MAX(1, NUMCOL)
   50       CONTINUE
         END IF

      ELSE IF (VERB .EQ. 'TWIST') THEN
         ITRANT = 8
         CALL FFREAL (IFLD, INTYP, RFIELD,
     *      'twist angle', TWANGL, TWANGL, *170)
         IF (.NOT. FFNUMB (IFLD, INTYP)) GO TO 60
         CALL FFREAL (IFLD, INTYP, RFIELD,
     *      'twist X center', TWXCEN, TWXCEN, *170)
         IF (.NOT. FFNUMB (IFLD, INTYP)) GO TO 60
         CALL FFREAL (IFLD, INTYP, RFIELD,
     *      'twist Y center', TWYCEN, TWYCEN, *170)

   60    CONTINUE

         CALL FFCHAR (IFLD, INTYP, CFIELD,
     *      'TRANSLAT', WORD)
         IF (MATSTR(WORD, 'TRANSLAT', 1)) THEN
            ITWTYP = 1
            CALL GETINT ('translation', IFLD, INTYP, IFIELD, RFIELD,
     *         NBLK, NRTRAN, D3TRAN, ZGRAD, NEREPL, NNREPL, DIM3,
     *         MAXINT, *170)
            NUMCOL = -999
            NUMCDM = MAX(1, NUMCOL)

         ELSE IF (MATSTR(WORD, 'ROTATE', 1)) THEN
            ITWTYP = 2
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &         'center of rotation', 0.0, CENTER, *170)

            CALL GETINT ('rotation', IFLD, INTYP, IFIELD, RFIELD,
     *         NBLK, NRTRAN, D3TRAN, ZGRAD, NEREPL, NNREPL, DIM3,
     *         MAXINT, *170)

         ELSE
            CALL PRTERR ('CMDERR', 'Invalid TWIST Option')
            GOTO 170
         END IF

      ELSE IF (VERB .EQ. 'PROJECT') THEN
         ITRANT = 16
         CALL GETINT ('translation', IFLD, INTYP, IFIELD, RFIELD,
     *      NBLK, NRTRAN, D3TRAN, ZGRAD, NEREPL, NNREPL, DIM3,
     *      MAXINT, *170)

         NUMCOL = -999
         NUMCDM = MAX(1, NUMCOL)

         CALL GETPRO (NEREPL, NNREPL, *170)

      ELSE IF (VERB .EQ. 'OFFSET' .OR. VERB .EQ. 'SHIFT') THEN

C     ... Originally offset just asked for NDIM values.  It was changed
C     to go by axis type (OFFSET Y 1.0). We do maintain compatibility
C     and check for both types of input.

         RMULT = 0.0
         IF (INTYP(IFLD) .EQ. 0) THEN
  930       CONTINUE
            IF (FFEXST (IFLD, INTYP)) THEN
               CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
               IF (MATSTR (WORD, 'RESET', 1)) THEN
                  XOFFS = 0.0
                  YOFFS = 0.0
                  ZOFFS = 0.0
                  RMULT = 0.0
               ELSE IF (MATSTR (WORD, 'ADD', 2)) THEN
C     ... Set for cumulative offsets/shifts
                  RMULT = 1.0
               ELSE IF (MATSTR (WORD, 'ALL', 2)) THEN
                  CALL FFREAL (IFLD, INTYP, RFIELD,
     &                 'model offset', 0.0, TOFFS, *170)
                  XOFFS = RMULT * XOFFS + TOFFS
                  YOFFS = RMULT * YOFFS + TOFFS
                  ZOFFS = RMULT * ZOFFS + TOFFS
               ELSE IF (WORD .EQ. 'X') THEN
                  CALL FFREAL (IFLD, INTYP, RFIELD,
     &                 'X coordinate offset', 0.0, TOFFS, *170)
                  XOFFS = RMULT * XOFFS + TOFFS
               ELSE IF (WORD .EQ. 'Y') THEN
                  CALL FFREAL (IFLD, INTYP, RFIELD,
     &                 'Y coordinate offset', 0.0, TOFFS, *170)
                  YOFFS = RMULT * YOFFS + TOFFS
               ELSE IF (WORD .EQ. 'Z') THEN
                  CALL FFREAL (IFLD, INTYP, RFIELD,
     &                 'Z coordinate offset', 0.0, TOFFS, *170)
                  ZOFFS = RMULT * ZOFFS + TOFFS
               ELSE
                  CALL PRTERR ('CMDERR',
     &              'Expected "X", "Y", "Z", "ALL", "ADD", or "RESET"')
                  GOTO 170
               END IF
               GOTO 930
            END IF
         ELSE
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &           'X coordinate offset', XOFFS, XOFFS, *170)
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &           'Y coordinate offset', YOFFS, YOFFS, *170)
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &           'Z coordinate offset', ZOFFS, ZOFFS, *170)
         END IF

      ELSE IF (VERB .EQ. 'SCALE') THEN
   70    CONTINUE
         IF (FFEXST (IFLD, INTYP)) THEN
            CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
            IF (MATSTR (WORD, 'RESET', 1)) THEN
               XSCAL = 1.0
               YSCAL = 1.0
               ZSCAL = 1.0
            ELSE IF (MATSTR (WORD, 'ALL', 1)) THEN
               CALL FFREAL (IFLD, INTYP, RFIELD,
     &            'model scale factor', 1.0, XSCAL, *80)
               XSCAL = ABS(XSCAL)
               YSCAL = XSCAL
               IF (NDIM .EQ. 3) ZSCAL = XSCAL
            ELSE IF (WORD .EQ. 'X') THEN
               CALL FFREAL (IFLD, INTYP, RFIELD,
     &            'X scale factor', 1.0, XSCAL, *80)
               XSCAL = ABS(XSCAL)
            ELSE IF (WORD .EQ. 'Y') THEN
               CALL FFREAL (IFLD, INTYP, RFIELD,
     &            'Y scale factor', 1.0, YSCAL, *80)
               YSCAL = ABS(YSCAL)
            ELSE IF (WORD .EQ. 'Z') THEN
               CALL FFREAL (IFLD, INTYP, RFIELD,
     &            'Z scale factor', 1.0, ZSCAL, *80)
               ZSCAL = ABS(ZSCAL)
            ELSE
               CALL PRTERR ('CMDERR',
     &            'Expected "X", "Y", "Z", "ALL", or "RESET"')
               GOTO 170
            END IF
            GOTO 70
         END IF

   80    CONTINUE

      ELSE IF (VERB .EQ. 'ZERO') THEN
   90    CONTINUE
         IF (FFEXST (IFLD, INTYP)) THEN
            CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
            IF (MATSTR (WORD, 'RESET', 1)) THEN
               XZERO = 0.0
               YZERO = 0.0
               ZZERO = 0.0
            ELSE IF (WORD .EQ. 'X') THEN
               CALL FFREAL (IFLD, INTYP, RFIELD,
     &            'minimum X coordinate', 0.0, XZERO, *100)
            ELSE IF (WORD .EQ. 'Y') THEN
               CALL FFREAL (IFLD, INTYP, RFIELD,
     &            'minimum Y coordinate', 0.0, YZERO, *100)
            ELSE IF (WORD .EQ. 'Z') THEN
               CALL FFREAL (IFLD, INTYP, RFIELD,
     &            'minimum Z coordinate', 0.0, ZZERO, *100)
            ELSE
               CALL PRTERR ('CMDERR',
     &            'Expected "X", "Y", "Z" or "RESET"')
               GOTO 100
            END IF
            GOTO 90
         END IF

  100    CONTINUE

      ELSE IF (VERB .EQ. 'MIRROR') THEN

  110    CONTINUE
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
               ZMIRR = -1.
            ELSE
               CALL PRTERR ('CMDERR',
     &            'Expected "X", "Y", "Z" or "RESET"')
               GOTO 120
            END IF
            GOTO 110
         END IF

  120    CONTINUE

      ELSE IF (VERB .EQ. 'REVOLVE') THEN
         DEGANG = ATAN2(0.0, -1.0) / 180.0

  130    CONTINUE
         IF (FFEXST (IFLD, INTYP)) THEN
            CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
            IF (MATSTR (WORD, 'RESET', 1)) THEN
               ROT3D = .FALSE.
               CALL INIREA (3*3, 0.0, ROTMAT)
               DO 140 I = 1, 3
                  ROTMAT(I,I) = 1.0
  140          CONTINUE
            ELSE IF ((WORD .EQ. 'X') .OR. (WORD .EQ. 'Y')
     &         .OR. (WORD .EQ. 'Z')) THEN
               CALL FFREAL (IFLD, INTYP, RFIELD,
     &            'angle of rotation', 0.0, DEG, *150)
               ROT3D = .TRUE.
               CALL ROTXYZ (WORD, DEG * DEGANG, ROTMAT)
            ELSE
               CALL PRTERR ('CMDERR',
     &            'Expected "X", "Y", "Z" or "RESET"')
               GOTO 150
            END IF
            GOTO 130
         END IF

  150    CONTINUE

      ELSE IF (VERB .EQ. 'REVCEN') THEN
         CALL MINMAX (NUMNP, XN, XMIN, XMAX)
         CALL MINMAX (NUMNP, YN, YMIN, YMAX)

         CALL FFREAL (IFLD, INTYP, RFIELD,
     &      'X revolution center', XMIN, ROTCEN(1), *170)
         CALL FFREAL (IFLD, INTYP, RFIELD,
     &      'Y revolution center', YMIN, ROTCEN(2), *170)
         CALL FFREAL (IFLD, INTYP, RFIELD,
     &      'Z revolution center', 0.0, ROTCEN(3), *170)

      ELSE IF (VERB .EQ. 'BLOCK') THEN
         CALL USBLK (IFLD, INTYP, CFIELD, IFIELD,
     &      ' ', NELBLK, IDELB, BLKTYP, *170)

         NUMCOL = -999
         NUMCDM = MAX(1, NUMCOL)

      ELSE IF (VERB .EQ. 'CENTER') THEN
         IF (ITRANT .NE. 2) THEN
            CALL PRTERR ('CMDERR', 'Rotation must be selected')
            GOTO 170
         END IF

         if (rgrad .ne. 1.0) then
            rgrad = 1.0
            call prterr ('CMDWARN',
     $           'Gradient must equal 1.0 for center rotations')
         end if

         IF (.NOT. CPOINT) THEN
            IF (MOD (NDEGR, 90) .NE. 0) THEN
               CALL PRTERR ('CMDERR', 'Rotation must be 90, 180, 270'
     &            // ' or 360 degrees')
               GOTO 170
            END IF

            N45 = NINT (DIM3 / 45.0)
            IF (MOD (NEREPL, N45) .NE. 0) THEN
               CALL PRTERR ('CMDERR',
     &            'Number of rotations is incorrect')
               WRITE (*, 250) '(multiple of 2 for 90 degrees,'
     &            , ' 4 for 180, 6 for 270, 8 for 360)'
               GOTO 170
            END IF
         END IF

         CALL USBLK (IFLD, INTYP, CFIELD, IFIELD,
     &      'C', NELBLK, IDELB, BLKTYP, *170)

         ISCENT = .TRUE.
         NUMCOL = -999
         NUMCDM = MAX(1, NUMCOL)

      ELSE IF (VERB .EQ. 'TUNNEL') THEN
c         IF (ITRANT .NE. 1) THEN
c            CALL PRTERR ('CMDERR', 'Translation must be selected')
c            GOTO 170
c         END IF

         CALL FFINTG (IFLD, INTYP, IFIELD,
     &      'block id', 0, ID, *170)
         IELB = LOCINT (ID, NELBLK, IDELB)
         IF (IELB .LE. 0) THEN
            CALL PRTERR ('CMDERR', 'Invalid block id')
            GOTO 170
         END IF

         CALL FFINTG (IFLD, INTYP, IFIELD,
     &      'first tunnel', 2, IPARM1, *170)
         IPARM1 = MAX (IPARM1, 1)
         CALL FFINTG (IFLD, INTYP, IFIELD,
     &      'ending tunnel', 0, IPARM2, *170)
         CALL FFINTG (IFLD, INTYP, IFIELD,
     &      'tunnel increment', 1, IPARM3, *170)
         IF ((IPARM2 .GT. 0) .AND. (IPARM1 .GT. IPARM2)) THEN
            CALL PRTERR ('CMDERR', 'Starting level after ending level')
            GOTO 170
         END IF
         IF (IPARM3 .LE. 0) THEN
            CALL PRTERR ('CMDERR', 'Tunnel increment must be positive')
            GOTO 170
         END IF

         BLKTYP(IELB) = 'T'
         IBPARM(1,IELB) = IPARM1
         IBPARM(2,IELB) = IPARM2
         IBPARM(3,IELB) = IPARM3

         IF (IPARM1 .LE. 1) IPARM1 = IPARM3 + 1
         IF (IPARM1 .GT. NEREPL) THEN
            CALL PRTERR ('CMDWARN', 'Not enough levels are defined')
         END IF

         NUMCOL = -999
         NUMCDM = MAX(1, NUMCOL)

      ELSE IF (VERB .EQ. 'SPECIAL') THEN
         IF (ITRANT .NE. 1) THEN
            CALL PRTERR ('CMDERR', 'Translation must be selected')
            GOTO 170
         END IF
         IF (NEREPL .EQ. NRTRAN(1)) THEN
            CALL PRTERR ('CMDERR',
     &         'Multiple translations must be selected')
            GOTO 170
         END IF

         CALL USBLK (IFLD, INTYP, CFIELD, IFIELD,
     &      'S', NELBLK, IDELB, BLKTYP, *170)

         NUMCOL = -999
         NUMCDM = MAX(1, NUMCOL)

      ELSE IF ((VERB .EQ. 'NSETS') .OR. (VERB .EQ. 'NODESETS')
     &   .OR. (VERB .EQ. 'SSETS') .OR. (VERB .EQ. 'SIDESETS')) THEN
         CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
         IF ((.NOT. MATSTR (WORD, 'FRONT', 1)) .AND.
     &      (.NOT. MATSTR (WORD, 'BACK', 1))) THEN
            CALL PRTERR ('CMDERR', 'Expected FRONT or BACK')
            GOTO 170
         END IF
         ISFRON = MATSTR (WORD, 'FRONT', 1)
         IF ((.NOT. ISFRON)
     &      .AND. (ITRANT .EQ. 2) .AND. (NDEGR .EQ. 360)) THEN
            CALL PRTERR ('CMDERR', 'Back sets are not allowed'
     &         // ' for 360-degree rotation')
            GOTO 170
         END IF

         IF ((VERB .EQ. 'NSETS') .OR. (VERB .EQ. 'NODESETS')) THEN
            IF (ISFRON) THEN
               CALL USIDS (IFLD, INTYP, CFIELD, IFIELD,
     &            NUMNPS, IDNPS, IDNSET(0,2), IDNSET(1,2),
     &            IDNSET(0,1), IDNSET(1,1), MAXSET, *170)
            ELSE
               CALL USIDS (IFLD, INTYP, CFIELD, IFIELD,
     &            NUMNPS, IDNPS, IDNSET(0,1), IDNSET(1,1),
     &            IDNSET(0,2), IDNSET(1,2), MAXSET, *170)
            END IF
         ELSE
            IF (ISFRON) THEN
               CALL USIDS (IFLD, INTYP, CFIELD, IFIELD,
     &            NUMESS, IDESS, IDESET(0,2), IDESET(1,2),
     &            IDESET(0,1), IDESET(1,1), MAXSET, *170)
            ELSE
               CALL USIDS (IFLD, INTYP, CFIELD, IFIELD,
     &            NUMESS, IDESS, IDESET(0,1), IDESET(1,1),
     &            IDESET(0,2), IDESET(1,2), MAXSET, *170)
            END IF
         END IF

      ELSE IF (VERB .EQ. 'LIMITS') THEN
         WRITE (*, *) ' Input Mesh Limits:'
         WRITE (*, 45) 'X', XMIN, 'X', XMAX, XMAX-XMIN
         WRITE (*, 45) 'Y', YMIN, 'Y', YMAX, YMAX-YMIN
         VERB = ' '
 45      FORMAT( '  Minimum ',A1,' = ',1PE12.5,', Maximum ',A1,' = ',
     &        1PE12.5,', Range = ',1PE12.5)
      ELSE IF (VERB .EQ. 'SUMMARY') THEN
         PRINT 160
  160    FORMAT (
     &      /' TRANSLATE,   num_trans, tot_trans, gradient'
     &      /' ROTATE,      num_rotat, tot_rotat, gradient'
     &      , ', center_of_rotation'
     &      /' WARP, POINT, num_trans, tot_trans, gradient'
     &      , ', distance, edge'
     &      /' WARP, axis,  num_trans, tot_trans, gradient'
     &      , ', distance, edge'
     &      /' OFFSET,  xoffset, yoffset, zoffset'
     &      /' MIRROR,  axis, ... '
     &      /' REVOLVE, axis, number_of_degrees, ...'
     &      /' REVCEN,  xcenter, ycenter, zcenter'
     &      /' BLOCK,   id1, id2, ...'
     &      /' CENTER,  id1, id2, ...'
     &      /' TUNNEL,  id, starting_level, ending_level'
     &      , ', level_increment'
     &      /' SPECIAL, id1, id2, ...'
     &      /' NSETS,   FRONT_or_BACK, id1, id2, ...'
     &      /' SSETS,   FRONT_or_BACK, id1, id2, ...'
     &      )
         VERB = ' '

      ELSE IF ((VERB .EQ. 'SHOW') .OR. (VERB .EQ. 'LIST')) THEN
         CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
         IF (MATSTR (WORD, 'COMMANDS', 3)) THEN
            CALL SHOCMD ('COMMANDS', CMDTBL)
         ELSE
            CALL ABRSTR (VERB, WORD, CMDTBL)
            CALL SHOW (VERB, WORD, IDNPS, IDESS, IDNSET, IDESET,
     &         BLKTYP, IBPARM, IDELB, NUMELB, NUMLNK)
         END IF
         VERB = ' '

      ELSE IF (VERB .EQ. 'ATTRIBUT') THEN
         CALL FFINTG (IFLD, INTYP, IFIELD,
     *      'block ID', 0, ID, *170)
         CALL FFREAL (IFLD, INTYP, RFIELD,
     *      'element attribute', 0.0, ATTRIB, *170)
         if (id .lt. 0) then
            call  prterr( 'ERROR', 'Invalid Block ID')
            go to 170
         end if

         imat = locint (id, nelblk, idelb)
         IF (IMAT .EQ. 0) THEN
            CALL PRTERR ('ERROR', 'Invalid Block ID')
         else if (namelb(imat)(:3) .eq. 'BAR'
     &        .or. namelb(imat)(:4) .eq. 'BEAM'
     *        .or. namelb(imat)(:4) .eq. 'TRUS') then
            ELATTR(imat) = ATTRIB
            write (*, 165) id, attrib
 165        FORMAT(1x, 'Attribute for block ',i5,' set to ',1pe10.3)
         else
            CALL PRTERR ('ERROR', 'Block is not a beam, bar, or truss')
         end if
         verb = ' '

      ELSE IF (VERB .EQ. 'CHANGE') THEN
         CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
         CALL FFINTG (IFLD, INTYP, IFIELD,
     *      'old ID', 0, IDOLD, *170)
         CALL FFINTG (IFLD, INTYP, IFIELD,
     *      'new ID', 0, IDNEW, *170)

         IF (MATSTR(WORD,'MATERIAL',1)) THEN
            CALL NEWID ('M', IDELB, NELBLK, IDNEW, IDOLD)
         ELSE IF (MATSTR(WORD,'NODESET',1)) THEN
            CALL NEWID ('N', IDNPS, NUMNPS, IDNEW, IDOLD)
         ELSE IF (MATSTR(WORD,'SIDESET',1)) THEN
            CALL NEWID ('S', IDESS, NUMESS, IDNEW, IDOLD)
         END IF
         VERB = ' '

      ELSE IF (VERB .EQ. 'HELP') THEN
         ISHELP = HELP (MYNAME(:LNAM), 'COMMANDS', CFIELD(IFLD))
         IF (.NOT. ISHELP) CALL SHOCMD ('COMMANDS', CMDTBL)
         VERB = ' '

      ELSE IF ((VERB .EQ. 'END') .OR. (VERB .EQ. 'EXIT')) THEN
         CALL SCNEOF
         GOTO 220

      ELSE IF (VERB .EQ. 'QUIT') THEN
         CALL SCNEOF
         RETURN 1

      ELSE
         CALL PRTERR ('CMDERR', '"' // VERB(:LENSTR(VERB))
     &      // '" is an invalid command')
         VERB = ' '
      END IF

  170 CONTINUE

C   --Reset block types, ANYCEN, and NUMCOL

      IF (NUMCOL .LT. 0) THEN
         IF (ITRANT .NE. 2) THEN
            IF (NEREPL .EQ. NRTRAN(1)) THEN
               DO 180 IELB = 1, NELBLK
                  IF (BLKTYP(IELB) .EQ. 'S') BLKTYP(IELB) = ' '
  180          CONTINUE
            END IF

            DO 190 IELB = 1, NELBLK
               IF (BLKTYP(IELB) .EQ. 'C') BLKTYP(IELB) = ' '
  190       CONTINUE
            ANYCEN = .FALSE.
            NUMCOL = 0
            NUMCDM = MAX(1, NUMCOL)

         ELSE
            DO 200 IELB = 1, NELBLK
               IF (BLKTYP(IELB) .EQ. 'S') BLKTYP(IELB) = ' '
  200       CONTINUE

            ANYCEN = .FALSE.
            DO 210 IELB = 1, NELBLK
               IF (BLKTYP(IELB) .EQ. 'C') ANYCEN = .TRUE.
  210       CONTINUE

            IF (.NOT. ANYCEN) THEN
               NUMCOL = 0
            ELSE IF (CPOINT) THEN
               NUMCOL = 1
            ELSE IF (NUMCOL .NE. -1) THEN
               N45 = NINT (DIM3 / 45.0)
               NUMCOL = NEREPL / N45 + 1
            END IF
            NUMCDM = MAX(1, NUMCOL)
         END IF
      END IF

      IF (VERB .NE. ' ') THEN
         CALL SHOW (VERB, ' ', IDNPS, IDESS, IDNSET, IDESET,
     &      BLKTYP, IBPARM, IDELB, NUMELB, NUMLNK)
      END IF

      GOTO 20

  220 CONTINUE
      IF (ITRANT .EQ. 0) ITRANT = 1

C   --Delete back sets on 360-degree rotations

      IF ((ITRANT .EQ. 2) .AND. (NDEGR .EQ. 360)) THEN
         IF (IDNSET(0,2) .GT. 0) IDNSET(0,2) = 0
         IF (IDESET(0,2) .GT. 0) IDESET(0,2) = 0
      END IF

C   --Delete center block types if invalid

      IF (NUMCOL .EQ. -1) THEN
         DO 230 IELB = 1, NELBLK
            IF (BLKTYP(IELB) .EQ. 'C') BLKTYP(IELB) = ' '
  230    CONTINUE
      END IF

C   --Center of rotation is meaningless if center blocks defined

      IF (NUMCOL .GT. 0) THEN
         CENTER = 0.0
      END IF

C   --Fix up tunnel block type parameters

      DO 240 IELB = 1, NELBLK
         IF (BLKTYP(IELB) .EQ. 'T') THEN
            IF (IBPARM(1,IELB) .GT. NEREPL) THEN
               BLKTYP(IELB) = ' '
            ELSE
               IF (IBPARM(1,IELB) .LE. 1)
     &            IBPARM(1,IELB) = IBPARM(3,IELB) + 1
               IF ((IBPARM(2,IELB) .LE. 0)
     &            .OR. (IBPARM(2,IELB) .GT. NEREPL))
     &            IBPARM(2,IELB) = NEREPL
            END IF
         END IF
  240 CONTINUE

      RETURN
  250 FORMAT (5X, 5A)
      END
