C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C -*- Mode: fortran -*-
C=======================================================================
      SUBROUTINE COMAND (IDNPS, IDESS, IDNSET, IDESET,
     &   IDELB, NUMELB, NUMLNK, NAMELB, ELATTR,
     &  XN, YN, A, *)
C=======================================================================

C   --*** COMAND *** (GENSHELL) Input and process commands
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
C   --   IDELB - IN - the ids for each block
C   --   NUMELB - IN - the number of elements for each block
C   --   NUMLNK - IN - the number of nodes per element for each block
C   --   XN, YN - IN - the 2D nodal coordinates
C   --   * - return statement for QUIT
C   --
C   --Common Variables:
C   --   Sets NDIM, NUMNP, NUMEL, NELBLK,
C   --      NUMNPS, LNPSNL, NUMESS, LESSEL, LESSNL of /DBNUMS/
CCC   --   Sets ITRANT, NNREPL, NEREPL, DIM3, NRTRAN, D3TRAN, ZGRAD,
CCC   --      CENTER, NUMCOL of /PARAMS/
C   --   Sets XOFFS, YOFFS, ZOFFS of /XYZOFF/
C   --   Sets ROT3D, ROTMAT, ROTCEN of /XYZROT/

      PARAMETER (MAXFLD = 32)
      PARAMETER (MAXSET = 10)

      INCLUDE 'gs_dbnums.blk'
      INCLUDE 'gs_params.blk'
      INCLUDE 'gs_xyzoff.blk'
      INCLUDE 'gs_xyzrot.blk'
      INCLUDE 'gs_xyzmir.blk'
      INCLUDE 'gs_xyzero.blk'
      INCLUDE 'gs_xyzscl.blk'
      INCLUDE 'gs_splxyz.blk'

      INTEGER IDNPS(*), IDESS(*)
      INTEGER IDNSET(0:MAXSET,2), IDESET(0:MAXSET,2)
      INTEGER IDELB(*)
      CHARACTER*32 NAMELB(*)
      REAL ELATTR(7, NELBLK)
      INTEGER NUMELB(*)
      INTEGER NUMLNK(*)
      REAL XN(*), YN(*)
      REAL A(*)

      LOGICAL FFEXST, MATSTR, HELP

      CHARACTER*8 WORD, VERB, TMPWRD
      INTEGER     INTYP(MAXFLD+1)
      CHARACTER*8 CFIELD(MAXFLD)
      INTEGER     IFIELD(MAXFLD)
      REAL        RFIELD(MAXFLD)
      LOGICAL ISFRON
      LOGICAL ISHELP

c$$$      CHARACTER*8 CMDTBL(34)
c$$$      SAVE CMDTBL
c$$$C      --CMDTBL - the valid commands table
c$$$
c$$$C   --Command table follows.  Remember to change the dimensioned size when
c$$$C   --changing the table.
c$$$      DATA CMDTBL /
c$$$     1   'TRANSLAT', 'ROTATE  ', 'OFFSET  ', 'REVOLVE ', 'REVCEN  ',
c$$$     2   'WARP    ', 'BLOCK   ', 'CENTER  ', 'TUNNEL  ', 'SPECIAL ',
c$$$     3   'NSETS   ', 'NODESETS', 'SSETS   ', 'SIDESETS', 'MIRROR  ',
c$$$     4   'LIST    ', 'SHOW    ', 'HELP    ', 'ZERO    ', 'SCALE   ',
c$$$     5   'END     ', 'EXIT    ', 'QUIT    ', 'TWIST   ', 'PROJECT ',
c$$$     6   'SUMMARY ', 'INTERVAL', 'ROTCEN  ', 'SPAWN   ', 'SPLINE  ',
c$$$     7   'CHANGE  ', 'TRANSPLI', 'SHIFT   ', '        ' /
c$$$

      CHARACTER*8 CMDTBL(27)
      SAVE CMDTBL
C      --CMDTBL - the valid commands table

C   --Command table follows.  Remember to change the dimensioned size when
C   --changing the table.
      DATA CMDTBL /
     1   'TRANSLAT', 'OFFSET  ', 'REVOLVE ', 'REVCEN  ',
     2   'WARP    ',
     3   'NSETS   ', 'NODESETS', 'SSETS   ', 'SIDESETS', 'MIRROR  ',
     4   'LIST    ', 'SHOW    ', 'HELP    ', 'ZERO    ', 'SCALE   ',
     5   'END     ', 'EXIT    ', 'QUIT    ', 'RANDOMIZ',
     6   'SUMMARY ', 'SPAWN   ', 'SPLINE  ', 'ATTRIBUT',
     7   'CHANGE  ', 'SHIFT   ', 'THICKNES', '        ' /

C   --Initialize

      ITRANT = 0
      DIM3  = 1.0
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

      ROT3D = .FALSE.
      CALL INIREA (3*3, 0.0, ROTMAT)
      DO 10 I = 1, 3
         ROTMAT(I,I) = 1.0
   10 CONTINUE

C     -- Initialize the elattr array here.  These will be used for the
C     -- attributes of any 3d BEAM blocks. User can change using
C     -- ATTRIBUTE command.
      do 15 i = 1, nelblk
        elattr(1,i) = 1.0
        elattr(2,i) = 1.0
        elattr(3,i) = 1.0
        elattr(4,i) = 1.0
        elattr(5,i) = 0.0
        elattr(6,i) = 0.0
        elattr(7,i) = 1.0
 15   continue

      CALL MINMAX (NUMNP, XN, XMIN, XMAX)
      ROTCEN(1) = XMIN
      CALL MINMAX (NUMNP, YN, YMIN, YMAX)
      ROTCEN(2) = YMIN
      ROTCEN(3) = 0.0

      IDNSET(0,1) = 0
      IDNSET(0,2) = 0
      IDESET(0,1) = 0
      IDESET(0,2) = 0

   20 CONTINUE

C   --Read command line

      WRITE (*, *)
      CALL FREFLD (0, 0, 'GenShell> ', MAXFLD,
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
         CALL FFREAL (IFLD, INTYP, RFIELD,
     *      'translation thickness', 1.0, DIM3, *170)

      ELSE IF (VERB .EQ. 'THICKNES') THEN
         if (itrant .eq. 0) then
           call prterr('CMDERR',
     $       'A Mesh Transformation command must be specified')
            goto 170
         else
            CALL FFREAL (IFLD, INTYP, RFIELD,
     *           'thickness', 1.0, DIM3, *170)
         end if

      ELSE IF (VERB .EQ. 'SPLINE') THEN
         ITRANT = 64
         CALL FFREAL (IFLD, INTYP, RFIELD,
     *      'spline thickness', 1.0, DIM3, *170)

         CALL GETSPL (A)

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

         IF (INTYP(IFLD) .EQ. 0) THEN
            CALL FFCHAR (IFLD, INTYP, CFIELD, 'MAP', TMPWRD)
         ELSE
            TMPWRD = 'MAP'
         END IF

         CALL FFREAL (IFLD, INTYP, RFIELD,
     *      'warping thickness', 1.0, DIM3, *170)

         IF (MATSTR(WORD, 'LINE', 1)) THEN
            IWARP = 2

C ... For XAXIS and YAXIS warps, we can either map vertically onto
C     the cylinder, keeping both X and Y coordinates the same as in
C     the 2D mesh.  Or, we can map onto the surface keeping the length
C     of the mesh the same in the 2D and 3D meshes.

         ELSE IF (MATSTR(WORD, 'XAXIS', 1)) THEN
            IF (MATSTR(TMPWRD, 'VERTICAL', 1)) THEN
               IWARP = -3
            ELSE
               IWARP = -1
            END IF

         ELSE IF (MATSTR(WORD, 'YAXIS', 1)) THEN
            IF (MATSTR(TMPWRD, 'VERTICAL', 1)) THEN
               IWARP = -4
            ELSE
               IWARP = -2
            END IF

         ELSE IF (MATSTR(WORD, 'POINT', 1)) THEN
            IWARP = 1
         ELSE IF (MATSTR(WORD, 'ELLIPSE', 1)) THEN
            IWARP = 2
            CALL FFREAL (IFLD, INTYP, RFIELD,
     *         'ellipse Y axis', 0.0, XRAD, *170)
            CALL FFREAL (IFLD, INTYP, RFIELD,
     *         'ellipse Z axis', 0.0, YRAD, *170)
         ELSE
            IWARP = 0
            GO TO 170
         END IF
         CALL FFREAL (IFLD, INTYP, RFIELD,
     *      'warping distance', 0.0, DWARP, *170)

C --- END OF WARP

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
               ZSCAL = XSCAL
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
            ELSE IF (MATSTR (WORD, 'ALL', 2)) THEN
               CALL FFREAL (IFLD, INTYP, RFIELD,
     &              'minimum coordinate', 0.0, XZERO, *100)
               YZERO = XZERO
               ZZERO = XZERO
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
     &            'Expected "X", "Y", "Z", "ALL" or "RESET"')
               GOTO 100
            END IF
            GOTO 90
         END IF

 100     CONTINUE

      ELSE IF (VERB .EQ. 'RANDOMIZ') THEN
 101     CONTINUE
         IF (FFEXST (IFLD, INTYP)) THEN
            CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
            IF (MATSTR (WORD, 'RESET', 1)) THEN
               XRAND = 0.0
               YRAND = 0.0
               ZRAND = 0.0
            ELSE IF (MATSTR (WORD, 'ALL', 1)) THEN
               CALL FFREAL (IFLD, INTYP, RFIELD,
     &              'model random factor', 0.0, TRAND, *141)
               XRAND = ABS(TRAND)
               YRAND = ABS(TRAND)
               ZRAND = ABS(TRAND)
            ELSE IF (WORD .EQ. 'X') THEN
               CALL FFREAL (IFLD, INTYP, RFIELD,
     &              'X random factor', 0.0, TRAND, *141)
               XRAND = ABS(TRAND)
            ELSE IF (WORD .EQ. 'Y') THEN
               CALL FFREAL (IFLD, INTYP, RFIELD,
     &              'Y random factor', 0.0, TRAND, *141)
               YRAND = ABS(TRAND)
            ELSE IF (WORD .EQ. 'Z') THEN
               CALL FFREAL (IFLD, INTYP, RFIELD,
     &              'Z random factor', 0.0, TRAND, *141)
               ZRAND = ABS(TRAND)
            ELSE
               CALL PRTERR ('CMDERR',
     &              'Expected "X", "Y", "Z", "ALL", or "RESET"')
               GOTO 141
            END IF
            GOTO 101
         END IF

 141     CONTINUE

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
               CALL ROTXYZ (WORD(:1), DEG * DEGANG, ROTMAT)
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

      ELSE IF ((VERB .EQ. 'NSETS') .OR. (VERB .EQ. 'NODESETS')
     &   .OR. (VERB .EQ. 'SSETS') .OR. (VERB .EQ. 'SIDESETS')) THEN
         CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
         IF ((.NOT. MATSTR (WORD, 'FRONT', 1)) .AND.
     &      (.NOT. MATSTR (WORD, 'BACK', 1))) THEN
            CALL PRTERR ('CMDERR', 'Expected FRONT or BACK')
            GOTO 170
         END IF
         ISFRON = MATSTR (WORD, 'FRONT', 1)

         IF ((VERB .EQ. 'NSETS') .OR. (VERB .EQ. 'NODESETS')) THEN
            IF (ISFRON) THEN
               CALL USIDS (IFLD, INTYP, CFIELD, IFIELD,
     &            NUMNPS, IDNPS, IDNSET(0,2), IDNSET(1,2),
     &            IDNSET(0,1), IDNSET(1,1), *170)
            ELSE
               CALL USIDS (IFLD, INTYP, CFIELD, IFIELD,
     &            NUMNPS, IDNPS, IDNSET(0,1), IDNSET(1,1),
     &            IDNSET(0,2), IDNSET(1,2), *170)
            END IF
         ELSE
            IF (ISFRON) THEN
               CALL USIDS (IFLD, INTYP, CFIELD, IFIELD,
     &            NUMESS, IDESS, IDESET(0,2), IDESET(1,2),
     &            IDESET(0,1), IDESET(1,1), *170)
            ELSE
               CALL USIDS (IFLD, INTYP, CFIELD, IFIELD,
     &            NUMESS, IDESS, IDESET(0,1), IDESET(1,1),
     &            IDESET(0,2), IDESET(1,2), *170)
            END IF
         END IF

      ELSE IF (VERB .EQ. 'SUMMARY') THEN
         PRINT 160
  160    FORMAT (
     &        /' TRANSLATE,    thickness'
     &        /' WARP, POINT,  thickness, radius'
     &        /' WARP, axis,   thickness, radius'
     $        /' WARP, ELLIPSE,thickness, major, minor'
     $        /' SPLINE,       thickness'
     $        /' THICKNESS,    thickness'
     &        /' MIRROR,       axis, ... '
     &        /' OFFSET,       xoffset, yoffset, zoffset'
     $        /' OFFSET,       axis, offset, ...'
     &        /' REVCEN,       xcenter, ycenter, zcenter'
     &        /' REVOLVE,      axis, number_of_degrees, ...'
     $        /' SCALE,        axis, scale, ...'
     &        /' SHIFT,        axis, shift, ...'
     $        /' ZERO,         axis, minimum, ...'
     &        /' NSETS,        FRONT_or_BACK, id1, id2, ...'
     &        /' SSETS,        FRONT_or_BACK, id1, id2, ...'
     $        /' CHANGE,       MATERIAL|NODESET|SIDESET old_id, new_id'
     &        /' ATTRIBUTE,    id which_attribute new_value'
     &        )
         VERB = ' '

      ELSE IF ((VERB .EQ. 'SHOW') .OR. (VERB .EQ. 'LIST')) THEN
         CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
         IF (MATSTR (WORD, 'COMMANDS', 3)) THEN
            CALL SHOCMD ('COMMANDS', CMDTBL)
         ELSE IF (MATSTR (WORD, 'BLOCK', 3)) THEN
            CALL SHOW ('BLOCK', 'BLOCK', IDNPS, IDESS, IDNSET, IDESET,
     &         IDELB, NAMELB, NUMELB, NUMLNK, ELATTR)
         ELSE
            CALL ABRSTR (VERB, WORD, CMDTBL)
            CALL SHOW (VERB, WORD, IDNPS, IDESS, IDNSET, IDESET,
     &         IDELB, NAMELB, NUMELB, NUMLNK, ELATTR)
         END IF
         VERB = ' '

C ... ATTRIBUTE {id} {which } {value}
      ELSE IF (VERB .EQ. 'ATTRIBUT') THEN
         CALL FFINTG (IFLD, INTYP, IFIELD,
     *      'block ID', 0, ID, *170)
         CALL FFINTG (IFLD, INTYP, IFIELD,
     *      'which attribute', 0, IWHICH, *170)
         CALL FFREAL (IFLD, INTYP, RFIELD,
     *      'element attribute', 0.0, ATTRIB, *170)
         if (id .lt. 0) then
            call  prterr( 'ERROR', 'Invalid Block ID')
            go to 170
         end if
         if (iwhich .lt. 1 .or. iwhich .gt. 7) then
            call  prterr( 'ERROR',
     *       'Invalid Attribute Location. Range = 1..7')
            go to 170
         end if

         imat = locint (id, nelblk, idelb)
         IF (IMAT .EQ. 0) THEN
            CALL PRTERR ('ERROR', 'Invalid Block ID')
         else if (namelb(imat)(:3) .eq. 'BAR'
     &        .or. namelb(imat)(:4) .eq. 'BEAM'
     *        .or. namelb(imat)(:4) .eq. 'TRUS') then
            ELATTR(iwhich, imat) = ATTRIB
            write (*, 165) iwhich, id, attrib
 165        FORMAT(1x, 'Attribute number ',i5,' for block ',i5,
     *        ' set to ',1pe10.3)
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
         ELSE IF (MATSTR(WORD,'NODESETS',1) .OR.
     $           MATSTR(WORD,'NSETS',1)) THEN
            CALL NEWID ('N', IDNPS, NUMNPS, IDNEW, IDOLD)
         ELSE IF (MATSTR(WORD,'SIDESETS',1) .OR.
     $           MATSTR(WORD,'SSETS',1)) THEN
            CALL NEWID ('S', IDESS, NUMESS, IDNEW, IDOLD)
         END IF
         VERB = ' '

      ELSE IF (VERB .EQ. 'HELP') THEN
         ISHELP = HELP ('GENSHELL', 'COMMANDS', CFIELD(IFLD))
         IF (.NOT. ISHELP) CALL SHOCMD ('COMMANDS', CMDTBL)
         VERB = ' '

      ELSE IF (VERB .EQ. 'SPAWN') THEN
          CALL PRTERR ('CMDSPEC', 'Spawn not implemented')

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

      IF (VERB .NE. ' ') THEN
         CALL SHOW (VERB, '        ', IDNPS, IDESS, IDNSET, IDESET,
     &      IDELB, NAMELB, NUMELB, NUMLNK, ELATTR)
      END IF

      GOTO 20

  220 CONTINUE
      IF (ITRANT .EQ. 0) ITRANT = 1

      RETURN
      END
