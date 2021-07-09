C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SHOW (STYP, INTYP, IDNPS, IDESS, IDNSET, IDESET,
     &   IDELB, NAMELB, NUMELB, NUMLNK, ELATTR)
C=======================================================================

C   --*** SHOW *** (GEN3D) Display information
C   --   Written by Amy Gilkey - revised 03/07/88
C   --
C   --SHOW displays information about the database and the plot parameters.
C   --
C   --The SHOW options with the items they display are:
C   --   TRANSLAT - number of translations or rotations
C   --   ROTATE   -
C   --   EXPROTAT -
C   --   WARP     -
C   --   TWIST    -
C   --   PROJECT  -
C   --   OFFSET   - generated mesh offsets
C   --   REVOLVE  - generated mesh rotation matrix
C   --   REVCEN   - generated mesh rotation center
C   --   MIRROR   - axes about which mesh is reflected
C   --   BLOCK    - information on all element blocks
C   --   CENTER   -    (number of elements, defined type, etc)
C   --   TUNNEL   -
C   --   NSETS    - IDs of existing and front and back node sets
C   --   SSETS    - IDs of existing and front and back side sets
C   --   VARS     - database title and initial variables
C   --
C   --Parameters:
C   --   STYP - IN - the exact SHOW option string
C   --   INTYP - IN - the abbreviated SHOW option string, ' ' if exact
C   --   IDNPS - IN - the IDs of existing node sets; (0) = length
C   --   IDESS - IN - the IDs of existing side sets; (0) = length
C   --   IDNSET - OUT - the IDs of the front and back surface node sets;
C   --      (0) = length
C   --   IDESET - OUT - the IDs of the front and back surface side sets;
C   --      (0) = length
C   --   IDELB - IN - the ids for each block
C   --   NUMELB - IN - the number of elements for each block
C   --   NUMLNK - IN - the number of nodes per element for each block
C   --
C   --Common Variables:
C   --   Uses NDBIN, NDBOUT of /DBASE/
C   --   Uses TITLE of /DBTITL/
C   --   Uses NDIM, NUMNP, NUMEL, NELBLK,
C   --      NUMNPS, LNPSNL, NUMESS, LESSEL, LESSNL of /DBNUMS/
C   --   Uses DOTRAN, NNREPL, NEREPL, DIM3, NRTRAN, D3TRAN, ZGRAD,
C   --      CENTER, NUMCOL of /PARAMS/
C   --   Uses XOFFS, YOFFS, ZOFFS of /XYZOFF/
C   --   Uses ROT3D, ROTMAT of /XYZROT/

      include 'exodusII.inc'
      INCLUDE 'gs_dbase.blk'
      INCLUDE 'gs_dbtitl.blk'
      INCLUDE 'gs_dbnums.blk'
      INCLUDE 'gs_params.blk'
      INCLUDE 'gs_xyzoff.blk'
      INCLUDE 'gs_xyzrot.blk'
      INCLUDE 'gs_xyzmir.blk'
      INCLUDE 'gs_xyzscl.blk'
      INCLUDE 'gs_xyzero.blk'

      CHARACTER*8 STYP
      CHARACTER*8 INTYP
      INTEGER IDNPS(*), IDESS(*)
      INTEGER IDNSET(0:10,2), IDESET(0:10,2)
      INTEGER IDELB(*)
      INTEGER NUMELB(*)
      CHARACTER*32 NAMELB(*)
      INTEGER NUMLNK(*)
      REAL    ELATTR(7,*)

      CHARACTER*8 SHOTYP
      CHARACTER*80 STRING
      REAL RNUM(9)
      CHARACTER*20 RSTR(9)
      CHARACTER*32 STRB

      CHARACTER*8 SHOTBL(4)
      SAVE SHOTBL
C      --SHOTBL - the special save option table

      DATA SHOTBL /
     &   'VARS    ', 'BLOCK   ', 'ATTRIBUT', '        ' /

C   --Determine the show option

      IF ((INTYP .EQ. ' ') .OR. (STYP .EQ. INTYP)) THEN
         SHOTYP = STYP
      ELSE
         CALL ABRSTR (SHOTYP, INTYP, SHOTBL)
         IF ((STYP .NE. ' ')
     &      .AND. ((SHOTYP .EQ. ' ') .OR. (SHOTYP .NE. INTYP)))
     &      SHOTYP = STYP
      END IF

      IF ((    SHOTYP .EQ. 'TRANSLAT')
     &   .OR. (SHOTYP .EQ. 'INTERVAL')
     &   .OR. (SHOTYP .EQ. 'ROTATE')
     &   .OR. (SHOTYP .EQ. 'EXPROTAT')
     &   .OR. (SHOTYP .EQ. 'WARP')
     &   .OR. (SHOTYP .EQ. 'TWIST')
     &   .OR. (SHOTYP .EQ. 'PROJECT')
     &   .OR. (SHOTYP .EQ. 'SPLINE')
     &   .OR. (SHOTYP .EQ. 'CPOINT')) THEN

         IF (ITRANT .NE. 2) THEN
            CALL SHWINT (ITRANT, DIM3)
         END IF

         IF (ITRANT .EQ. 1) THEN
            CONTINUE
         ELSE IF (ITRANT .EQ. 4) THEN
            CALL NUMSTR1 (3, DWARP, RSTR(1), LR1)
            IF (IWARP .EQ.  1) STRB = 'Point'
            IF (IWARP .EQ. -1) STRB = 'X Axis, Map to surface'
            IF (IWARP .EQ. -2) STRB = 'Y Axis, Map to surface'
            IF (IWARP .EQ. -3) STRB = 'X Axis, Vertical projection'
            IF (IWARP .EQ. -4) STRB = 'Y Axis, Vertical projection'
            IF (IWARP .EQ. 1) THEN
               WRITE (*, 130) '   At a distance of ', RSTR(1)(:LR1),
     *            ' from the ',STRB(:LENSTR(STRB)),
     *            ' X = 0.0, Y = 0.0, Z = ',RSTR(1)(:LR1)
            ELSE
               WRITE (*, 130) '   At a distance of ', RSTR(1)(:LR1),
     *            ' from the ',STRB(:LENSTR(STRB))
            END IF
         ELSE IF (ITRANT .EQ. 8) THEN
            CONTINUE
         ELSE IF (ITRANT .EQ. 16) THEN
            CONTINUE
         ELSE IF (ITRANT .EQ. 64) THEN
            CONTINUE
         END IF

      ELSE IF (SHOTYP .EQ. 'MIRROR') THEN
         IF (XMIRR .EQ. -1.) THEN
            STRING = 'New X = - Old X'
            WRITE (*, 130) '   ', STRING(:LENSTR(STRING))
         END IF
         IF (YMIRR .EQ. -1.) THEN
            STRING = 'New Y = - Old Y'
            WRITE (*, 130) '   ', STRING(:LENSTR(STRING))
         END IF
         IF (ZMIRR .EQ. -1.) THEN
            STRING = 'New Z = - Old Z'
            WRITE (*, 130) '   ', STRING(:LENSTR(STRING))
         END IF

      ELSE IF (SHOTYP .EQ. 'OFFSET' .OR. SHOTYP .EQ. 'SHIFT') THEN
         RNUM(1) = XOFFS
         RNUM(2) = YOFFS
         RNUM(3) = ZOFFS
         CALL NUMSTR (3, 4, RNUM, RSTR, LR)
         WRITE (*, 130)
     &      'Coordinate offsets =', (' ', RSTR(I)(:LR), I=1,3)

      ELSE IF (SHOTYP .EQ. 'ZERO') THEN
         RNUM(1) = XZERO
         RNUM(2) = YZERO
         RNUM(3) = ZZERO
         CALL NUMSTR (3, 4, RNUM, RSTR, LR)
         WRITE (*, 130) 'Minimum nonzero coordinates =',
     *      (' ', RSTR(I)(:LR), I=1,3)

      ELSE IF (SHOTYP .EQ. 'SCALE') THEN
         RNUM(1) = XSCAL
         RNUM(2) = YSCAL
         RNUM(3) = ZSCAL
         CALL NUMSTR (3, 4, RNUM, RSTR, LR)
         WRITE (*, 130)
     &      'Coordinate scale factors =', (' ', RSTR(I)(:LR), I=1,3)

      ELSE IF (SHOTYP .EQ. 'RANDOMIZ') THEN
         RNUM(1) = XRAND
         RNUM(2) = YRAND
         RNUM(3) = ZRAND
         CALL NUMSTR (3, 4, RNUM, RSTR, LR)
         WRITE (*, 130)
     &      'Coordinate random factors =', (' ', RSTR(I)(:LR),I=1,3)

      ELSE IF (SHOTYP .EQ. 'REVOLVE') THEN
         IF (ROT3D) THEN
            WRITE (*, 130) 'Rotation matrix for generated mesh:'
            DO 30 I = 1, 3
               IX = (I-1) * 3
               DO 20 J = 1, 3
                  RNUM(IX+J) = ROTMAT(I,J)
   20          CONTINUE
   30       CONTINUE
            CALL NUMSTR (9, 4, RNUM, RSTR, LR)
            DO 40 I = 1, 3
               IX = (I-1) * 3
               WRITE (*, 130) ('   ', RSTR(IX+J)(:LR), J=1,3)
   40       CONTINUE
         ELSE
            WRITE (*, 130) 'No rotation defined for generated mesh'
         END IF

      ELSE IF (SHOTYP .EQ. 'REVCEN') THEN
         CALL NUMSTR (3, 4, ROTCEN, RSTR, LR)
         WRITE (*, 130)
     &      'Center of revolution =', (' ', RSTR(I)(:LR), I=1,3)
         IF (.NOT. ROT3D) THEN
            WRITE (*, 130) 'No revolution defined for generated mesh'
         END IF

      ELSE IF (SHOTYP .EQ. 'ROTCEN') THEN
         CALL NUMSTR1 (4, CENTER, RSTR, LR)
         WRITE (*, 130)
     &      'Center of rotation = ',RSTR(1)(:LR)

      ELSE IF (SHOTYP(:5) .EQ. 'BLOCK') THEN
         DO 70 IELB = 1, NELBLK
           WRITE (*, 60, IOSTAT=IDUM)
     &       IDELB(IELB), NUMELB(IELB),
     *       NUMLNK(IELB), NAMELB(IELB)(:LENSTR(NAMELB(IELB)))
 60        FORMAT (1X, 'Block', I6, ':',
     &       I6, ' elements', I4, '-node ', A)
 70      CONTINUE

      ELSE IF (SHOTYP .EQ. 'ATTRIBUT') THEN
         DO 90 IELB = 1, NELBLK
           WRITE (*, 80, IOSTAT=IDUM)
     &       IDELB(IELB), (ELATTR(I,IELB),I=1,7)
 80        FORMAT (1X, 'Block', I6, ':',
     &       7(1pE10.3))
 90      CONTINUE

      ELSE IF ((SHOTYP .EQ. 'NSETS') .OR. (SHOTYP .EQ. 'NODESETS')) THEN
         WRITE (*, 140, IOSTAT=IDUM)
     &      'IDs for INPUT node sets',
     &      (IDNPS(I), I=1,NUMNPS)
         WRITE (*, 140, IOSTAT=IDUM)
     &      'IDs for FRONT node sets',
     &      (IDNSET(I,1), I=1,IDNSET(0,1))
         WRITE (*, 140, IOSTAT=IDUM)
     &      'IDs for BACK node sets',
     &      (IDNSET(I,2), I=1,IDNSET(0,2))

      ELSE IF ((SHOTYP .EQ. 'SSETS') .OR. (SHOTYP .EQ. 'SIDESETS')) THEN
         WRITE (*, 140, IOSTAT=IDUM)
     &      'IDs for INPUT side sets',
     &      (IDESS(I), I=1,NUMESS)
         WRITE (*, 140, IOSTAT=IDUM)
     &      'IDs for FRONT side sets',
     &      (IDESET(I,1), I=1,IDESET(0,1))
         WRITE (*, 140, IOSTAT=IDUM)
     &      'IDs for BACK side sets',
     &      (IDESET(I,2), I=1,IDESET(0,2))

      ELSE IF (SHOTYP .EQ. 'VARS') THEN
         CALL DBPINI ('NTIS', NDBIN, TITLE, NDIM, NUMNP, NUMEL, NELBLK,
     &      NUMNPS, LNPSNL, LNPSNL, NUMESS, LESSEL, LESSNL, LESSNL,
     &      IDUM, IDUM, IDUM)

      ELSE
         CALL PRTERR ('CMDWARN', 'Invalid SHOW option')
      END IF

      RETURN
  130 FORMAT (1X, 10A)
  140 FORMAT (1X, A, ':', T36, 5I6, :, /, (1X, 15I6))
      END
