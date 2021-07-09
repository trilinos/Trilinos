C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SHOW (STYP, INTYP, IDNPS, IDESS, IDNSET, IDESET,
     &   BLKTYP, IBPARM, IDELB, NUMELB, NUMLNK)
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
C   --   BLKTYP - IN - the element block type
C   --   IBPARM - IN - the block parameters (defined by the block type)
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

      INCLUDE 'g3_dbase.blk'
      INCLUDE 'g3_dbtitl.blk'
      INCLUDE 'g3_dbnums.blk'
      INCLUDE 'g3_params.blk'
      INCLUDE 'g3_xyzoff.blk'
      INCLUDE 'g3_xyzrot.blk'
      INCLUDE 'g3_cmdsho.blk'
      INCLUDE 'g3_xyzmir.blk'
      INCLUDE 'g3_xyzscl.blk'
      INCLUDE 'g3_xyzero.blk'

      CHARACTER*(*) STYP
      CHARACTER*(*) INTYP
      INTEGER IDNPS(*), IDESS(*)
      INTEGER IDNSET(0:MAXSET,2), IDESET(0:MAXSET,2)
      CHARACTER BLKTYP(*)
      INTEGER IBPARM(4,*)
      INTEGER IDELB(*)
      INTEGER NUMELB(*)
      INTEGER NUMLNK(*)

      CHARACTER*8 SHOTYP
      CHARACTER*80 STRING
      REAL RNUM(9)
      CHARACTER*20 RSTR(9)
      CHARACTER*20 STRA, STRB

      CHARACTER*8 SHOTBL(28)
      SAVE SHOTBL
C      --SHOTBL - the special save option table

      DATA SHOTBL /
     &  'TRANSLAT', 'INTERVAL', 'ROTATE  ', 'EXPROTAT', 'WARP    ',
     &  'TWIST   ', 'PROJECT ', 'SPLINE  ', 'CPOINT  ', 'MIRROR  ',
     &  'OFFSET  ', 'SHIFT   ', 'ZERO    ', 'SCALE   ', 'REVOLVE ',
     &  'REVCEN  ', 'ROTCEN  ', 'ROTAXIS ', 'BLOCK   ', 'TUNNEL  ',
     &  'CENTER  ', 'SPECIAL ', 'NSETS   ', 'NODESETS', 'SSETS   ',
     &  'SIDESETS', 'VARS    ', '        ' /

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
            CALL SHWINT (ITRANT, NEREPL, DIM3, NRTRAN, D3TRAN, ZGRAD)
         END IF

         IF (ITRANT .EQ. 1) THEN
            CONTINUE
         ELSE IF (ITRANT .EQ. 2) THEN
            CALL INTSTR (1, 0, NEREPL, STRA, LSTRA)
            RNUM(1) = DIM3
            RNUM(2) = RGRAD
            CALL NUMSTR (2, 4, RNUM, RSTR(1), LR)

            IF (ABS (RGRAD - 1.0) .LE. 1.0E-6) THEN
               WRITE (*, 130) 'Rotate mesh ', STRA(:LSTRA),
     &            ' times for a total of ', RSTR(1)(:LR), ' degrees'
            ELSE
               WRITE (*, 130) 'Rotate mesh ', STRA(:LSTRA),
     &            ' times for a total of ', RSTR(1)(:LR), ' degrees',
     &            ' with a gradient of ', RSTR(2)(:LR)
            END IF

            if (rotax .eq. 0) then
              write (*, 130) '   About the Y Axis'
            else
              write (*, 130) '   About the X Axis'
            end if

            IF (CPOINT) THEN
               IF (NUMCOL .EQ. 1) THEN
                  WRITE (*, 130) '   Single point of rotation'
               ELSE
                  WRITE (*, 130) '   Single point of rotation'
     &               , ' (undefined)'
               END IF
            ELSE IF (NUMCOL .GT. 0) THEN
               CALL INTSTR (1, 0, NUMCOL, STRA, LSTRA)
               WRITE (*, 130) '   Center of rotation in ',
     &            STRA(:LSTRA), ' columns'
            ELSE
               CALL NUMSTR1 (4, CENTER, RSTR(1), LR)
               WRITE (*, 130) '   Center of rotation = ', RSTR(1)(:LR)
            END IF

         ELSE IF (ITRANT .EQ. 4) THEN
            CALL NUMSTR1 (3, DWARP, RSTR(1), LR1)
            IF (IWARP .EQ.  1) STRB = 'Point'
            IF (IWARP .EQ. -1) STRB = 'X Axis'
            IF (IWARP .EQ. -2) STRB = 'Y Axis'
            IF (IWARP .EQ. 1) THEN
               WRITE (*, 130) '   At a distance of ', RSTR(1)(:LR1),
     *            ' from the ',STRB(:LENSTR(STRB)),
     *            ' X = 0.0, Y = 0.0, Z = ',RSTR(1)(:LR1)
            ELSE
               WRITE (*, 130) '   At a distance of ', RSTR(1)(:LR1),
     *            ' from the ',STRB(:LENSTR(STRB))
            END IF
            IF (VEDGE) THEN
               WRITE (*, 130) '   With vertical edges'
            ELSE
               WRITE (*, 130) '   With radial edges'
            END IF
         ELSE IF (ITRANT .EQ. 8) THEN
            CONTINUE
         ELSE IF (ITRANT .EQ. 16) THEN
            CONTINUE
         ELSE IF (ITRANT .EQ. 32) THEN
            IF (CPOINT) THEN
               IF (NUMCOL .EQ. 1) THEN
                  WRITE (*, 130) '   Single point of rotation'
               ELSE
                  WRITE (*, 130) '   Single point of rotation'
     &               , ' (undefined)'
               END IF
            ELSE IF (NUMCOL .GT. 0) THEN
               CALL INTSTR (1, 0, NUMCOL, STRA, LSTRA)
               WRITE (*, 130) '   Center of rotation in ',
     &            STRA(:LSTRA), ' columns'
            ELSE
               CALL NUMSTR1 (4, CENTER, RSTR(1), LR)
               WRITE (*, 130) '   Center of rotation = ', RSTR(1)(:LR)
            END IF
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

      ELSE IF (SHOTYP .EQ. 'ROTAXIS') THEN
            if (rotax .eq. 0) then
              write (*, 130) 'Rotate about the Y Axis'
            else
              write (*, 130) 'Rotate about the X Axis'
            end if

      ELSE IF ((SHOTYP .EQ. 'BLOCK') .OR. (SHOTYP .EQ. 'CENTER')
     &   .OR. (SHOTYP .EQ. 'TUNNEL') .OR. (SHOTYP .EQ. 'SPECIAL'))
     &   THEN

         NCEN = 0
         NTUN = 0
         NSPEC = 0
         DO 70 IELB = 1, NELBLK
            WRITE (STRA, 50, IOSTAT=IDUM) IELB
 50         FORMAT ('(#', I5, ')')
            CALL PCKSTR (1, STRA)
            LSTRA = LENSTR (STRA)
            WRITE (*, 60, IOSTAT=IDUM)
     &         IDELB(IELB), STRA(:LSTRA), NUMELB(IELB), NUMLNK(IELB)
   60       FORMAT (1X, 'Block', I6, 1X, A, ':',
     &         I6, ' elements', I4, '-node')
            IF (BLKTYP(IELB) .EQ. ' ') THEN
               CONTINUE
            ELSE IF (BLKTYP(IELB) .EQ. 'C') THEN
               NCEN = NCEN + 1
            ELSE IF (BLKTYP(IELB) .EQ. 'T') THEN
               NTUN = NTUN + 1
            ELSE IF (BLKTYP(IELB) .EQ. 'S') THEN
               NSPEC = NSPEC + 1
            END IF
   70    CONTINUE

         IF ((NCEN .GT. 0) .AND. (NUMCOL .GT. 0)) THEN
            WRITE (*, *)
            WRITE (*, 130) 'Rotation Center Block IDs:'
            N = 0
            STRING = ' '
            DO 80 IELB = 1, NELBLK
               IF (BLKTYP(IELB) .EQ. 'C') THEN
                  N = N + 1
                  WRITE (STRING((N-1)*6+1:N*6), '(I6)', IOSTAT=IDUM)
     &               IDELB(IELB)
                  IF (N .GE. 12) THEN
                     WRITE (*, 130) '   ', STRING(:LENSTR(STRING))
                     N = 0
                     STRING = ' '
                  END IF
               END IF
   80       CONTINUE
            IF (N .GT. 0) THEN
               WRITE (*, 130) '   ', STRING(:LENSTR(STRING))
            END IF
         END IF

         IF (NTUN .GT. 0) THEN
            WRITE (*, *)
            DO 110 IELB = 1, NELBLK
               IF (BLKTYP(IELB) .EQ. 'T') THEN
                  IF (IBPARM(2,IELB) .GT. 0) THEN
                     WRITE (STRB, '(I5)', IOSTAT=IDUM) IBPARM(2,IELB)
                     WRITE (STRING, 90)
     &                  IBPARM(1,IELB), STRB(:5), IBPARM(3,IELB)
                  ELSE
                     WRITE (STRING, 90)
     &                  IBPARM(1,IELB), 'end', IBPARM(3,IELB)
                  END IF
   90             FORMAT ('Tunnel from step ', I5, ' to ', A,
     &               ', changes every ', I5, ' step')
                  CALL SQZSTR (STRING, LSTR)
                  WRITE (*, 100) IDELB(IELB), STRING(:LSTR)
  100             FORMAT (1X, 'Block', I5, ' (ID):', 2X, A)
               END IF
  110       CONTINUE
         END IF

         IF (NSPEC .GT. 0) THEN
            WRITE (*, *)
            WRITE (*, 130) 'Special Block IDs:'
            N = 0
            STRING = ' '
            DO 120 IELB = 1, NELBLK
               IF (BLKTYP(IELB) .EQ. 'S') THEN
                  N = N + 1
                  WRITE (STRING((N-1)*6+1:N*6), '(I6)', IOSTAT=IDUM)
     &               IDELB(IELB)
                  IF (N .GE. 12) THEN
                     WRITE (*, 130) '   ', STRING(:LENSTR(STRING))
                     N = 0
                     STRING = ' '
                  END IF
               END IF
  120       CONTINUE
            IF (N .GT. 0) THEN
               WRITE (*, 130) '   ', STRING(:LENSTR(STRING))
            END IF
         END IF

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
         CALL DBPINI ('TIS', NDBIN, TITLE, NDIM, NUMNP, NUMEL, NELBLK,
     &      NUMNPS, LNPSNL, LNPSNL, NUMESS, LESSEL, LESSNL, LESSNL,
     &      IDUM, IDUM, IDUM, ' ')

      ELSE IF (SHOTYP .EQ. ' ') THEN
        CALL SHOCMD ('LIST/SHOW Options', SHOTBL)
      ELSE
         CALL PRTERR ('CMDWARN', 'Invalid SHOW option')
      END IF

      RETURN
  130 FORMAT (1X, 10A)
  140 FORMAT (1X, A, ':', T36, 5I6, :, /, (1X, 15I6))
      END
