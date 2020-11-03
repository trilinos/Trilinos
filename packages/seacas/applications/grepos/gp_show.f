C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C -*- Mode: fortran -*-
C=======================================================================
      SUBROUTINE SHOW (STYP, INTYP)
C=======================================================================

C     --*** SHOW *** (GEN3D) Display information
C     --   Written by Amy Gilkey - revised 03/07/88
C     --
C     --SHOW displays information about the database and the plot parameters.
C     --
C     --The SHOW options with the items they display are:
C     --   OFFSET   - generated mesh offsets
C     --   REVOLVE  - generated mesh rotation matrix
C     --   REVCEN   - generated mesh rotation center
C     --   MIRROR   - axes about which mesh is reflected
C     --   VARS     - database title and initial variables
C     --
C     --Parameters:
C     --   STYP - IN - the exact SHOW option string
C     --   INTYP - IN - the abbreviated SHOW option string, ' ' if exact
C     --
C     --Common Variables:
C     --   Uses NDBIN, NDBOUT of /DBASE/
C     --   Uses TITLE of /DBTITL/
C     --   Uses NDIM, NUMNP, NUMEL, NELBLK,
C     --      NUMNPS, LNPSNL, NUMESS, LESSEL, LESSNL of /DBNUMS/
C     --   Uses XOFFS, YOFFS, ZOFFS of /XYZOFF/
C     --   Uses ROT3D, ROTMAT of /XYZROT/

      include 'gp_dbase.blk'
      include 'gp_dbtitl.blk'
      include 'gp_dbnums.blk'
      include 'gp_xyzoff.blk'
      include 'gp_xyzrot.blk'
      include 'gp_xyzmir.blk'
      include 'gp_xyzero.blk'
      include 'gp_xyzscl.blk'
      include 'gp_xyzwrp.blk'
      include 'gp_smooth.blk'
      include 'gp_snap.blk'
      include 'gp_combine.blk'
      include 'gp_deform.blk'

      CHARACTER*(*) STYP
      CHARACTER*(*) INTYP

      CHARACTER*32 SHOTYP, SMTYP
      CHARACTER*80 STRING
      REAL RNUM(9)
      CHARACTER*20 RSTR(12)
      CHARACTER*20 STRA, STRB

      CHARACTER*32 SHOTBL(2)
      SAVE SHOTBL
C     --SHOTBL - the special save option table

      DATA SHOTBL /
     &  'VARS    ', '        ' /

C     --Determine the show option

      IF ((INTYP .EQ. ' ') .OR. (STYP .EQ. INTYP)) THEN
        SHOTYP = STYP
      ELSE
        CALL ABRSTR (SHOTYP, INTYP, SHOTBL)
        IF ((STYP .NE. ' ')
     &    .AND. ((SHOTYP .EQ. ' ') .OR. (SHOTYP .NE. INTYP)))
     &    SHOTYP = STYP
      END IF

      IF (SHOTYP .EQ. 'EXECUTE') THEN
        WRITE (*, 40) 'Executing currently specified transformations'

      ELSE IF (SHOTYP .EQ. 'MIRROR') THEN
        IF (XMIRR .EQ. -1.) THEN
          STRING = 'New X = - Old X'
          WRITE (*, 40) '   ', STRING(:LENSTR(STRING))
        ELSE
          STRING = 'New X =   Old X'
          WRITE (*, 40) '   ', STRING(:LENSTR(STRING))
        END IF
        IF (YMIRR .EQ. -1.) THEN
          STRING = 'New Y = - Old Y'
          WRITE (*, 40) '   ', STRING(:LENSTR(STRING))
        ELSE
          STRING = 'New Y =   Old Y'
          WRITE (*, 40) '   ', STRING(:LENSTR(STRING))
        END IF
        IF (NDIM .EQ. 3) THEN
          IF (ZMIRR .EQ. -1.) THEN
            STRING = 'New Z = - Old Z'
            WRITE (*, 40) '   ', STRING(:LENSTR(STRING))
          ELSE
            STRING = 'New Z =   Old Z'
            WRITE (*, 40) '   ', STRING(:LENSTR(STRING))
          END IF
        END IF
      ELSE IF (SHOTYP .EQ. 'OFFSET' .or. SHOTYP .EQ. 'SHIFT') THEN
        IF (SPLOFF) THEN
          WRITE (*,*) 'Spline offset selected'
        ELSE
          RNUM(1) = XOFFS
          RNUM(2) = YOFFS
          RNUM(3) = ZOFFS
          CALL NUMSTR (NDIM, 4, RNUM, RSTR, LR)
          WRITE (*, 40)
     &      'Coordinate offsets =', (' ', RSTR(I)(:LR), I=1,NDIM)
        END IF

      ELSE IF (SHOTYP .EQ. 'SCALE') THEN
        RNUM(1) = XSCAL
        RNUM(2) = YSCAL
        RNUM(3) = ZSCAL
        CALL NUMSTR (NDIM, 4, RNUM, RSTR, LR)
        WRITE (*, 40)
     &    'Coordinate scale factors =', (' ', RSTR(I)(:LR),I=1,NDIM)

      ELSE IF (SHOTYP .EQ. 'RANDOMIZE') THEN
        RNUM(1) = XRAND
        RNUM(2) = YRAND
        RNUM(3) = ZRAND
        CALL NUMSTR (NDIM, 4, RNUM, RSTR, LR)
        WRITE (*, 40)
     &    'Coordinate random factors =', (' ', RSTR(I)(:LR),I=1,NDIM)

      ELSE IF (SHOTYP .EQ. 'ZERO') THEN
        RNUM(1) = XZERO
        RNUM(2) = YZERO
        RNUM(3) = ZZERO
        CALL NUMSTR (NDIM, 4, RNUM, RSTR, LR)
        WRITE (*, 40)
     &    'Minimum nonzero coordinates =', (' ', RSTR(I)(:LR),
     &    I=1,NDIM)

      ELSE IF (SHOTYP .EQ. 'REVOLVE' .OR. SHOTYP .EQ. 'ROTATE') THEN
        IF (ROT3D) THEN
          WRITE (*, 40) 'Rotation matrix for generated mesh:'
          DO I = 1, 3
            write (*,50) (rotmat(i,j),j=1,3)
          end do
        ELSE
          WRITE (*, 40) 'No rotation defined for generated mesh'
        END IF

      ELSE IF (SHOTYP .EQ. 'REVCEN') THEN
        CALL NUMSTR (NDIM, 4, ROTCEN, RSTR, LR)
        WRITE (*, 40)
     &    'Center of revolution =', (' ', RSTR(I)(:LR), I=1,NDIM)
        IF (.NOT. ROT3D) THEN
          WRITE (*, 40) 'No rotation defined for generated mesh'
        END IF

      ELSE IF (SHOTYP .EQ. 'DEFORM') THEN
        if (idefst .gt. 0) then
          call intstr(1, 0, IDEFST, STRA,  LR1)
          WRITE (*, 40)
     &      'Deform at step ', STRA(:LR1),
     *      '. Displacements at that step will be set to zero.'
        else
          write (*,*) 'Deformation turned off (RESET)'
        end if

      ELSE IF (SHOTYP .EQ. 'VARS') THEN
        CALL DBPINI ('TIS', NDBIN, TITLE, NDIM, NUMNP, NUMEL, NELBLK,
     &    NUMNPS, LNPSNL, LNPSNL, NUMESS, LESSEL, LESSNL, LESSNL,
     &    IDUM, IDUM, IDUM, ' ')

      ELSE IF (SHOTYP .EQ. 'SMOOTH') THEN
        write (*, 40)
     $    'Smoothing Type = LAPLACIAN'
        call numstr1(4, TOLER, RSTR(1), LR1)
        call intstr (1, 0, NIT,   STRA,    LR2)
        call numstr1(4, R0,    RSTR(3), LR3)
        write (*, 40)
     $    'Tolerance      = ',RSTR(1)(:LR1)
        write (*, 40)
     $    'Iterations     = ',STRA(:LR2)
        write (*, 40)
     $    'Relaxation Par = ',RSTR(3)(:LR3)

      ELSE IF (SHOTYP .EQ. 'USERSUBROUTINE') THEN
        write (*, 40) 'Coordinate modification by ',
     *    'user supplied subroutine (xyzmod).'

      ELSE IF (SHOTYP .EQ. 'CENTROIDS') THEN
        write (*, 40) 'Element centroids will be calculated and output.'

      ELSE IF (SHOTYP .EQ. 'EQUIVALENCE') THEN
         IF (EQUIV) THEN
            CALL NUMSTR1(6, EQTOLER, RSTR, LR)
            WRITE (*, 40)
     &           'Node Equivalence Tolerance = ', RSTR(1)(:LR)
         ELSE
            WRITE (*, 40)
     &           'Node Equivalencing turned off (RESET)'
         END IF

      ELSE IF (SHOTYP .EQ. 'SNAP' .OR. SHOTYP .EQ. 'MOVE') THEN
        if (numsnp .gt. 0) then
          do i=1, numsnp
            call intstr(1, 0, IDSSSL(i), STRA,  LR1)
            call intstr(1, 0, IDSSMA(i), STRB,  LR2)
            call numstr1(4, snptol(i), RSTR(4), LR4)
            call numstr1(4, delmax(i), RSTR(5), LR5)
            if (usnorm(i) .eq. PNORM) then
              string = 'normal to slave surf'
            else if (usnorm(i) .eq. PRAD) then
              string = 'radially from '
            else if (usnorm(i) .eq. PVECT) then
              string = 'along vector'
            else if (usnorm(i) .eq. PEDGE) then
              string = 'along element edge'
            end if

            if (ismtyp(i) .eq. IMOVE) then
              SMTYP = 'Move'
            else if (ismtyp(i) .eq. ISNAP) then
              SMTYP = 'Snap'
            else
              SMTYP = 'ERROR'
            end if

            if (usnorm(i) .eq. PNORM .or. usnorm(i) .eq. PEDGE) then
              write (*, 40) SMTYP(:4), ' Sideset ', STRA(:LR1),
     *          ' to ', STRB(:LR2),' ',STRING(:LENSTR(STRING)),
     *          ' tol. ', RSTR(4)(:LR4),
     *          ' max delta ', RSTR(5)(:LR5)
            else
              call numstr(3, 4, VECTOR(1,i), RSTR,  LR3)
              call numstr1(4, gap(i), RSTR(6), LR6)
              write (*, 40) SMTYP(:4), ' Sideset ', STRA(:LR1),
     *          ' to ', STRB(:LR2),' ',
     *          STRING(:LENSTR(STRING)),' ',
     *          RSTR(1)(:LR3), 'X ', RSTR(2)(:LR3), 'Y ',
     *          RSTR(3)(:LR3), 'Z ',
     *          ' tol. ', RSTR(4)(:LR4),
     *          ' max delta ', RSTR(5)(:LR5),
     *          ' gap ', RSTR(6)(:LR6)
            end if
          end do
        else
          write (*, 40) 'No sideset snapping or moving specified'
        end if
      ELSE IF (SHOTYP .EQ. 'WARP') THEN
        STRING = 'ERROR'
        IF (IWARP .EQ. 1) THEN
          STRING = 'Origin'
        ELSE IF (IWARP .EQ. -1) THEN
          STRING = 'X axis'
        ELSE IF (IWARP .EQ. -2) THEN
          STRING = 'Y axis'
        ELSE IF (IWARP .EQ. -3) THEN
          STRING = 'Z axis'
        END IF

        RSTR(2) = 'ERROR'
        IF (NRMWRP .EQ. 1) THEN
          RSTR(2) = 'X axis'
        ELSE IF (NRMWRP .EQ. 2) THEN
          RSTR(2) = 'Y axis'
        ELSE IF (NRMWRP .EQ. 3) THEN
          RSTR(2) = 'Z axis'
        END IF

        CALL NUMSTR1(4, WRPDIS, RSTR, LR)

        WRITE (*, 40) 'Warp mesh about the ', STRING(:LENSTR(STRING)),
     *    ', Reference Radius = ', RSTR(1)(:LR), ', Normal Vector = ',
     *    RSTR(2)(:LENSTR(RSTR(2)))
      ELSE
        CALL PRTERR ('CMDWARN', 'Invalid SHOW option')
      END IF

      RETURN
 40   FORMAT (1X, 20A)
 50   format (8x, 3(1pE12.5,4x))
      END
