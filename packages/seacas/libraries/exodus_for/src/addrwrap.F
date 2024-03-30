C Copyright(C) 1999-2020, 2022 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

        SUBROUTINE ADRDUM
        END

#ifdef Build64
        SUBROUTINE EXCLOS (IDEXO, IERR)
C> \sa ex_close()

C       CLOSE EXODUS II FILE

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXCLOS4 (IDEXO4, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXCN2S (IDEXO, NSESS, NDESS, IXEESS, IXNESS,
     &          LTEESS, LTNESS, LTSESS, IERR)
C> \sa ex_cvt_nodes_to_sides()

C       CONVERT SIDE SET NODES TO SIDES

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER NSESS (*)       ! (R)
        INTEGER NDESS (*)       ! (R)
        INTEGER IXEESS (*)      ! (R)
        INTEGER IXNESS (*)      ! (R)
        INTEGER LTEESS (*)      ! (R)
        INTEGER LTNESS (*)      ! (R)
        INTEGER LTSESS (*)      ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXCN2S4 (IDEXO4, NSESS, NDESS, IXEESS, IXNESS,
     &          LTEESS, LTNESS, LTSESS, IERR4)
        IERR = IERR4
        END

        INTEGER FUNCTION EXCRE (PATH, ICMODE, ICOMPWS, IOWS, IERR)
C> \sa ex_create()

C       CREATE EXODUS II FILE

        IMPLICIT NONE
        CHARACTER* (*) PATH     ! (R)
        INTEGER ICMODE          ! (R)
        INTEGER ICOMPWS         ! (RW)
        INTEGER IOWS            ! (RW)
        INTEGER IERR            ! (W)

        INTEGER*4 EXCRE4
        EXTERNAL EXCRE4
        INTEGER*4 IRET4
        INTEGER IRET
        INTEGER*4 ICMODE4       ! (R)
        INTEGER*4 ICOMPWS4      ! (RW)
        INTEGER*4 IOWS4         ! (RW)
        INTEGER*4 IERR4         ! (W)
        INTEGER*4 IFLAG4

C ... This is hardwired instead of getting from exodusII.inc due to compiler errors about
C     exopen/excre already having a type...
        integer EX_MAPS_INT64_API
        integer EX_BULK_INT64_API
        integer EX_IDS_INT64_API
        integer EX_INQ_INT64_API
        parameter (EX_MAPS_INT64_API=   8192) ! All maps (id, order, ...) store int64_t values
        parameter (EX_IDS_INT64_API =  16384) ! All entity ids
        parameter (EX_BULK_INT64_API=  32768) ! All integer bulk data
        parameter (EX_INQ_INT64_API =  65536) ! All integer inquire arguments

C ... All wrapper routines are passing INTEGER*8, but at this
C     time, only the maps are supported...
        ICMODE4 = ICMODE
        IFLAG4 = EX_MAPS_INT64_API+EX_IDS_INT64_API+EX_BULK_INT64_API+
     *    EX_INQ_INT64_API
        ICMODE4 = IOR(ICMODE4,  IFLAG4)

C ... This needs to be set explicitly since all wrapper routines
C     are passing REAL*8 and the exo_jack.c routines all take double*
        ICOMPWS4 = 8
        IOWS4 = IOWS

        IRET4 = EXCRE4 (PATH, ICMODE4, ICOMPWS4, IOWS4, IERR4)
        ICOMPWS = ICOMPWS4
        IOWS = IOWS4
        IERR = IERR4
        IRET = IRET4
        EXCRE = IRET
        END

        SUBROUTINE EXERR (MODNAM, MSG, ERRNUM)
C> \sa ex_error()

C       ERROR REPORTING

        CHARACTER* (*) MODNAM   ! (R)
        CHARACTER* (*) MSG      ! (R)
        INTEGER ERRNUM          ! (R)

        INTEGER*4 ERRNUM4       ! (R)

        ERRNUM4 = ERRNUM
        CALL EXERR4 (MODNAM, MSG, ERRNUM4)
        END

        SUBROUTINE EXGATM (IDEXO, TIME, IERR)
C> \sa ex_get_all_times()

C       READ ALL TIME VALUES

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        REAL TIME (*)           ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXGATM4 (IDEXO4, TIME, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGCNS (IDEXO, IDNPSS, NNNPS, NDNPS, IXNNPS,
     &          IXDNPS, LTNNPS, FACNPS, IERR)
C> \sa ex_get_concat_node_sets()

C       READ CONCATENATED NODE SETS

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDNPSS (*)      ! (W)
        INTEGER NNNPS (*)       ! (W)
        INTEGER NDNPS (*)       ! (W)
        INTEGER IXNNPS (*)      ! (W)
        INTEGER IXDNPS (*)      ! (W)
        INTEGER LTNNPS (*)      ! (W)
        REAL FACNPS (*)         ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO

        CALL EXGCNS4 (IDEXO4, IDNPSS, NNNPS, NDNPS, IXNNPS,
     &                IXDNPS, LTNNPS, FACNPS, IERR4)

        IERR = IERR4
        END

        SUBROUTINE EXGFRM (IDEXO, NFRAME, IDSCF, COORD, ITAGS, IERR)
C> \sa ex_get_coordinate_frames()

C       READ COORDINATE FRAMES

        implicit none
        include 'exodusII.inc'

        INTEGER IDEXO           ! (R)
        INTEGER NFRAME          ! (R)
        INTEGER IDSCF(*)        ! (W)
        INTEGER ITAGS(*)        ! (W)
        INTEGER IERR            ! (W)

        REAL COORD(*)           ! (W)

        INTEGER*4 IDEXO4        ! (R)
        integer*4 nframe4       ! (R)
        INTEGER*4 IERR4         ! (W)

        idexo4 = idexo

        NFRAME = EXINQI (IDEXO, EXNCF)
        NFRAME4 = NFRAME

        call exgfrm4(idexo4, nframe4, idscf, coord, itags, ierr4)
        ierr = ierr4

        END

        SUBROUTINE EXPFRM (IDEXO, NFRAME, IDSCF, COORD, ITAGS, IERR)
C> \sa ex_put_coordinate_frames()

C       DEFINE/WRITE COORDINATE FRAMES

        implicit none
        include 'exodusII.inc'

        INTEGER IDEXO           ! (R)
        INTEGER NFRAME          ! (R)
        INTEGER IDSCF(*)        ! (R)
        INTEGER ITAGS(*)        ! (R)
        INTEGER IERR            ! (W)

        REAL COORD(*)           ! (R)

        INTEGER*4 IDEXO4        ! (R)
        integer*4 nframe4       ! (R)
        INTEGER*4 IERR4         ! (W)

        idexo4 = idexo
        NFRAME4 = NFRAME

        call expfrm4(idexo4, nframe4, idscf, coord, itags, ierr4)
        ierr = ierr4

        END

        SUBROUTINE EXGCON (IDEXO, NAMECO, IERR)
C> \sa ex_get_coord_names()

C       READ COORDINATE NAMES

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        CHARACTER* (*) NAMECO (*) ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXGCON4 (IDEXO4, NAMECO, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGCOR (IDEXO, XN, YN, ZN, IERR)
C> \sa ex_get_coord()

C       READ NODAL COORDINATES

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        REAL XN (*)             ! (W)
        REAL YN (*)             ! (W)
        REAL ZN (*)             ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXGCOR4 (IDEXO4, XN, YN, ZN, IERR4)

        IERR = IERR4
        END

        SUBROUTINE EXGCSS (IDEXO, IDESSS, NSESS, NDESS, IXEESS,
     &          IXDESS, LTEESS, LTSESS, FACESS, IERR)
C> \sa ex_get_concat_side_sets()

C       READ CONCATENATED SIDE SETS

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDESSS (*)      ! (W)
        INTEGER NSESS (*)       ! (W)
        INTEGER NDESS (*)       ! (W)
        INTEGER IXEESS (*)      ! (W)
        INTEGER IXDESS (*)      ! (W)
        INTEGER LTEESS (*)      ! (W)
        INTEGER LTSESS (*)      ! (W)
        REAL FACESS (*)         ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO

        CALL EXGCSS4 (IDEXO4, IDESSS, NSESS, NDESS, IXEESS,
     &    IXDESS, LTEESS, LTSESS, FACESS, IERR4)

        IERR = IERR4
        END

        SUBROUTINE EXGCSSF (IDEXO, IDESSS, NSESS, NDESS, IXEESS,
     &          IXDESS, LTEESS, LTSESS, IERR)
C> \sa ex_get_concat_side_sets()

C>       READ CONCATENATED SIDE SETS, no distribution factors

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDESSS (*)      ! (W)
        INTEGER NSESS  (*)      ! (W)
        INTEGER NDESS  (*)      ! (W)
        INTEGER IXEESS (*)      ! (W)
        INTEGER IXDESS (*)      ! (W)
        INTEGER LTEESS (*)      ! (W)
        INTEGER LTSESS (*)      ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO

        CALL EXGCSSF4 (IDEXO4, IDESSS, NSESS, NDESS, IXEESS,
     &          IXDESS, LTEESS, LTSESS, IERR4)

        IERR = IERR4
        END

        SUBROUTINE EXGEAT (IDEXO, IDELB, ATRIB, IERR)
C> \sa ex_get_elem_attr()

C>       READ ELEMENT BLOCK ATTRIBUTES

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDELB           ! (R)
        REAL ATRIB (*)          ! (NUMATR,NUMELB) ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXGEAT4 (IDEXO4, IDELB, ATRIB, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGEAN (IDEXO, IDELB, NATTR, NAMES, IERR)
C> \sa ex_get_elem_attr_names()

C>       READ ELEMENT BLOCK ATTRIBUTES

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDELB           ! (R)
        INTEGER NATTR           ! (R)
        CHARACTER* (*) NAMES (*) ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 NATTR4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        NATTR4 = NATTR
        CALL EXGEAN4 (IDEXO4, IDELB, NATTR4, NAMES, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGEBI (IDEXO, IDELBS, IERR)
C> \sa ex_get_elem_blk_ids

C>       READ ELEMENT BLOCK IDS

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDELBS (*)      ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXGEBI4 (IDEXO4, IDELBS, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGELB (IDEXO, IDELB, NAMELB, NUMELB, NUMLNK,
     &          NUMATR, IERR)
C> \sa ex_get_elem_block()

C>       READ ELEMENT BLOCK PARAMETERS

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDELB           ! (R)
        CHARACTER* (*) NAMELB   ! (W)
        INTEGER NUMELB          ! (W)
        INTEGER NUMLNK          ! (W)
        INTEGER NUMATR          ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXGELB4 (IDEXO4, IDELB, NAMELB, NUMELB, NUMLNK,
     &          NUMATR, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGELC (IDEXO, IDELB, LINK, IERR)
C> \sa ex_get_elem_conn()

C>       READ ELEMENT BLOCK CONNECTIVITY

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDELB           ! (R)
        INTEGER LINK (*)        ! (NUMLNK,NUMELB) ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXGELC4 (IDEXO4, IDELB, LINK, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGECPP (IDEXO, ITYPE, IDELB, COUNTS, IERR)
C> \sa ex_get_entity_count_per_polyhedra()

C>       READ ENTITY_COUNT-PER-POLYHEDRA

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDELB           ! (R)
        INTEGER ITYPE           ! (R)
        INTEGER COUNTS (*)      ! (NUMELB) ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 ITYPE4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        ITYPE4 = ITYPE
        CALL EXGECPP4 (IDEXO4, ITYPE4, IDELB, COUNTS, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGENM (IDEXO, MAPEL, IERR)
C> \sa ex_get_id_map()

C>       READ ELEMENT NUMBER MAP

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER MAPEL (*)       ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXGENM4 (IDEXO4, MAPEL, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGEM (IDEXO, MAPID, MAPEL, IERR)
C> \sa ex_get_num_map()

C>       READ ELEMENT ORDER MAP

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER MAPID           ! (R)
        INTEGER MAPEL (*)       ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXGEM4 (IDEXO4, MAPID, MAPEL, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGPEM (IDEXO, MAPID, ISTART, ICOUNT, MAPEL, IERR)
C> \sa ex_get_partial_num_map()

C>       READ ELEMENT ORDER MAP

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER MAPID           ! (R)
        INTEGER ISTART          ! (R)
        INTEGER ICOUNT          ! (R)
        INTEGER MAPEL (*)       ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXGPEM4 (IDEXO4, MAPID, ISTART, ICOUNT, MAPEL, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGEV (IDEXO, ISTEP, IXELEV, IDELB, NUMELB,
     &          VALEV, IERR)
C> \sa ex_get_elem_var()

C>       READ ELEMENT VARIABLE VALUES AT A TIME STEP

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER ISTEP           ! (R)
        INTEGER IXELEV          ! (R)
        INTEGER IDELB           ! (R)
        INTEGER NUMELB          ! (R)
        REAL VALEV (*)          ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 ISTEP4        ! (R)
        INTEGER*4 IXELEV4       ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        ISTEP4 = ISTEP
        IXELEV4 = IXELEV
        CALL EXGEV4 (IDEXO4, ISTEP4, IXELEV4, IDELB, NUMELB,
     &          VALEV, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGEVT (IDEXO, IXELEV, IELNUM, ISTPB, ISTPE,
     &          VALEV, IERR)
C> \sa ex_get_elem_var_time()

C>       READ ELEMENT VARIABLE VALUES THROUGH TIME

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IXELEV          ! (R)
        INTEGER IELNUM          ! (R)
        INTEGER ISTPB           ! (R)
        INTEGER ISTPE           ! (R)
        REAL VALEV (*)          ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IXELEV4       ! (R)
        INTEGER*4 ISTPB4        ! (R)
        INTEGER*4 ISTPE4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4  = IDEXO
        IXELEV4 = IXELEV
        ISTPB4  = ISTPB
        ISTPE4  = ISTPE

        CALL EXGEVT4 (IDEXO4, IXELEV4, IELNUM, ISTPB4, ISTPE4,
     &          VALEV, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGNSV (IDEXO, ISTEP, IXV, ID, NUM, VAL, IERR)
C> \sa ex_get_nset_var()

C>       READ NODESET VARIABLE VALUES AT A TIME STEP

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER ISTEP           ! (R)
        INTEGER IXV             ! (R)
        INTEGER ID              ! (R)
        INTEGER NUM             ! (R)
        REAL VAL (*)            ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 ISTEP4        ! (R)
        INTEGER*4 IXV4  ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        ISTEP4 = ISTEP
        IXV4   = IXV
        CALL EXGNSV4 (IDEXO4, ISTEP4, IXV4, ID, NUM, VAL, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGSSV (IDEXO, ISTEP, IXV, ID, NUM, VAL, IERR)
C> \sa ex_get_sset_var()

C>       READ SIDESET VARIABLE VALUES AT A TIME STEP

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER ISTEP           ! (R)
        INTEGER IXV             ! (R)
        INTEGER ID              ! (R)
        INTEGER NUM             ! (R)
        REAL VAL (*)            ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 ISTEP4        ! (R)
        INTEGER*4 IXV4          ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        ISTEP4 = ISTEP
        IXV4   = IXV
        CALL EXGSSV4 (IDEXO4, ISTEP4, IXV4, ID, NUM, VAL, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGGV (IDEXO, ISTEP, NVARGL, VALGV, IERR)
C> \sa ex_get_glob_vars()

C>       READ GLOBAL VARIABLE VALUES AT A TIME STEP

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER ISTEP           ! (R)
        INTEGER NVARGL          ! (R)
        REAL VALGV (*)          ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 ISTEP4        ! (R)
        INTEGER*4 NVARGL4       ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        ISTEP4 = ISTEP
        NVARGL4 = NVARGL
        CALL EXGGV4 (IDEXO4, ISTEP4, NVARGL4, VALGV, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGGVT (IDEXO, IXGLOV, ISTPB, ISTPE, VALGV, IERR)
C> \sa ex_get_glob_var_time()

C>       READ GLOBAL VARIABLE VALUES THROUGH TIME

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IXGLOV          ! (R)
        INTEGER ISTPB           ! (R)
        INTEGER ISTPE           ! (R)
        REAL VALGV (*)          ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IXGLOV4       ! (R)
        INTEGER*4 ISTPB4        ! (R)
        INTEGER*4 ISTPE4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4  = IDEXO
        IXGLOV4 = IXGLOV
        ISTPB4  = ISTPB
        ISTPE4  = ISTPE
        CALL EXGGVT4 (IDEXO4, IXGLOV4, ISTPB4, ISTPE4, VALGV, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGINF (IDEXO, INFO, IERR)
C> \sa ex_get_info()

C>       READ INFORMATION RECORDS

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        CHARACTER* (*) INFO (*) ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXGINF4 (IDEXO4, INFO, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGINI (IDEXO, TITLE, NDIM, NUMNP, NUMEL, NELBLK,
     &          NUMNPS, NUMESS, IERR)
C> \sa ex_get_init()

C>       READ INITIALIZATION PARAMETERS

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        CHARACTER* (*) TITLE    ! (W)
        INTEGER NDIM            ! (W)
        INTEGER NUMNP           ! (W)
        INTEGER NUMEL           ! (W)
        INTEGER NELBLK          ! (W)
        INTEGER NUMNPS          ! (W)
        INTEGER NUMESS          ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXGINI4 (IDEXO4, TITLE, NDIM, NUMNP, NUMEL, NELBLK,
     &          NUMNPS, NUMESS, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGMP (IDEXO, NNMAP, NEMAP, IERR)
C> \sa ex_get_map_param()

C>       GET MAP PARAMETERS

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER NNMAP           ! (W)
        INTEGER NEMAP           ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 NNMAP4        ! (W)
        INTEGER*4 NEMAP4        ! (W)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXGMP4 (IDEXO4, NNMAP4, NEMAP4, IERR4)
        NNMAP = NNMAP4
        NEMAP = NEMAP4
        IERR = IERR4
        END

        SUBROUTINE EXPMP (IDEXO, NNMAP, NEMAP, IERR)
C> \sa ex_put_map_param()

C       PUT MAP PARAMETERS

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER NNMAP           ! (W)
        INTEGER NEMAP           ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 NNMAP4        ! (W)
        INTEGER*4 NEMAP4        ! (W)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        NNMAP4 = NNMAP
        NEMAP4 = NEMAP
        CALL EXPMP4 (IDEXO4, NNMAP4, NEMAP4, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGMAP (IDEXO, MAPEL, IERR)
C> \sa ex_get_map()

C>       READ ELEMENT ORDER MAP

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER MAPEL (*)       ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXGMAP4 (IDEXO4, MAPEL, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGNNM (IDEXO, MAPNOD, IERR)
C> \sa ex_get_id_map()

C>       READ NODE NUMBER MAP

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER MAPNOD (*)      ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXGNNM4 (IDEXO4, MAPNOD, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGNM (IDEXO, MAPID, MAPNOD, IERR)
C> \sa ex_get_num_map()

C>       READ NODE ORDER MAP

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER MAPID           ! (R)
        INTEGER MAPNOD (*)      ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXGNM4 (IDEXO4, MAPID, MAPNOD, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGNP (IDEXO, IDNPS, NNNPS, NDNPS, IERR)
C> \sa ex_get_node_set_param()

C>       READ NODE SET PARAMETERS

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDNPS           ! (R)
        INTEGER NNNPS           ! (W)
        INTEGER NDNPS           ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXGNP4 (IDEXO4, IDNPS, NNNPS, NDNPS, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGNS (IDEXO, IDNPS, LTNNPS, IERR)
C> \sa ex_get_node_set()

C>       READ NODE SET

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDNPS           ! (R)
        INTEGER LTNNPS (*)      ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXGNS4 (IDEXO4, IDNPS, LTNNPS, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGNSD (IDEXO, IDNPS, FACNPS, IERR)
C> \sa ex_get_node_set_dist_fact()

C>       READ NODE SET DISTRIBUTION FACTORS

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDNPS           ! (R)
        REAL FACNPS (*)         ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO

        CALL EXGNSD4 (IDEXO4, IDNPS, FACNPS, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGNSI (IDEXO, IDNPSS, IERR)
C> \sa ex_get_node_set_ids()

C>       READ NODE SET IDS

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDNPSS (*)      ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXGNSI4 (IDEXO4, IDNPSS, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGNV (IDEXO, ISTEP, IXNODV, NUMNP, VALNV, IERR)
C> \sa ex_get_nodal_var()

C>       READ NODAL VARIABLE VALUES AT A TIME STEP

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER ISTEP           ! (R)
        INTEGER IXNODV          ! (R)
        INTEGER NUMNP           ! (R)
        REAL VALNV (*)          ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 ISTEP4        ! (R)
        INTEGER*4 IXNODV4       ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4  = IDEXO
        ISTEP4  = ISTEP
        IXNODV4 = IXNODV
        CALL EXGNV4 (IDEXO4, ISTEP4, IXNODV4, NUMNP, VALNV, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGNVT (IDEXO, IXNODV, NODNUM, ISTPB, ISTPE,
     &          VALNV, IERR)
C> \sa ex_get_nodal_var_time()

C>       READ NODAL VARIABLE VALUES THROUGH TIME

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IXNODV          ! (R)
        INTEGER NODNUM          ! (R)
        INTEGER ISTPB           ! (R)
        INTEGER ISTPE           ! (R)
        REAL VALNV (*)          ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IXNODV4       ! (R)
        INTEGER*4 ISTPB4        ! (R)
        INTEGER*4 ISTPE4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4  = IDEXO
        IXNODV4 = IXNODV
        ISTPB4  = ISTPB
        ISTPE4  = ISTPE
        CALL EXGNVT4 (IDEXO4, IXNODV4, NODNUM, ISTPB4, ISTPE4,
     &          VALNV, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGP (IDEXO, ITYPE, ID, NAMEPR, IVAL, IERR)
C> \sa ex_get_prop()

C>       READ OBJECT PROPERTY

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER ITYPE           ! (R)
        INTEGER ID              ! (R)
        CHARACTER* (*) NAMEPR   ! (R)
        INTEGER IVAL            ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 ITYPE4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        ITYPE4 = ITYPE
        CALL EXGP4 (IDEXO4, ITYPE4, ID, NAMEPR, IVAL, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGPA (IDEXO, ITYPE, NAMEPR, IVAL, IERR)
C> \sa ex_get_prop_array()

C>       READ OBJECT PROPERTY ARRAY

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER ITYPE           ! (R)
        CHARACTER* (*) NAMEPR   ! (R)
        INTEGER IVAL (*)        ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 ITYPE4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        ITYPE4 = ITYPE
        CALL EXGPA4 (IDEXO4, ITYPE4, NAMEPR, IVAL, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGPN (IDEXO, ITYPE, NAMEPR, IERR)
C> \sa ex_get_prop_names()

C>       READ PROPERTY ARRAY NAMES

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER ITYPE           ! (R)
        CHARACTER* (*) NAMEPR (*) ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 ITYPE4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        ITYPE4 = ITYPE
        CALL EXGPN4 (IDEXO4, ITYPE4, NAMEPR, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGQA (IDEXO, QAREC, IERR)
C> \sa ex_get_qa()

C>       READ QA RECORDS

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        CHARACTER* (*) QAREC (4,*) ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXGQA4 (IDEXO4, QAREC, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGSP (IDEXO, IDESS, NSESS, NDESS, IERR)
C> \sa ex_get_side_set_param()

C>       READ SIDE SET PARAMETERS

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDESS           ! (R)
        INTEGER NSESS           ! (W)
        INTEGER NDESS           ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXGSP4 (IDEXO4, IDESS, NSESS, NDESS, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGSNL (IDEXO, IDESS, NNESS, IERR)
C> \sa ex_get_side_set_node_list_len()

C>       GET SIDE SET NODELIST LENGTH

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDESS           ! (R)
        INTEGER NNESS           ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXGSNL4 (IDEXO4, IDESS, NNESS, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGSS (IDEXO, IDESS, LTEESS, LTSESS, IERR)
C> \sa ex_get_side_set()

C>       READ SIDE SET

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDESS           ! (R)
        INTEGER LTEESS (*)      ! (W)
        INTEGER LTSESS (*)      ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXGSS4 (IDEXO4, IDESS, LTEESS, LTSESS, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGSSD (IDEXO, IDESS, FACESS, IERR)
C> \sa ex_get_side_set_dist_fact()

C>       READ SIDE SET DISTRIBUTION FACTORS

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDESS           ! (R)
        REAL FACESS (*)         ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXGSSD4 (IDEXO4, IDESS, FACESS, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGSSI (IDEXO, IDESSS, IERR)
C> \sa ex_get_side_set_ids()

C>       READ SIDE SET IDS

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDESSS (*)      ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXGSSI4 (IDEXO4, IDESSS, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGSSN (IDEXO, IDESS, INCNT, LTNESS, IERR)
C> \sa ex_get_side_set_node_list()

C>       READ SIDE SET NODE LIST

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDESS           ! (R)
        INTEGER INCNT (*)       ! (W)
        INTEGER LTNESS (*)      ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXGSSN4 (IDEXO4, IDESS, INCNT, LTNESS, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGSSC (IDEXO, IDESS, INCNT, IERR)
C> \sa ex_get_side_set_node_count()

C>       READ SIDE SET NODE LIST COUNT

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDESS           ! (R)
        INTEGER INCNT (*)       ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXGSSC4 (IDEXO4, IDESS, INCNT, IERR4)
        IERR = IERR4
        
        END

        SUBROUTINE EXGCSSC (IDEXO, INCNT, IERR)
C> \sa ex_get_concat_side_set_node_count()

C>       READ CONCATENATED SIDE SET NODE COUNT LIST
C>       (nodes/face for all faces)

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER INCNT (*)       ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXGCSSC4 (IDEXO4, INCNT, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGTIM (IDEXO, NSTEP, TIME, IERR)
C> \sa ex_get_time()

C>       READ TIME VALUE FOR A TIME STEP

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER NSTEP           ! (R)
        REAL TIME               ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 NSTEP4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        NSTEP4 = NSTEP
        CALL EXGTIM4 (IDEXO4, NSTEP4, TIME, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGVAN (IDEXO, VARTYP, NVAR, NAMES, IERR)
C> \sa ex_get_var_names()

C>       READ RESULTS VARIABLE NAMES

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        CHARACTER*1 VARTYP      ! (R)
        INTEGER NVAR            ! (R)
        CHARACTER* (*) NAMES (*) ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 NVAR4         ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        NVAR4 = NVAR
        CALL EXGVAN4 (IDEXO4, VARTYP, NVAR4, NAMES, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGVNM (IDEXO, VARTYP, NVAR, NAME, IERR)
C> \sa ex_get_var_name()

C>       READ RESULTS VARIABLE NAME

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        CHARACTER*1 VARTYP      ! (R)
        INTEGER NVAR            ! (R)
        CHARACTER* (*) NAME     ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 NVAR4         ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        NVAR4 = NVAR
        CALL EXGVNM4 (IDEXO4, VARTYP, NVAR4, NAME, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPVNM (IDEXO, VARTYP, NVAR, NAME, IERR)
C> \sa ex_put_var_name()

C>       WRITE RESULTS VARIABLE NAME

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        CHARACTER*1 VARTYP      ! (R)
        INTEGER NVAR            ! (R)
        CHARACTER* (*) NAME     ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 NVAR4         ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        NVAR4 = NVAR
        CALL EXPVNM4 (IDEXO4, VARTYP, NVAR4, NAME, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGVP (IDEXO, VARTYP, NVAR, IERR)
C> \sa ex_get_var_param()

C>       READ RESULTS VARIABLE PARAMETERS

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        CHARACTER*1 VARTYP      ! (R)
        INTEGER NVAR            ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 NVAR4         ! (W)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXGVP4 (IDEXO4, VARTYP, NVAR4, IERR4)
        NVAR = NVAR4
        IERR = IERR4
        END

        SUBROUTINE EXGVTT (IDEXO, NELBLK, NVAREL, ISEVOK, IERR)
C> \sa ex_get_elem_var_tab()

C>       READ ELEMENT VARIABLE TRUTH TABLE

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER NELBLK          ! (R)
        INTEGER NVAREL          ! (R)
        INTEGER ISEVOK(*)       ! (NVAREL,NELBLK) (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 NELBLK4       ! (R)
        INTEGER*4 NVAREL4       ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        NELBLK4 = NELBLK
        NVAREL4 = NVAREL
        CALL EXGVTT4 (IDEXO4, NELBLK4, NVAREL4, ISEVOK, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGNSTT (IDEXO, NBLK, NVAR, ISVOK, IERR)
C> \sa ex_get_nset_var_tab()

C>       READ NODESET VARIABLE TRUTH TABLE

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER NBLK            ! (R)
        INTEGER NVAR            ! (R)
        INTEGER ISVOK(*)        ! (NVAR,NBLK) (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 NBLK4         ! (R)
        INTEGER*4 NVAR4         ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        NBLK4  = NBLK
        NVAR4  = NVAR
        CALL EXGNSTT4 (IDEXO4, NBLK4, NVAR4, ISVOK, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGSSTT (IDEXO, NBLK, NVAR, ISVOK, IERR)
C> \sa ex_get_nset_var_tab()

C>       READ SIDESET VARIABLE TRUTH TABLE

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER NBLK            ! (R)
        INTEGER NVAR            ! (R)
        INTEGER ISVOK(*)        ! (NVAR,NBLK) (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 NBLK4         ! (R)
        INTEGER*4 NVAR4         ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        NBLK4  = NBLK
        NVAR4  = NVAR
        CALL EXGSSTT4 (IDEXO4, NBLK4, NVAR4, ISVOK, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXINQ (IDEXO, INFREQ, INTRET, RELRET, CHRRET, IERR)
C> \sa ex_inquire()

C>       INQUIRE EXODUS PARAMETERS

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER INFREQ          ! (R)
        INTEGER INTRET          ! (W)
        REAL RELRET             ! (W)
        CHARACTER* (*) CHRRET   ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 INFREQ4       ! (R)
        REAL*4 RELRET4          ! (W)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        INFREQ4 = INFREQ
        CALL EXINQ4 (IDEXO4, INFREQ4, INTRET, RELRET4, CHRRET, IERR4)
        RELRET = RELRET4
        IERR = IERR4
        END

        INTEGER FUNCTION EXINQI (IDEXO, INFREQ)
C> \sa ex_inquire_int()

C>       INQUIRE EXODUS PARAMETERS -- Integer only

        IMPLICIT NONE
        INTEGER IDEXO           ! (R)
        INTEGER INFREQ          ! (R)

        INTEGER*4 EXINQI4
        EXTERNAL EXINQI4

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 INFREQ4       ! (R)

        IDEXO4 = IDEXO
        INFREQ4 = INFREQ
        EXINQI = EXINQI4 (IDEXO4, INFREQ4)
        END

        INTEGER FUNCTION EXOPEN (PATH, IMODE, ICOMPWS, IOWS, VERS, IERR)
C> \sa ex_open()

C>       OPEN EXODUS II FILE

        IMPLICIT NONE
        CHARACTER* (*) PATH     ! (R)
        INTEGER IMODE           ! (R)
        INTEGER ICOMPWS         ! (RW)
        INTEGER IOWS            ! (RW)
        REAL VERS               ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 EXOPEN4
        EXTERNAL EXOPEN4
        INTEGER*4 IRET4
        INTEGER IRET
        INTEGER*4 IMODE4        ! (R)
        INTEGER*4 ICOMPWS4      ! (RW)
        INTEGER*4 IOWS4         ! (RW)
        REAL*4 VERS4            ! (W)
        INTEGER*4 IERR4         ! (W)
        INTEGER*4 IFLAG4

C ... This is hardwired instead of getting from exodusII.inc due to compiler errors about
C     exopen/excre already having a type...
        integer EX_MAPS_INT64_API
        integer EX_BULK_INT64_API
        integer EX_IDS_INT64_API
        integer EX_INQ_INT64_API
        parameter (EX_MAPS_INT64_API=   8192) ! All maps (id, order, ...) store int64_t values
        parameter (EX_IDS_INT64_API =  16384) ! All entity ids
        parameter (EX_BULK_INT64_API=  32768) ! All integer bulk data
        parameter (EX_INQ_INT64_API =  65536) ! All integer inquire arguments

        IMODE4 = IMODE

C ... This needs to be set explicitly since all wrapper routines
C     are passing REAL*8 and the exo_jack.c routines all take double*
        ICOMPWS4 = 8
        IOWS4 = IOWS

C ... All wrapper routines are passing INTEGER*8, but at this
C     time, only the maps are supported...
        IFLAG4 = EX_MAPS_INT64_API+EX_IDS_INT64_API+EX_BULK_INT64_API+
     *    EX_INQ_INT64_API
        IMODE4 = IOR(IMODE4, IFLAG4)

        IRET4 = EXOPEN4 (PATH, IMODE4, ICOMPWS4, IOWS4, VERS4, IERR4)
        ICOMPWS = ICOMPWS4
        IOWS = IOWS4
        VERS = VERS4
        IERR = IERR4
        IRET = IRET4
        EXOPEN = IRET
        END

        INTEGER FUNCTION EXLGMD (IDEXO)
C> \sa ex_large_model()
        IMPLICIT NONE
        INTEGER IDEXO           ! (R)
        INTEGER IRET

        INTEGER*4 EXLGMD4
        EXTERNAL  EXLGMD4

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IRET4

        IDEXO4 = IDEXO
        IRET4  = EXLGMD4(IDEXO4)
        IRET   = IRET4
        EXLGMD = IRET
        END

        SUBROUTINE EXOPTS (OPTVAL, IERR)
C> \sa ex_opts()

C>       SET ERROR REPORTING LEVEL

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER OPTVAL          ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 OPTVAL4       ! (R)
        INTEGER*4 IERR4         ! (W)

        OPTVAL4 = OPTVAL
        CALL EXOPTS4 (OPTVAL4, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXMXNM (IDEXO, LENGTH, IERR)
C> \sa ex_set_max_name_length()

C>       SET MAXIMUM NAME LENGTH

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER LENGTH          ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 LENGTH4       ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        LENGTH4 = LENGTH
        CALL EXMXNM4 (IDEXO4, LENGTH4, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXSETOPT (IDEXO, OPTION, VALUE, IERR)
C> \sa ex_set_opt()

C>       SET SPECIFIED OPTION

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER OPTION          ! (R)
        INTEGER VALUE           ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 OPTION4       ! (R)
        INTEGER*4 VALUE4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        OPTION4 = OPTION
        VALUE4 = VALUE
        CALL EXSETOPT4 (IDEXO4, OPTION4, VALUE4, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPCNS (IDEXO, IDNPSS, NNNPS, NDNPS, IXNNPS,
     &          IXDNPS, LTNNPS, FACNPS, IERR)
C> \sa ex_put_concat_node_sets()

C>       WRITE CONCATENATED NODE SETS

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDNPSS (*)      ! (R)
        INTEGER NNNPS (*)       ! (R)
        INTEGER NDNPS (*)       ! (R)
        INTEGER IXNNPS (*)      ! (R)
        INTEGER IXDNPS (*)      ! (R)
        INTEGER LTNNPS (*)      ! (R)
        REAL FACNPS (*)         ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO

        CALL EXPCNS4 (IDEXO4, IDNPSS, NNNPS, NDNPS, IXNNPS,
     &          IXDNPS, LTNNPS, FACNPS, IERR4)

        IERR = IERR4
        END

        SUBROUTINE EXPCON (IDEXO, NAMECO, IERR)
C> \sa ex_put_coord_names()

C>       WRITE COORDINATE NAMES

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        CHARACTER* (*) NAMECO (*) ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO

        CALL EXPCON4 (IDEXO4, NAMECO, IERR4)

        IERR = IERR4
        END

        SUBROUTINE EXPCOR (IDEXO, XN, YN, ZN, IERR)
C> \sa ex_put_coord()

C>       WRITE NODAL COORDINATES

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        REAL XN (*)             ! (R)
        REAL YN (*)             ! (R)
        REAL ZN (*)             ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO

        CALL EXPCOR4 (IDEXO4, XN, YN, ZN, IERR4)

        IERR = IERR4
        END

        SUBROUTINE EXPCSS (IDEXO, IDESSS, NSESS, NDESS, IXEESS,
     &          IXDESS, LTEESS, LTSESS, FACESS, IERR)
C> \sa ex_put_concat_side_sets()

C>       WRITE CONCATENATED SIDE SETS

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDESSS (*)      ! (R)
        INTEGER NSESS (*)       ! (R)
        INTEGER NDESS (*)       ! (R)
        INTEGER IXEESS (*)      ! (R)
        INTEGER IXDESS (*)      ! (R)
        INTEGER LTEESS (*)      ! (R)
        INTEGER LTSESS (*)      ! (R)
        REAL FACESS (*)         ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO

        CALL EXPCSS4 (IDEXO4, IDESSS, NSESS, NDESS, IXEESS,
     &          IXDESS, LTEESS, LTSESS, FACESS, IERR4)

        IERR = IERR4
        END

        SUBROUTINE EXPEAT (IDEXO, IDELB, ATRIB, IERR)
C> \sa ex_put_elem_attr()

C>       WRITE ELEMENT BLOCK ATTRIBUTES

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDELB           ! (R)
        REAL ATRIB(*)           ! (NUMATR,NUMELB) ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO

        CALL EXPEAT4 (IDEXO4, IDELB, ATRIB, IERR4)

        IERR = IERR4
        END

        SUBROUTINE EXPEAN (IDEXO, IDELB, NATTR, NAMES, IERR)
C> \sa ex_put_elem_attr_names()

C>       WRITE ELEMENT BLOCK ATTRIBUTE NAMES

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDELB           ! (R)
        INTEGER NATTR           ! (R)
        CHARACTER* (*) NAMES (*) ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 NATTR4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        NATTR4 = NATTR

        CALL EXPEAN4 (IDEXO4, IDELB, NATTR4, NAMES, IERR4)

        IERR = IERR4
        END

        SUBROUTINE EXPCLB (IDEXO, IDELB, NAMELB, NUMELB, NUMLNK,
     &          NUMATR, MKMAP, IERR)
C> \sa ex_put_concat_elem_block()

C>       WRITE CONCATENATED ELEMENT BLOCK PARAMETERS

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDELB(*)        ! (R)
        CHARACTER*(*) NAMELB(*) ! (R)
        INTEGER NUMELB(*)       ! (R)
        INTEGER NUMLNK(*)       ! (R)
        INTEGER NUMATR(*)       ! (R)
        INTEGER MKMAP           ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 MKMAP4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        MKMAP4  = MKMAP

        CALL EXPCLB4 (IDEXO4, IDELB, NAMELB, NUMELB, NUMLNK,
     &          NUMATR, MKMAP4, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPELB (IDEXO, IDELB, NAMELB, NUMELB, NUMLNK,
     &          NUMATR, IERR)
C> \sa ex_put_elem_block()

C>       WRITE ELEMENT BLOCK PARAMETERS

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDELB           ! (R)
        CHARACTER* (*) NAMELB   ! (R)
        INTEGER NUMELB          ! (R)
        INTEGER NUMLNK          ! (R)
        INTEGER NUMATR          ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO

        CALL EXPELB4 (IDEXO4, IDELB, NAMELB, NUMELB, NUMLNK,
     &          NUMATR, IERR4)

        IERR = IERR4
        END

        SUBROUTINE EXPELC (IDEXO, IDELB, LINK, IERR)
C> \sa ex_put_elem_conn()

C>       WRITE ELEMENT BLOCK CONNECTIVITY

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDELB           ! (R)
        INTEGER LINK (*)        ! (NUMLNK,NUMELB) ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO

        CALL EXPELC4 (IDEXO4, IDELB, LINK, IERR4)

        IERR = IERR4
        END

        SUBROUTINE EXPECPP (IDEXO, ITYPE, IDELB, COUNTS, IERR)
C> ex_put_entity_count_per_polyhedra()

C>       WRITE ELEMENT BLOCK ELEMENT COUNT PER POLYHEDRA

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER ITYPE           ! (R)
        INTEGER IDELB           ! (R)
        INTEGER COUNTS (*)      ! (NUMELB) ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 ITYPE4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        ITYPE4 = ITYPE

        CALL EXPECPP4 (IDEXO4, ITYPE4, IDELB, COUNTS, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPENM (IDEXO, MAPEL, IERR)
C> ex_put_id_map()

C>       WRITE ELEMENT NUMBER MAP

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER MAPEL (*)       ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXPENM4 (IDEXO4, MAPEL, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPEM (IDEXO, MAPID, MAPEL, IERR)
C> ex_put_processor_elem_maps()

C>       WRITE ELEMENT NUMBER MAP

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER MAPID           ! (R)
        INTEGER MAPEL (*)       ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXPEM4 (IDEXO4, MAPID, MAPEL, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPPEM (IDEXO, MAPID, ISTART, ICOUNT, MAPEL, IERR)
C> ex_put_partial_num_map()

C>       WRITE ELEMENT NUMBER MAP

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER MAPID           ! (R)
        INTEGER ISTART          ! (R)
        INTEGER ICOUNT          ! (R)
        INTEGER MAPEL (*)       ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO

        CALL EXPPEM4 (IDEXO4, MAPID, ISTART, ICOUNT, MAPEL, IERR4)

        IERR = IERR4
        END

        SUBROUTINE EXPEV (IDEXO, ISTEP, IXELEV, IDELB, NUMELB,
     &          VALEV, IERR)
C> ex_put_elem_var()

C>       WRITE ELEMENT VARIABLE VALUES AT A TIME STEP

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER ISTEP           ! (R)
        INTEGER IXELEV          ! (R)
        INTEGER IDELB           ! (R)
        INTEGER NUMELB          ! (R)
        REAL VALEV (*)          ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 ISTEP4        ! (R)
        INTEGER*4 IXELEV4       ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4  = IDEXO
        ISTEP4  = ISTEP
        IXELEV4 = IXELEV
        CALL EXPEV4 (IDEXO4, ISTEP4, IXELEV4, IDELB, NUMELB,
     &          VALEV, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPNSV (IDEXO, ISTEP, IXV, ID, NUM, VAL, IERR)
C> ex_put_nset_var()

C>       WRITE NODESET VARIABLE VALUES AT A TIME STEP

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER ISTEP           ! (R)
        INTEGER IXV             ! (R)
        INTEGER ID              ! (R)
        INTEGER NUM             ! (R)
        REAL VAL (*)            ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 ISTEP4        ! (R)
        INTEGER*4 IXV4  ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4  = IDEXO
        ISTEP4  = ISTEP
        IXV4 = IXV
        CALL EXPNSV4 (IDEXO4, ISTEP4, IXV4, ID, NUM, VAL, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPSSV (IDEXO, ISTEP, IXV, ID, NUM, VAL, IERR)
C> ex_put_sset_var()

C>       WRITE SIDESET VARIABLE VALUES AT A TIME STEP

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER ISTEP           ! (R)
        INTEGER IXV             ! (R)
        INTEGER ID              ! (R)
        INTEGER NUM             ! (R)
        REAL VAL (*)            ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 ISTEP4        ! (R)
        INTEGER*4 IXV4  ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4  = IDEXO
        ISTEP4  = ISTEP
        IXV4 = IXV
        CALL EXPSSV4 (IDEXO4, ISTEP4, IXV4, ID, NUM, VAL, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPGV (IDEXO, ISTEP, NVARGL, VALGV, IERR)
C> ex_put_glob_vars()

C>       WRITE GLOBAL VARIABLE VALUES AT A TIME STEP

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER ISTEP           ! (R)
        INTEGER NVARGL          ! (R)
        REAL VALGV (*)          ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 ISTEP4        ! (R)
        INTEGER*4 NVARGL4       ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        ISTEP4 = ISTEP
        NVARGL4 = NVARGL
        CALL EXPGV4 (IDEXO4, ISTEP4, NVARGL4, VALGV, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPINF (IDEXO, NINFO, INFO, IERR)
C> ex_put_info()

C>       WRITE INFORMATION RECORDS

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER NINFO           ! (R)
        CHARACTER* (*) INFO (*) ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 NINFO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        NINFO4 = NINFO
        CALL EXPINF4 (IDEXO4, NINFO4, INFO, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPINI (IDEXO, TITLE, NDIM, NUMNP, NUMEL, NELBLK,
     &          NUMNPS, NUMESS, IERR)
C> ex_put_init()

C>       WRITE INITIALIZATION PARAMETERS

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        CHARACTER* (*) TITLE    ! (R)
        INTEGER NDIM            ! (R)
        INTEGER NUMNP           ! (R)
        INTEGER NUMEL           ! (R)
        INTEGER NELBLK          ! (R)
        INTEGER NUMNPS          ! (R)
        INTEGER NUMESS          ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXPINI4 (IDEXO4, TITLE, NDIM, NUMNP, NUMEL, NELBLK,
     &          NUMNPS, NUMESS, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPMAP (IDEXO, MAPEL, IERR)
C> ex_put_map()

C>       WRITE ELEMENT ORDER MAP

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER MAPEL (*)       ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXPMAP4 (IDEXO4, MAPEL, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPNNM (IDEXO, MAPNOD, IERR)
C> ex_put_id_map()

C>       WRITE NODE NUMBER MAP

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER MAPNOD (*)      ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXPNNM4 (IDEXO4, MAPNOD, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPNM (IDEXO, MAPID, MAPNOD, IERR)
C> ex_put_num_map()

C>       WRITE NODE NUMBER MAP

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER MAPID           ! (R)
        INTEGER MAPNOD (*)      ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXPNM4 (IDEXO4, MAPID, MAPNOD, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPNP (IDEXO, IDNPS, NNNPS, NDNPS, IERR)
C> ex_put_node_set_param()

C>       WRITE NODE SET PARAMETERS

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDNPS           ! (R)
        INTEGER NNNPS           ! (R)
        INTEGER NDNPS           ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXPNP4 (IDEXO4, IDNPS, NNNPS, NDNPS, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPNS (IDEXO, IDNPS, LTNNPS, IERR)
C> ex_put_node_set()

C>       WRITE NODE SET

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDNPS           ! (R)
        INTEGER LTNNPS (*)      ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXPNS4 (IDEXO4, IDNPS, LTNNPS, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPNSD (IDEXO, IDNPS, FACNPS, IERR)
C> ex_put_node_set_dist_fact()

C>       WRITE NODE SET DISTRIBUTION FACTORS

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDNPS           ! (R)
        REAL FACNPS (*)         ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXPNSD4 (IDEXO4, IDNPS, FACNPS, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPNV (IDEXO, ISTEP, IXNODV, NUMNP, VALNV, IERR)
C> ex_put_nodal_var()

C>       WRITE NODAL VARIABLE VALUES AT A TIME STEP

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER ISTEP           ! (R)
        INTEGER IXNODV          ! (R)
        INTEGER NUMNP           ! (R)
        REAL VALNV (*)          ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 ISTEP4        ! (R)
        INTEGER*4 IXNODV4       ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4  = IDEXO
        ISTEP4  = ISTEP
        IXNODV4 = IXNODV
        CALL EXPNV4 (IDEXO4, ISTEP4, IXNODV4, NUMNP, VALNV, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPP (IDEXO, ITYPE, ID, NAMEPR, IVAL, IERR)
C> ex_put_prop()

C>       WRITE OBJECT PROPERTY

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER ITYPE           ! (R)
        INTEGER ID              ! (R)
        CHARACTER* (*) NAMEPR   ! (R)
        INTEGER IVAL            ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 ITYPE4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        ITYPE4 = ITYPE
        CALL EXPP4 (IDEXO4, ITYPE4, ID, NAMEPR, IVAL, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPPA (IDEXO, ITYPE, NAMEPR, IVAL, IERR)
C> ex_put_prop_array()

C>       WRITE OBJECT PROPERTY ARRAY

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER ITYPE           ! (R)
        CHARACTER* (*) NAMEPR   ! (R)
        INTEGER IVAL (*)        ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 ITYPE4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        ITYPE4 = ITYPE
        CALL EXPPA4 (IDEXO4, ITYPE4, NAMEPR, IVAL, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPPN (IDEXO, ITYPE, NPROPS, NAMEPR, IERR)
C> ex_put_prop_names()

C>       WRITE PROPERTY ARRAY NAMES

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER ITYPE           ! (R)
        INTEGER NPROPS          ! (R)
        CHARACTER* (*) NAMEPR (*) ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 ITYPE4        ! (R)
        INTEGER*4 NPROPS4       ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        ITYPE4 = ITYPE
        NPROPS4 = NPROPS
        CALL EXPPN4 (IDEXO4, ITYPE4, NPROPS4, NAMEPR, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPQA (IDEXO, NQAREC, QAREC, IERR)
C> ex_put_qa()

C>       WRITE QA RECORDS

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER NQAREC          ! (R)
        CHARACTER* (*) QAREC (4,*) ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 NQAREC4       ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        NQAREC4 = NQAREC
        CALL EXPQA4 (IDEXO4, NQAREC4, QAREC, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPSP (IDEXO, IDESS, NSESS, NDESS, IERR)
C> ex_put_side_set_param()

C>       WRITE SIDE SET PARAMETERS

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDESS           ! (R)
        INTEGER NSESS           ! (R)
        INTEGER NDESS           ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXPSP4 (IDEXO4, IDESS, NSESS, NDESS, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPSS (IDEXO, IDESS, LTEESS, LTSESS, IERR)
C> ex_put_side_set()

C>       WRITE SIDE SET

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDESS           ! (R)
        INTEGER LTEESS (*)      ! (R)
        INTEGER LTSESS (*)      ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXPSS4 (IDEXO4, IDESS, LTEESS, LTSESS, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPSSD (IDEXO, IDESS, FACESS, IERR)
C> ex_put_side_set_dist_fact()

C>       WRITE SIDE SET DISTRIBUTION FACTORS

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDESS           ! (R)
        REAL FACESS (*)         ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXPSSD4 (IDEXO4, IDESS, FACESS, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPTIM (IDEXO, NSTEP, TIME, IERR)
C> ex_put_time()

C>       WRITE TIME VALUE FOR A TIME STEP

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER NSTEP           ! (R)
        REAL TIME               ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 NSTEP4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        NSTEP4 = NSTEP
        CALL EXPTIM4 (IDEXO4, NSTEP4, TIME, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPVAN (IDEXO, VARTYP, NVAR, NAMES, IERR)
C> ex_put_var_names()

C>       WRITE RESULTS VARIABLE NAMES

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        CHARACTER*1 VARTYP      ! (R)
        INTEGER NVAR            ! (R)
        CHARACTER* (*) NAMES (*) ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 NVAR4         ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        NVAR4 = NVAR
        CALL EXPVAN4 (IDEXO4, VARTYP, NVAR4, NAMES, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPVP (IDEXO, VARTYP, NVAR, IERR)
C> ex_put_var_param()

C>       WRITE RESULTS VARIABLE PARAMETERS

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        CHARACTER*1 VARTYP      ! (R)
        INTEGER NVAR            ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 NVAR4         ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        NVAR4 = NVAR
        CALL EXPVP4 (IDEXO4, VARTYP, NVAR4, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPVTT (IDEXO, NELBLK, NVAREL, ISEVOK, IERR)
C> ex_put_elem_var_tab()

C>       WRITE ELEMENT VARIABLE TRUTH TABLE

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER NELBLK          ! (R)
        INTEGER NVAREL          ! (R)
        INTEGER ISEVOK (*)      ! (NVAREL,NELBLK) (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 NELBLK4       ! (R)
        INTEGER*4 NVAREL4       ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        NELBLK4 = NELBLK
        NVAREL4 = NVAREL
        CALL EXPVTT4 (IDEXO4, NELBLK4, NVAREL4, ISEVOK, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPNSTT (IDEXO, NBLK, NVAR, ISVOK, IERR)
C> ex_put_nset_var_tab()

C>       WRITE NODESET VARIABLE TRUTH TABLE

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER NBLK            ! (R)
        INTEGER NVAR            ! (R)
        INTEGER ISVOK (*)       ! (NVAR,NBLK) (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 NBLK4 ! (R)
        INTEGER*4 NVAR4 ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        NBLK4 = NBLK
        NVAR4 = NVAR
        CALL EXPNSTT4 (IDEXO4, NBLK4, NVAR4, ISVOK, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPSSTT (IDEXO, NBLK, NVAR, ISVOK, IERR)
C> ex_put_sset_var_tab()

C>       WRITE SIDESET VARIABLE TRUTH TABLE

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER NBLK            ! (R)
        INTEGER NVAR            ! (R)
        INTEGER ISVOK (*)       ! (NVAR,NBLK) (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 NBLK4 ! (R)
        INTEGER*4 NVAR4 ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        NBLK4 = NBLK
        NVAR4 = NVAR
        CALL EXPSSTT4 (IDEXO4, NBLK4, NVAR4, ISVOK, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGNAMS (IDEXO, ITYPE, NVAR, NAMES, IERR)
C> ex_get_names()

C>       READ OBJECT NAMES

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER NVAR            ! (R)
        INTEGER ITYPE           ! (R)
        CHARACTER* (*) NAMES (*) ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 ITYPE4        ! (R)
        INTEGER*4 NVAR4         ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        ITYPE4 = ITYPE
        NVAR4 = NVAR
        CALL EXGNAMS4 (IDEXO4, ITYPE4, NVAR4, NAMES, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGNAM (IDEXO, ITYPE, ID, NAME, IERR)
C> ex_get_name()

C>       READ OBJECT NAME

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER ITYPE           ! (R)
        INTEGER ID              ! (R)
        CHARACTER*(*) NAME      ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 ITYPE4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        ITYPE4 = ITYPE
        CALL EXGNAM4 (IDEXO4, ITYPE4, ID, NAME, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPNAMS (IDEXO, ITYPE, NVAR, NAMES, IERR)
C> ex_put_names()

C>       WRITE OBJECT NAMES

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER NVAR            ! (R)
        INTEGER ITYPE           ! (R)
        CHARACTER* (*) NAMES (*) ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 ITYPE4        ! (R)
        INTEGER*4 NVAR4         ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        ITYPE4 = ITYPE
        NVAR4 = NVAR
        CALL EXPNAMS4 (IDEXO4, ITYPE4, NVAR4, NAMES, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPNAM (IDEXO, ITYPE, ID, NAME, IERR)
C> ex_put_name()

C>       WRITE OBJECT NAME

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER ITYPE           ! (R)
        INTEGER ID              ! (R)
        CHARACTER* (*) NAME     ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 ITYPE4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        ITYPE4 = ITYPE
        CALL EXPNAM4 (IDEXO4, ITYPE4, ID, NAME, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXUPDA (IDEXO, IERR)
C> ex_update()

C>       UPDATE EXODUS II FILE

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXUPDA4 (IDEXO4, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXCOPY (IDEXO1, IDEXO2, IERR)
C> ex_copy()

C>       COPY EXODUS II FILE

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO1          ! (R)
        INTEGER IDEXO2          ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO14       ! (R)
        INTEGER*4 IDEXO24       ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO14 = IDEXO1
        IDEXO24 = IDEXO2
        CALL EXCOPY4 (IDEXO14, IDEXO24, IERR4)
        IERR = IERR4
        END

C-----------------------------------------------------------------------

C   Get initial information from nemesis file

      subroutine exgii (idexo, nproc, nproc_in_f, ftype, ierr)
C> ex_get_init_info()
      implicit none
      integer idexo
      integer nproc
      integer nproc_in_f
      character*(*) ftype
      integer ierr

      integer*4 idexo4
      integer*4 nproc4
      integer*4 nproc_in_f4
      integer*4 ierr4

      idexo4 = idexo
      call exgii4 (idexo4, nproc4, nproc_in_f4, ftype, ierr4)
      nproc = nproc4
      nproc_in_f = nproc_in_f4
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C   Write initial information from nemesis file

      subroutine expii (idexo, nproc, nproc_in_f, ftype, ierr)
C> ex_put_init_info()
      implicit none
      integer idexo
      integer nproc
      integer nproc_in_f
      character*(*) ftype
      integer ierr

      integer*4 idexo4
      integer*4 nproc4
      integer*4 nproc_in_f4
      integer*4 ierr4

      idexo4 = idexo
      nproc4 = nproc
      nproc_in_f4 = nproc_in_f
      call expii4 (idexo4, nproc4, nproc_in_f4, ftype, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read initial global information

      subroutine exgig (idexo, nnodes_g, nelems_g, nelem_blks_g,
     $  nnode_sets_g, nside_sets_g, ierr)
C> ex_get_init_global()
      implicit none
      integer idexo
      integer nnodes_g
      integer nelems_g
      integer nelem_blks_g
      integer nnode_sets_g
      integer nside_sets_g
      integer ierr

      integer*4 idexo4
      integer*4 ierr4

      idexo4 = idexo
      call exgig4 (idexo4, nnodes_g, nelems_g, nelem_blks_g,
     $  nnode_sets_g, nside_sets_g, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write initial global information

      subroutine expig (idexo, nnodes_g, nelems_g, nelem_blks_g,
     $  nnode_sets_g, nside_sets_g, ierr)
C> ex_put_init_global()
      implicit none
      integer idexo
      integer nnodes_g
      integer nelems_g
      integer nelem_blks_g
      integer nnode_sets_g
      integer nside_sets_g
      integer ierr

      integer*4 idexo4
      integer*4 ierr4

      idexo4 = idexo
      call expig4 (idexo4, nnodes_g, nelems_g, nelem_blks_g,
     $  nnode_sets_g, nside_sets_g, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read load balance parameters

      subroutine exglbp (idexo, nint_nodes, nbor_nodes, next_nodes,
     $  nint_elems, nbor_elems, nnode_cmaps, nelem_cmaps,
     $  processor, ierr)
C> ex_get_loadbal_param()
      implicit none
      integer idexo
      integer nint_nodes
      integer nbor_nodes
      integer next_nodes
      integer nint_elems
      integer nbor_elems
      integer nnode_cmaps
      integer nelem_cmaps
      integer processor
      integer ierr

      integer*4 idexo4
      integer*4 processor4
      integer*4 ierr4

      idexo4 = idexo
      processor4 = processor
      call exglbp4 (idexo4, nint_nodes, nbor_nodes, next_nodes,
     $  nint_elems, nbor_elems, nnode_cmaps, nelem_cmaps,
     $  processor4, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write load balance parameters

      subroutine explbp (idexo, nint_nodes, nbor_nodes, next_nodes,
     $  nint_elems, nbor_elems, nnode_cmaps, nelem_cmaps,
     $  processor, ierr)
C> ex_put_loadbal_param()
      implicit none
      integer idexo
      integer nint_nodes
      integer nbor_nodes
      integer next_nodes
      integer nint_elems
      integer nbor_elems
      integer nnode_cmaps
      integer nelem_cmaps
      integer processor
      integer ierr

      integer*4 idexo4
      integer*4 processor4
      integer*4 ierr4

      idexo4 = idexo
      processor4 = processor
      call explbp4 (idexo4, nint_nodes, nbor_nodes, next_nodes,
     $  nint_elems, nbor_elems, nnode_cmaps, nelem_cmaps,
     $  processor4, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write concatenated load balance parameters

      subroutine explbpc (idexo, nint_nodes, nbor_nodes, next_nodes,
     $  nint_elems, nbor_elems, nnode_cmaps, nelem_cmaps, ierr)
C> ex_put_loadbal_param_cc()
      implicit none
      integer idexo
      integer nint_nodes(*)
      integer nbor_nodes(*)
      integer next_nodes(*)
      integer nint_elems(*)
      integer nbor_elems(*)
      integer nnode_cmaps(*)
      integer nelem_cmaps(*)
      integer ierr

      integer*4 idexo4
      integer*4 ierr4

      idexo4 = idexo
      call explbpc4 (idexo4, nint_nodes, nbor_nodes, next_nodes,
     $  nint_elems, nbor_elems, nnode_cmaps, nelem_cmaps,
     $  ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read global node set parameters

      subroutine exgnspg (idexo, ns_ids_glob, ns_n_cnt_glob,
     $  ns_df_cnt_glob, ierr)
C> ex_get_ns_param_global()
      implicit none
      integer idexo
      integer ns_ids_glob(*)
      integer ns_n_cnt_glob(*)
      integer ns_df_cnt_glob(*)
      integer ierr

      integer*4 idexo4
      integer*4 ierr4

      idexo4 = idexo
      call exgnspg4 (idexo4, ns_ids_glob, ns_n_cnt_glob,
     $  ns_df_cnt_glob, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write global node set parameters

      subroutine expnspg (idexo, global_ids, global_n_cnts,
     $  global_df_cnts, ierr)
C> ex_put_ns_param_global()
      integer idexo
      integer global_ids(*)
      integer global_n_cnts(*)
      integer global_df_cnts(*)
      integer ierr

      integer*4 idexo4
      integer*4 ierr4

      idexo4 = idexo
      call expnspg4 (idexo4, global_ids, global_n_cnts,
     $  global_df_cnts, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read global side set parameters

      subroutine exgsspg (idexo, ss_ids_glob, ss_n_cnt_glob,
     $  ss_df_cnt_glob, ierr)
C> ex_get_ss_param_global()
      implicit none
      integer idexo
      integer ss_ids_glob(*)
      integer ss_n_cnt_glob(*)
      integer ss_df_cnt_glob(*)
      integer ierr

      integer*4 idexo4
      integer*4 ierr4

      idexo4 = idexo
      call exgsspg4 (idexo4, ss_ids_glob, ss_n_cnt_glob,
     $  ss_df_cnt_glob, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write global side set parameters

      subroutine expsspg (idexo, global_ids, global_el_cnts,
     $  global_df_cnts, ierr)
C> ex_put_ss_param_global()
      implicit none
      integer idexo
      integer global_ids(*)
      integer global_el_cnts(*)
      integer global_df_cnts(*)
      integer ierr

      integer*4 idexo4
      integer*4 ierr4

      idexo4 = idexo
      call expsspg4 (idexo4, global_ids, global_el_cnts,
     $  global_df_cnts, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read global element block information

      subroutine exgebig (idexo, el_blk_ids, el_blk_cnts, ierr)
C> ex_get_eb_info_global()
      implicit none
      integer idexo
      integer el_blk_ids(*)
      integer el_blk_cnts(*)
      integer ierr

      integer*4 idexo4
      integer*4 ierr4

      idexo4 = idexo
      call exgebig4 (idexo4, el_blk_ids, el_blk_cnts, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write global element block information

      subroutine expebig (idexo, el_blk_ids, el_blk_cnts, ierr)
C> ex_put_eb_info_global()
      implicit none
      integer idexo
      integer el_blk_ids(*)
      integer el_blk_cnts(*)
      integer ierr

      integer*4 idexo4
      integer*4 ierr4

      idexo4 = idexo
      call expebig4 (idexo4, el_blk_ids, el_blk_cnts, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read side set element list and side set side list

      subroutine exgnss (idexo, ss_id, start_side_num, num_sides,
     $  ss_elem_list, ss_side_list, ierr)
C> ex_get_n_side_set()
      implicit none
      integer idexo
      integer ss_id
      integer start_side_num
      integer num_sides
      integer ss_elem_list(num_sides)
      integer ss_side_list(num_sides)
      integer ierr

      integer*4 idexo4
      integer*4 ierr4

      idexo4 = idexo
      call exgnss4 (idexo4, ss_id, start_side_num, num_sides,
     $  ss_elem_list, ss_side_list, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write side set element list and side set side list

      subroutine expnss (idexo, ss_id, start_side_num, num_sides,
     $  ss_elem_list, ss_side_list, ierr)
C> ex_put_n_side_set()
      implicit none
      integer idexo
      integer ss_id
      integer start_side_num
      integer num_sides
      integer ss_elem_list(num_sides)
      integer ss_side_list(num_sides)
      integer ierr

      integer*4 idexo4
      integer*4 ierr4

      idexo4 = idexo
      call expnss4 (idexo4, ss_id, start_side_num, num_sides,
     $  ss_elem_list, ss_side_list, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read side set distribution factor

      subroutine exgnssd (idexo, ss_id, start_num, num_df_to_get,
     $  ss_df, ierr)
C> ex_get_n_side_set_df()
      implicit none
      integer idexo
      integer ss_id
      integer start_num
      integer num_df_to_get
      real ss_df(num_df_to_get)
      integer ierr

      integer*4 idexo4
      integer*4 ierr4

      idexo4 = idexo
      call exgnssd4 (idexo4, ss_id, start_num, num_df_to_get,
     $  ss_df, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write side set distribution factor

      subroutine expnssd (idexo, ss_id, start_num, num_df_to_get,
     $  ss_df, ierr)
C> ex_put_n_side_set_df()
      implicit none
      integer idexo
      integer ss_id
      integer start_num
      integer num_df_to_get
      real ss_df(num_df_to_get)
      integer ierr

      integer*4 idexo4
      integer*4 ierr4

      idexo4 = idexo
      call expnssd4 (idexo4, ss_id, start_num, num_df_to_get,
     $  ss_df, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read node set list for a single node set

      subroutine exgnns (idexo, ns_id, start_node_num, num_node,
     $  ns_node_list, ierr)
C> ex_get_n_node_set()
      implicit none
      integer idexo
      integer ns_id
      integer start_node_num
      integer num_node
      integer ns_node_list(num_node)
      integer ierr

      integer*4 idexo4
      integer*4 ierr4

      idexo4 = idexo
      call exgnns4 (idexo4, ns_id, start_node_num, num_node,
     $  ns_node_list, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write node set list for a single node set

      subroutine expnns (idexo, ns_id, start_node_num, num_node,
     $  ns_node_list, ierr)
C> ex_put_n_node_set()
      implicit none
      integer idexo
      integer ns_id
      integer start_node_num
      integer num_node
      integer ns_node_list(num_node)
      integer ierr

      integer*4 idexo4
      integer*4 ierr4

      idexo4 = idexo
      call expnns4 (idexo4, ns_id, start_node_num, num_node,
     $  ns_node_list, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read node set distribution factor

      subroutine exgnnsd (idexo, ns_id, start_num, num_df_to_get,
     $  ns_df, ierr)
C> ex_get_n_node_set_df()
      implicit none
      integer idexo
      integer ns_id
      integer start_num
      integer num_df_to_get
      real ns_df(num_df_to_get)
      integer ierr

      integer*4 idexo4
      integer*4 ierr4

      idexo4 = idexo
      call exgnnsd4 (idexo4, ns_id, start_num, num_df_to_get,
     $  ns_df, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write node set distribution factor

      subroutine expnnsd (idexo, ns_id, start_num, num_df_to_get,
     $  ns_df, ierr)
C> ex_put_n_node_set_df()
      implicit none
      integer idexo
      integer ns_id
      integer start_num
      integer num_df_to_get
      real ns_df(num_df_to_get)
      integer ierr

      integer*4 idexo4
      integer*4 ierr4

      idexo4 = idexo
      call expnnsd4 (idexo4, ns_id, start_num, num_df_to_get,
     $  ns_df, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read coordinates of the nodes

      subroutine exgncor (idexo, start_node_num, num_nodes, x_coor,
     $  y_coor, z_coor, ierr)
C> ex_get_n_coord()
      implicit none
      integer idexo
      integer start_node_num
      integer num_nodes
      real x_coor(num_nodes)
      real y_coor(num_nodes)
      real z_coor(num_nodes)
      integer ierr

      integer*4 idexo4
      integer*4 ierr4

      idexo4 = idexo
      call exgncor4 (idexo4, start_node_num, num_nodes, x_coor,
     $  y_coor, z_coor, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write coordinates of the nodes

      subroutine expncor (idexo, start_node_num, num_nodes, x_coor,
     $  y_coor, z_coor, ierr)
C> ex_put_n_coord()
      implicit none
      integer idexo
      integer start_node_num
      integer num_nodes
      real x_coor(num_nodes)
      real y_coor(num_nodes)
      real z_coor(num_nodes)
      integer ierr

      integer*4 idexo4
      integer*4 ierr4

      idexo4 = idexo
      call expncor4 (idexo4, start_node_num, num_nodes, x_coor,
     $  y_coor, z_coor, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read an element blocks connectivity list

      subroutine exgnec (idexo, elem_blk_id, start_elem_num,
     $  num_elems, connect, ierr)
C> ex_get_n_elem_conn()
      implicit none
      integer idexo
      integer elem_blk_id
      integer start_elem_num
      integer num_elems
      integer connect(num_elems)
      integer ierr

      integer*4 idexo4
      integer*4 ierr4

      idexo4 = idexo
      call exgnec4 (idexo4, elem_blk_id, start_elem_num,
     $  num_elems, connect, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write an element blocks connectivity list

      subroutine expnec (idexo, elem_blk_id, start_elem_num,
     $  num_elems, connect, ierr)
C> ex_put_n_elem_conn()
      implicit none
      integer idexo
      integer elem_blk_id
      integer start_elem_num
      integer num_elems
      integer connect(num_elems)
      integer ierr

      integer*4 idexo4
      integer*4 ierr4

      idexo4 = idexo
      call expnec4 (idexo4, elem_blk_id, start_elem_num,
     $  num_elems, connect, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read an element blocks attributes

      subroutine exgneat (idexo, elem_blk_id, start_elem_num,
     $  num_elems, attrib, ierr)
C> ex_get_n_elem_attr()
      implicit none
      integer idexo
      integer elem_blk_id
      integer start_elem_num
      integer num_elems
      real attrib(num_elems)
      integer ierr

      integer*4 idexo4
      integer*4 ierr4

      idexo4 = idexo
      call exgneat4 (idexo4, elem_blk_id, start_elem_num,
     $  num_elems, attrib, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write an element blocks attributes

      subroutine expneat (idexo, elem_blk_id, start_elem_num,
     $  num_elems, attrib, ierr)
C> ex_put_n_elem_attr()
      implicit none
      integer idexo
      integer elem_blk_id
      integer start_elem_num
      integer num_elems
      real attrib(num_elems)
      integer ierr

      integer*4 idexo4
      integer*4 ierr4

      idexo4 = idexo
      call expneat4 (idexo4, elem_blk_id, start_elem_num,
     $  num_elems, attrib, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------
C  Read the element type for a specific element block

      subroutine exgelt (idexo, elem_blk_id, elem_type, ierr)
C> ex_get_elem_type()
      implicit none
      integer idexo
      integer elem_blk_id
      character*(*) elem_type
      integer ierr

      integer*4 idexo4
      integer*4 ierr4

      idexo4 = idexo
      call exgelt4 (idexo4, elem_blk_id, elem_type, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read a variable for an element block

      subroutine exgnev (idexo, time_step, elem_var_index,
     $  elem_blk_id, num_elem_this_blk, start_elem_num,
     $  num_elem, elem_var_vals, ierr)
C> ex_get_n_var()
      implicit none
      integer idexo
      integer time_step
      integer elem_var_index
      integer elem_blk_id
      integer num_elem_this_blk
      integer start_elem_num
      integer num_elem
      real elem_var_vals(num_elem)
      integer ierr

      integer*4 idexo4
      integer*4 time_step4
      integer*4 elem_var_index4
      integer*4 ierr4

      idexo4 = idexo
      time_step4 = time_step
      elem_var_index4 = elem_var_index
      call exgnev4 (idexo4, time_step4, elem_var_index4,
     $  elem_blk_id, num_elem_this_blk, start_elem_num,
     $  num_elem, elem_var_vals, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write a variable slab for an element block

      subroutine expevs (idexo, time_step, elem_var_index,
     $  elem_blk_id, start_pos, num_vals, elem_var_vals,
     $  ierr)
C> ex_put_elem_var_slab()
      implicit none
      integer idexo
      integer time_step
      integer elem_var_index
      integer elem_blk_id
      integer start_pos
      integer num_vals
      real elem_var_vals(num_vals)
      integer ierr

      integer*4 idexo4
      integer*4 time_step4
      integer*4 elem_var_index4
      integer*4 ierr4

      idexo4 = idexo
      time_step4 = time_step
      elem_var_index4 = elem_var_index
      call expevs4 (idexo4, time_step4, elem_var_index4,
     $  elem_blk_id, start_pos, num_vals, elem_var_vals,
     $  ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read the values of a single nodal variable for a single time step

      subroutine exgnnv (idexo, time_step, nodal_var_index,
     $  start_node_num, num_nodes, nodal_vars, ierr)
C> ex_get_n_var()
      implicit none
      integer idexo
      integer time_step
      integer nodal_var_index
      integer start_node_num
      integer num_nodes
      real nodal_vars(num_nodes)
      integer ierr

      integer*4 idexo4
      integer*4 time_step4
      integer*4 nodal_var_index4
      integer*4 ierr4

      idexo4 = idexo
      time_step4 = time_step
      nodal_var_index4 = nodal_var_index
      call exgnnv4 (idexo4, time_step4, nodal_var_index4,
     $  start_node_num, num_nodes, nodal_vars, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write nodal variable slab

      subroutine expnvs (idexo, time_step, nodal_var_index,
     $  start_pos, num_vals, nodal_var_vals, ierr)
C> ex_put_nodal_var_slab()
      implicit none
      integer idexo
      integer time_step
      integer nodal_var_index
      integer start_pos
      integer num_vals
      real nodal_var_vals(num_vals)
      integer ierr

      integer*4 idexo4
      integer*4 time_step4
      integer*4 nodal_var_index4
      integer*4 ierr4

      idexo4 = idexo
      time_step4 = time_step
      nodal_var_index4 = nodal_var_index
      call expnvs4 (idexo4, time_step4, nodal_var_index4,
     $  start_pos, num_vals, nodal_var_vals, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read the element numbering map

      subroutine exgnenm (idexo, starte, num_ent, elem_map, ierr)
C> ex_get_n_elem_num_map()
      implicit none
      integer idexo
      integer starte
      integer num_ent
      integer elem_map(num_ent)
      integer ierr

      integer*4 idexo4
      integer*4 ierr4

      idexo4 = idexo
      call exgnenm4 (idexo4, starte, num_ent, elem_map, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write the element numbering map

      subroutine expnenm (idexo, starte, num_ent, elem_map, ierr)
C> ex_put_partial_id_map()
      implicit none
      integer idexo
      integer starte
      integer num_ent
      integer elem_map(num_ent)
      integer ierr

      integer*4 idexo4
      integer*4 ierr4

      idexo4 = idexo
      call expnenm4 (idexo4, starte, num_ent, elem_map, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read the node numbering map

      subroutine exgnnnm (idexo, startn, num_ent, node_map, ierr)
C> ex_get_n_node_num_map()
      implicit none
      integer idexo
      integer startn
      integer num_ent
      integer node_map(num_ent)
      integer ierr

      integer*4 idexo4
      integer*4 ierr4

      idexo4 = idexo
      call exgnnnm4 (idexo4, startn, num_ent, node_map, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write the node numbering map

      subroutine expnnnm (idexo, startn, num_ent, node_map, ierr)
C> ex_put_partial_id_map()
      implicit none
      integer idexo
      integer startn
      integer num_ent
      integer node_map(num_ent)
      integer ierr

      integer*4 idexo4
      integer*4 ierr4

      idexo4 = idexo
      call expnnnm4 (idexo4, startn, num_ent, node_map, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read the node map for a processor

      subroutine exgnmp (idexo, node_mapi, node_mapb, node_mape,
     $  processor, ierr)
C> ex_get_processor_node_maps()
      implicit none
      integer idexo
      integer node_mapi(*)
      integer node_mapb(*)
      integer node_mape(*)
      integer processor
      integer ierr

      integer*4 idexo4
      integer*4 processor4
      integer*4 ierr4

      idexo4 = idexo
      processor4 = processor
      call exgnmp4 (idexo4, node_mapi, node_mapb, node_mape,
     $  processor4, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write a node map for a processor

      subroutine expnmp (idexo, node_mapi, node_mapb, node_mape,
     $  processor, ierr)
C> ex_put_processor_node_maps()
      implicit none
      integer idexo
      integer node_mapi(*)
      integer node_mapb(*)
      integer node_mape(*)
      integer processor
      integer ierr

      integer*4 idexo4
      integer*4 processor4
      integer*4 ierr4

      idexo4 = idexo
      processor4 = processor
      call expnmp4 (idexo4, node_mapi, node_mapb, node_mape,
     $  processor4, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read the element map for a processor

      subroutine exgemp (idexo, elem_mapi, elem_mapb, processor, ierr)
C> ex_get_processor_elem_maps()
      implicit none
      integer idexo
      integer elem_mapi(*)
      integer elem_mapb(*)
      integer processor
      integer ierr

      integer*4 idexo4
      integer*4 processor4
      integer*4 ierr4

      idexo4 = idexo
      processor4 = processor
      call exgemp4 (idexo4, elem_mapi, elem_mapb, processor4, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write the element map for a processor

      subroutine expemp (idexo, elem_mapi, elem_mapb, processor, ierr)
C> ex_put_processor_elem_maps()
      implicit none
      integer idexo
      integer elem_mapi(*)
      integer elem_mapb(*)
      integer processor
      integer ierr

      integer*4 idexo4
      integer*4 processor4
      integer*4 ierr4

      idexo4 = idexo
      processor4 = processor
      call expemp4 (idexo4, elem_mapi, elem_mapb, processor4, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read the communications map parameters for a single processor

      subroutine exgcmp (idexo, ncmap_ids, ncmap_node_cnts,
     $  ecmap_ids, ecmap_elem_cnts, processor, ierr)
C> ex_get_cmap_params()
      implicit none
      integer idexo
      integer ncmap_ids(*)
      integer ncmap_node_cnts(*)
      integer ecmap_ids(*)
      integer ecmap_elem_cnts(*)
      integer processor
      integer ierr

      integer*4 idexo4
      integer*4 processor4
      integer*4 ierr4

      idexo4 = idexo
      processor4 = processor
      call exgcmp4 (idexo4, ncmap_ids, ncmap_node_cnts,
     $  ecmap_ids, ecmap_elem_cnts, processor4, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write the communications map parameters for a single processor

      subroutine expcmp (idexo, nmap_ids, nmap_node_cnts, emap_ids,
     $  emap_elem_cnts, processor, ierr)
C> ex_put_cmap_params()
      implicit none
      integer idexo
      integer nmap_ids(*)
      integer nmap_node_cnts(*)
      integer emap_ids(*)
      integer emap_elem_cnts(*)
      integer processor
      integer ierr

      integer*4 idexo4
      integer*4 processor4
      integer*4 ierr4

      idexo4 = idexo
      processor4 = processor
      call expcmp4 (idexo4, nmap_ids, nmap_node_cnts, emap_ids,
     $  emap_elem_cnts, processor4, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write the communications map parameters for all processors

      subroutine expcmpc (idexo, nmap_ids, nmap_node_cnts,
     $  nproc_ptrs, emap_ids, emap_elem_cnts, eproc_ptrs, ierr)
C> ex_put_cmap_params()
      implicit none
      integer idexo
      integer nmap_ids(*)
      integer nmap_node_cnts(*)
      integer nproc_ptrs(*)
      integer emap_ids(*)
      integer emap_elem_cnts(*)
      integer eproc_ptrs(*)
      integer ierr

      integer*4 idexo4
      integer*4 ierr4

      idexo4 = idexo
      call expcmpc4 (idexo4, nmap_ids, nmap_node_cnts,
     $  nproc_ptrs, emap_ids, emap_elem_cnts, eproc_ptrs, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read the nodal communications map for a single processor

      subroutine exgncm (idexo, map_id, node_ids, proc_ids,
     $  processor, ierr)
C> ex_get_node_cmap()
      implicit none
      integer idexo
      integer map_id
      integer node_ids(*)
      integer proc_ids(*)
      integer processor
      integer ierr

      integer*4 idexo4
      integer*4 processor4
      integer*4 ierr4

      idexo4 = idexo
      processor4 = processor
      call exgncm4 (idexo4, map_id, node_ids, proc_ids,
     $  processor4, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write the nodal communications map for a single processor

      subroutine expncm (idexo, map_id, node_ids, proc_ids,
     $  processor, ierr)
C> ex_put_node_cmap()
      implicit none
      integer idexo
      integer map_id
      integer node_ids(*)
      integer proc_ids(*)
      integer processor
      integer ierr

      integer*4 idexo4
      integer*4 processor4
      integer*4 ierr4

      idexo4 = idexo
      processor4 = processor
      call expncm4 (idexo4, map_id, node_ids, proc_ids,
     $  processor4, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read the elemental communications map for a single processor

      subroutine exgecm (idexo, map_id, elem_ids, side_ids,
     $  proc_ids, processor, ierr)
C> ex_get_elem_cmap()
      implicit none
      integer idexo
      integer map_id
      integer elem_ids(*)
      integer side_ids(*)
      integer proc_ids(*)
      integer processor
      integer ierr

      integer*4 idexo4
      integer*4 processor4
      integer*4 ierr4

      idexo4 = idexo
      processor4 = processor
      call exgecm4 (idexo4, map_id, elem_ids, side_ids,
     $  proc_ids, processor4, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write the elemental communications map for a single processor

      subroutine expecm (idexo, map_id, elem_ids, side_ids,
     $  proc_ids, processor, ierr)
C> ex_put_elem_cmap()
      implicit none
      integer idexo
      integer map_id
      integer elem_ids(*)
      integer side_ids(*)
      integer proc_ids(*)
      integer processor
      integer ierr

      integer*4 idexo4
      integer*4 processor4
      integer*4 ierr4

      idexo4 = idexo
      processor4 = processor
      call expecm4 (idexo4, map_id, elem_ids, side_ids,
     $  proc_ids, processor4, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C> \sa ex_get_partial_coord_component()
       SUBROUTINE EXGPCC(exoid, start_node_num, num_nodes, component,
     & coor, ierr)
C> writes the coordinates of some of the nodes in the model for the specified component

       IMPLICIT NONE
       INCLUDE 'exodusII.inc'
       INTEGER exoid           ! (R)
       INTEGER start_node_num  ! (R)
       INTEGER num_nodes       ! (R)
       INTEGER component       ! (R)
       REAL coor(*)            ! (R)
       INTEGER ierr            ! (W)

       INTEGER*4 exoid4        ! (R)
       INTEGER*4 component4    ! (R)
       INTEGER*4 ierr4         ! (W)

       exoid4        = exoid
       component4    = component
       ierr4         = ierr

       CALL EXGPCC4(exoid4, start_node_num, num_nodes, component4,
     & coor, ierr4)
       IERR = IERR4
       END

C-----------------------------------------------------------------------

C> \sa ex_put_partial_coord_component()
       SUBROUTINE EXPPCC(exoid, start_node_num, num_nodes, component,
     & coor, ierr)
C> writes the coordinates of some of the nodes in the model for the specified component

       IMPLICIT NONE
       INCLUDE 'exodusII.inc'
       INTEGER exoid           ! (R)
       INTEGER start_node_num  ! (R)
       INTEGER num_nodes       ! (R)
       INTEGER component       ! (R)
       REAL coor(*)            ! (R)
       INTEGER ierr            ! (W)

       INTEGER*4 exoid4        ! (R)
       INTEGER*4 component4    ! (R)
       INTEGER*4 ierr4         ! (W)

       exoid4        = exoid
       component4    = component
       ierr4         = ierr

       CALL EXPPCC4(exoid4, start_node_num, num_nodes, component4,
     & coor, ierr4)
       IERR = IERR4
       END

C-----------------------------------------------------------------------

C> \sa ex_get_partial_var()
       SUBROUTINE EXGPV(idexo, time_step, var_type, var_index,
     &                  obj_id, start_index, num_entities,
     &                  var_vals, ierr)
C> writes the values of a single variable for a partial block at one time

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER idexo           ! (R)
        INTEGER time_step       ! (R)
        INTEGER var_type        ! (R)
        INTEGER var_index       ! (R)
        INTEGER obj_id          ! (R)
        INTEGER start_index     ! (R)
        INTEGER num_entities    ! (R)
        REAL    var_vals(*)     ! (W)
        INTEGER ierr            ! (W)

        INTEGER*4 idexo4        ! (R)
        INTEGER*4 time_step4    ! (R)
        INTEGER*4 var_type4     ! (R)
        INTEGER*4 var_index4    ! (R)
        INTEGER*4 ierr4         ! (W)

                idexo4        = idexo
                time_step4    = time_step
                var_type4     = var_type
                var_index4    = var_index
                ierr4         = ierr

        CALL EXGPV4(idexo4, time_step4, var_type4, var_index4,
     &                  obj_id, start_index, num_entities,
     &                  var_vals, ierr4)
        IERR = IERR4
        END

C-----------------------------------------------------------------------

C> \sa ex_put_partial_var()
       SUBROUTINE EXPPV(idexo, time_step, var_type, var_index,
     &                  obj_id, start_index, num_entities,
     &                  var_vals, ierr)
C> writes the values of a single variable for a partial block at one time

        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER idexo           ! (R)
        INTEGER time_step       ! (R)
        INTEGER var_type        ! (R)
        INTEGER var_index       ! (R)
        INTEGER obj_id          ! (R)
        INTEGER start_index     ! (R)
        INTEGER num_entities    ! (R)
        REAL    var_vals(*)     ! (W)
        INTEGER ierr            ! (W)

        INTEGER*4 idexo4        ! (R)
        INTEGER*4 time_step4    ! (R)
        INTEGER*4 var_type4     ! (R)
        INTEGER*4 var_index4    ! (R)
        INTEGER*4 ierr4         ! (W)

                idexo4        = idexo
                time_step4    = time_step
                var_type4     = var_type
                var_index4    = var_index
                ierr4         = ierr

        CALL EXPPV4(idexo4, time_step4, var_type4, var_index4,
     &                  obj_id, start_index, num_entities,
     &                  var_vals, ierr4)
        IERR = IERR4
        END

C-----------------------------------------------------------------------

#endif
