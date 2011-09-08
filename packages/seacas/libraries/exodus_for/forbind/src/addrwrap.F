C Copyright (c) 2005 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
C retains certain rights in this software.
C
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.
C
C     * Neither the name of Sandia Corporation nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C

        SUBROUTINE ADRDUM
        END

#ifdef Build64
        SUBROUTINE EXCLOS (IDEXO, IERR)
C
C       CLOSE EXODUS II FILE
C
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
C
C       CONVERT SIDE SET NODES TO SIDES
C
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
        INTEGER*4 NSESS4 (1)    ! (R)
        INTEGER*4 NDESS4 (1)    ! (R)
        INTEGER*4 IXEESS4 (1)   ! (R)
        INTEGER*4 IXNESS4 (1)   ! (R)
        INTEGER*4 LTEESS4 (1)   ! (R)
        INTEGER*4 LTNESS4 (1)   ! (R)
        INTEGER*4 LTSESS4 (1)   ! (W)
        POINTER (PNSESS4, NSESS4)
        POINTER (PNDESS4, NDESS4)
        POINTER (PIXEESS4, IXEESS4)
        POINTER (PIXNESS4, IXNESS4)
        POINTER (PLTEESS4, LTEESS4)
        POINTER (PLTNESS4, LTNESS4)
        POINTER (PLTSESS4, LTSESS4)
        INTEGER*4 IERR4         ! (W)

        INTEGER*4 INFREQ4
        INTEGER*4 INTRET4
        REAL*4 RELRET4
        CHARACTER* (4) CHRRET

        INTEGER*4 NEXSIDS
        INTEGER*4 NEXELEM
        INTEGER*4 NEXNODE
        INTEGER*4 NEXNSNL
        INTEGER*4 NEXSSEL
        INTEGER*4 I

        IDEXO4 = IDEXO
        INFREQ4 = EXSIDS
        CALL EXINQ4 (IDEXO4, INFREQ4, INTRET4, RELRET4, CHRRET, IERR4)
        NEXSIDS = INTRET4
        INFREQ4 = EXELEM
        CALL EXINQ4 (IDEXO4, INFREQ4, INTRET4, RELRET4, CHRRET, IERR4)
        NEXELEM = INTRET4
        INFREQ4 = EXNODE
        CALL EXINQ4 (IDEXO4, INFREQ4, INTRET4, RELRET4, CHRRET, IERR4)
        NEXNODE = INTRET4
        NEXNSNL = 0
        NEXSSEL = 0
        DO I = 1, NEXSIDS
                NEXNSNL = NEXNSNL + NDESS(I)
                NEXSSEL = NEXSSEL + NSESS(I)
        END DO
        CALL I8I4 (NEXSIDS, NSESS, PNSESS4)
        CALL I8I4 (NEXSIDS, NDESS, PNDESS4)
        CALL I8I4 (NEXSIDS, IXEESS, PIXEESS4)
        CALL I8I4 (NEXSIDS, IXNESS, PIXNESS4)
        CALL I8I4 (NEXSSEL, LTEESS, PLTEESS4)
        CALL I8I4 (NEXNSNL, LTNESS, PLTNESS4)
        CALL I4ALLOC (NEXSSEL, PLTSESS4)
        CALL EXCN2S4 (IDEXO4, NSESS4, NDESS4, IXEESS4, IXNESS4,
     &          LTEESS4, LTNESS4, LTSESS4, IERR4)
        CALL I4I8 (NEXSSEL, PLTSESS4, LTSESS)
        CALL I4FREE (PNSESS4)
        CALL I4FREE (PNDESS4)
        CALL I4FREE (PIXEESS4)
        CALL I4FREE (PIXNESS4)
        CALL I4FREE (PLTEESS4)
        CALL I4FREE (PLTNESS4)
        IERR = IERR4
        END

        INTEGER FUNCTION EXCRE (PATH, ICMODE, ICOMPWS, IOWS, IERR)
C
C       CREATE EXODUS II FILE
C
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

        ICMODE4 = ICMODE

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
C
C       ERROR REPORTING
C
        CHARACTER* (*) MODNAM   ! (R)
        CHARACTER* (*) MSG      ! (R)
        INTEGER ERRNUM          ! (R)

        INTEGER*4 ERRNUM4       ! (R)

        ERRNUM4 = ERRNUM
        CALL EXERR4 (MODNAM, MSG, ERRNUM4)
        END

        SUBROUTINE EXGATM (IDEXO, TIME, IERR)
C
C       READ ALL TIME VALUES
C
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
C
C       READ CONCATENATED NODE SETS
C
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
        INTEGER*4 IDNPSS4 (1)   ! (W)
        INTEGER*4 NNNPS4 (1)    ! (W)
        INTEGER*4 NDNPS4 (1)    ! (W)
        INTEGER*4 IXNNPS4 (1)   ! (W)
        INTEGER*4 IXDNPS4 (1)   ! (W)
        INTEGER*4 LTNNPS4 (1)   ! (W)
        POINTER (PIDNPSS4, IDNPSS4)
        POINTER (PNNNPS4, NNNPS4)
        POINTER (PNDNPS4, NDNPS4)
        POINTER (PIXNNPS4, IXNNPS4)
        POINTER (PIXDNPS4, IXDNPS4)
        POINTER (PLTNNPS4, LTNNPS4)
        INTEGER*4 IERR4         ! (W)

        INTEGER*4 INFREQ4
        INTEGER*4 INTRET4
        REAL*4 RELRET4
        CHARACTER* (4) CHRRET

        INTEGER*4 NEXNODS
        INTEGER*4 NEXNSNL

        IDEXO4 = IDEXO
        INFREQ4 = EXNODS
        CALL EXINQ4 (IDEXO4, INFREQ4, INTRET4, RELRET4, CHRRET, IERR4)
        NEXNODS = INTRET4
        INFREQ4 = EXNSNL
        CALL EXINQ4 (IDEXO4, INFREQ4, INTRET4, RELRET4, CHRRET, IERR4)
        NEXNSNL = INTRET4
        CALL I4ALLOC (NEXNODS, PIDNPSS4)
        CALL I4ALLOC (NEXNODS, PNNNPS4)
        CALL I4ALLOC (NEXNODS, PNDNPS4)
        CALL I4ALLOC (NEXNODS, PIXNNPS4)
        CALL I4ALLOC (NEXNODS, PIXDNPS4)
        CALL I4ALLOC (NEXNSNL, PLTNNPS4)
        CALL EXGCNS4 (IDEXO4, IDNPSS4, NNNPS4, NDNPS4, IXNNPS4,
     &                IXDNPS4, LTNNPS4, FACNPS, IERR4)
        CALL I4I8 (NEXNODS, PIDNPSS4, IDNPSS)
        CALL I4I8 (NEXNODS, PNNNPS4, NNNPS)
        CALL I4I8 (NEXNODS, PNDNPS4, NDNPS)
        CALL I4I8 (NEXNODS, PIXNNPS4, IXNNPS)
        CALL I4I8 (NEXNODS, PIXDNPS4, IXDNPS)
        CALL I4I8 (NEXNSNL, PLTNNPS4, LTNNPS)
        IERR = IERR4
        END

        SUBROUTINE EXGFRM (IDEXO, NFRAME, IDSCF, COORD, ITAGS, IERR)
C
C       READ COORDINATE FRAMES
C
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
        integer*4 idscf4(1)     ! (W)
        integer*4 itags4(1)     ! (W)
        INTEGER*4 IERR4         ! (W)

        POINTER (PIDSCF4, IDSCF4)
        POINTER (PITAGS4, ITAGS4)

        INTEGER*4 INFREQ4
        INTEGER*4 INTRET4
        REAL*4 RELRET4
        CHARACTER* (4) CHRRET

        idexo4 = idexo
        infreq4 = EXNCF

        CALL EXINQ4 (IDEXO4, INFREQ4, INTRET4, RELRET4, CHRRET, IERR4)

        NFRAME4 = INTRET4
        NFRAME = NFRAME4
        CALL I4ALLOC (NFRAME4, PIDSCF4)
        CALL I4ALLOC (NFRAME4, PITAGS4)

        call exgfrm4(idexo4, nframe4, idscf4, coord, itags4, ierr4)

        CALL I4I8 (NFRAME4, PIDSCF4, IDSCF)
        CALL I4I8 (NFRAME4, PITAGS4, ITAGS)
        ierr = ierr4

        END

        SUBROUTINE EXPFRM (IDEXO, NFRAME, IDSCF, COORD, ITAGS, IERR)
C
C       DEFINE/WRITE COORDINATE FRAMES
C
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
        integer*4 idscf4(1)     ! (R)
        integer*4 itags4(1)     ! (R)
        INTEGER*4 IERR4         ! (W)

        POINTER (PIDSCF4, IDSCF4)
        POINTER (PITAGS4, ITAGS4)

        idexo4 = idexo
        NFRAME4 = NFRAME

        CALL I8I4 (NFRAME4, IDSCF, PIDSCF4)
        CALL I8I4 (NFRAME4, ITAGS, PITAGS4)

        call expfrm4(idexo4, nframe4, idscf4, coord, itags4, ierr4)

        CALL I4FREE (PIDSCF4)
        CALL I4FREE (PITAGS4)
        ierr = ierr4

        END

        SUBROUTINE EXGCON (IDEXO, NAMECO, IERR)
C
C       READ COORDINATE NAMES
C
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
C
C       READ NODAL COORDINATES
C
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
C
C       READ CONCATENATED SIDE SETS
C
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
        INTEGER*4 IDESSS4 (1)   ! (W)
        INTEGER*4 NSESS4 (1)    ! (W)
        INTEGER*4 NDESS4 (1)    ! (W)
        INTEGER*4 IXEESS4 (1)   ! (W)
        INTEGER*4 IXDESS4 (1)   ! (W)
        INTEGER*4 LTEESS4 (1)   ! (W)
        INTEGER*4 LTSESS4 (1)   ! (W)
        POINTER (PIDESSS4, IDESSS4)
        POINTER (PNSESS4, NSESS4)
        POINTER (PNDESS4, NDESS4)
        POINTER (PIXEESS4, IXEESS4)
        POINTER (PIXDESS4, IXDESS4)
        POINTER (PLTEESS4, LTEESS4)
        POINTER (PLTSESS4, LTSESS4)
        INTEGER*4 IERR4         ! (W)

        INTEGER*4 INFREQ4
        INTEGER*4 INTRET4
        REAL*4 RELRET4
        CHARACTER* (4) CHRRET

        INTEGER*4 NEXSIDS
        INTEGER*4 NEXSSEL
        INTEGER*4 NEXSSDF

        IDEXO4 = IDEXO
        INFREQ4 = EXSIDS
        CALL EXINQ4 (IDEXO4, INFREQ4, INTRET4, RELRET4, CHRRET, IERR4)
        NEXSIDS = INTRET4
        INFREQ4 = EXSSEL
        CALL EXINQ4 (IDEXO4, INFREQ4, INTRET4, RELRET4, CHRRET, IERR4)
        NEXSSEL = INTRET4
        INFREQ4 = EXSSDF
        CALL EXINQ4 (IDEXO4, INFREQ4, INTRET4, RELRET4, CHRRET, IERR4)
        NEXSSDF = INTRET4
        CALL I4ALLOC (NEXSIDS, PIDESSS4)
        CALL I4ALLOC (NEXSIDS, PNSESS4)
        CALL I4ALLOC (NEXSIDS, PNDESS4)
        CALL I4ALLOC (NEXSIDS, PIXEESS4)
        CALL I4ALLOC (NEXSIDS, PIXDESS4)
        CALL I4ALLOC (NEXSSEL, PLTEESS4)
        CALL I4ALLOC (NEXSSEL, PLTSESS4)
        CALL EXGCSS4 (IDEXO4, IDESSS4, NSESS4, NDESS4, IXEESS4,
     &          IXDESS4, LTEESS4, LTSESS4, FACESS, IERR4)
        CALL I4I8 (NEXSIDS, PIDESSS4, IDESSS)
        CALL I4I8 (NEXSIDS, PNSESS4, NSESS)
        CALL I4I8 (NEXSIDS, PNDESS4, NDESS)
        CALL I4I8 (NEXSIDS, PIXEESS4, IXEESS)
        CALL I4I8 (NEXSIDS, PIXDESS4, IXDESS)
        CALL I4I8 (NEXSSEL, PLTEESS4, LTEESS)
        CALL I4I8 (NEXSSEL, PLTSESS4, LTSESS)

        IERR = IERR4
        END

        SUBROUTINE EXGCSSF (IDEXO, IDESSS, NSESS, NDESS, IXEESS,
     &          IXDESS, LTEESS, LTSESS, IERR)
C
C       READ CONCATENATED SIDE SETS
C
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
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IDESSS4 (1)   ! (W)
        INTEGER*4 NSESS4 (1)    ! (W)
        INTEGER*4 NDESS4 (1)    ! (W)
        INTEGER*4 IXEESS4 (1)   ! (W)
        INTEGER*4 IXDESS4 (1)   ! (W)
        INTEGER*4 LTEESS4 (1)   ! (W)
        INTEGER*4 LTSESS4 (1)   ! (W)
        POINTER (PIDESSS4, IDESSS4)
        POINTER (PNSESS4, NSESS4)
        POINTER (PNDESS4, NDESS4)
        POINTER (PIXEESS4, IXEESS4)
        POINTER (PIXDESS4, IXDESS4)
        POINTER (PLTEESS4, LTEESS4)
        POINTER (PLTSESS4, LTSESS4)
        INTEGER*4 IERR4         ! (W)

        INTEGER*4 INFREQ4
        INTEGER*4 INTRET4
        REAL*4 RELRET4
        CHARACTER* (4) CHRRET

        INTEGER*4 NEXSIDS
        INTEGER*4 NEXSSEL
        INTEGER*4 NEXSSDF

        IDEXO4 = IDEXO
        INFREQ4 = EXSIDS
        CALL EXINQ4 (IDEXO4, INFREQ4, INTRET4, RELRET4, CHRRET, IERR4)
        NEXSIDS = INTRET4
        INFREQ4 = EXSSEL
        CALL EXINQ4 (IDEXO4, INFREQ4, INTRET4, RELRET4, CHRRET, IERR4)
        NEXSSEL = INTRET4
        INFREQ4 = EXSSDF
        CALL EXINQ4 (IDEXO4, INFREQ4, INTRET4, RELRET4, CHRRET, IERR4)
        NEXSSDF = INTRET4
        CALL I4ALLOC (NEXSIDS, PIDESSS4)
        CALL I4ALLOC (NEXSIDS, PNSESS4)
        CALL I4ALLOC (NEXSIDS, PNDESS4)
        CALL I4ALLOC (NEXSIDS, PIXEESS4)
        CALL I4ALLOC (NEXSIDS, PIXDESS4)
        CALL I4ALLOC (NEXSSEL, PLTEESS4)
        CALL I4ALLOC (NEXSSEL, PLTSESS4)
        CALL EXGCSSF4 (IDEXO4, IDESSS4, NSESS4, NDESS4, IXEESS4,
     &          IXDESS4, LTEESS4, LTSESS4, IERR4)
        CALL I4I8 (NEXSIDS, PIDESSS4, IDESSS)
        CALL I4I8 (NEXSIDS, PNSESS4, NSESS)
        CALL I4I8 (NEXSIDS, PNDESS4, NDESS)
        CALL I4I8 (NEXSIDS, PIXEESS4, IXEESS)
        CALL I4I8 (NEXSIDS, PIXDESS4, IXDESS)
        CALL I4I8 (NEXSSEL, PLTEESS4, LTEESS)
        CALL I4I8 (NEXSSEL, PLTSESS4, LTSESS)

        IERR = IERR4
        END

        SUBROUTINE EXGEAT (IDEXO, IDELB, ATRIB, IERR)
C
C       READ ELEMENT BLOCK ATTRIBUTES
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDELB           ! (R)
        REAL ATRIB (*)          ! (NUMATR,NUMELB) ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IDELB4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        IDELB4 = IDELB
        CALL EXGEAT4 (IDEXO4, IDELB4, ATRIB, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGEAN (IDEXO, IDELB, NATTR, NAMES, IERR)
C
C       READ ELEMENT BLOCK ATTRIBUTES
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDELB           ! (R)
        INTEGER NATTR           ! (R)
        CHARACTER* (*) NAMES (*) ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IDELB4        ! (R)
        INTEGER*4 NATTR4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        IDELB4 = IDELB
        NATTR4 = NATTR
        CALL EXGEAN4 (IDEXO4, IDELB4, NATTR4, NAMES, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGEBI (IDEXO, IDELBS, IERR)
C
C       READ ELEMENT BLOCK IDS
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDELBS (*)      ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IDELBS4 (1)   ! (W)
        POINTER (PIDELBS4, IDELBS4)
        INTEGER*4 IERR4         ! (W)

        INTEGER*4 INFREQ4
        INTEGER*4 INTRET4
        REAL*4 RELRET4
        CHARACTER* (4) CHRRET

        INTEGER*4 NEXELBL

        IDEXO4 = IDEXO
        INFREQ4 = EXELBL
        CALL EXINQ4 (IDEXO4, INFREQ4, INTRET4, RELRET4, CHRRET, IERR4)
        NEXELBL = INTRET4
        CALL I4ALLOC (NEXELBL, PIDELBS4)
        CALL EXGEBI4 (IDEXO4, IDELBS4, IERR4)
        CALL I4I8 (NEXELBL, PIDELBS4, IDELBS)
        IERR = IERR4
        END

        SUBROUTINE EXGELB (IDEXO, IDELB, NAMELB, NUMELB, NUMLNK,
     &          NUMATR, IERR)
C
C       READ ELEMENT BLOCK PARAMETERS
C
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
        INTEGER*4 IDELB4        ! (R)
        INTEGER*4 NUMELB4       ! (W)
        INTEGER*4 NUMLNK4       ! (W)
        INTEGER*4 NUMATR4       ! (W)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        IDELB4 = IDELB
        CALL EXGELB4 (IDEXO4, IDELB4, NAMELB, NUMELB4, NUMLNK4,
     &          NUMATR4, IERR4)
        NUMELB = NUMELB4
        NUMLNK = NUMLNK4
        NUMATR = NUMATR4
        IERR = IERR4
        END

        SUBROUTINE EXGELC (IDEXO, IDELB, LINK, IERR)
C
C       READ ELEMENT BLOCK CONNECTIVITY
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDELB           ! (R)
        INTEGER LINK (*)        ! (NUMLNK,NUMELB) ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IDELB4        ! (R)
        INTEGER*4 LINK4 (1)     ! (NUMLNK,NUMELB) ! (W)
        POINTER (PLINK4, LINK4)
        INTEGER*4 IERR4         ! (W)

        CHARACTER* (MXSTLN) NAMELB
        INTEGER*4 NUMELB4
        INTEGER*4 NUMLNK4
        INTEGER*4 NUMATR4

        IDEXO4 = IDEXO
        IDELB4 = IDELB
        CALL EXGELB4 (IDEXO4, IDELB4, NAMELB, NUMELB4, NUMLNK4,
     &          NUMATR4, IERR4)
        IF (NAMELB .EQ. 'NSIDED' .OR. NAMELB .EQ. 'nsided') THEN
          CALL I4ALLOC (NUMLNK4, PLINK4)
          CALL EXGELC4 (IDEXO4, IDELB4, LINK4, IERR4)
          CALL I4I8 (NUMLNK4, PLINK4, LINK)
        ELSE
          CALL I4ALLOC (NUMLNK4*NUMELB4, PLINK4)
          CALL EXGELC4 (IDEXO4, IDELB4, LINK4, IERR4)
          CALL I4I8 (NUMLNK4*NUMELB4, PLINK4, LINK)
        ENDIF
        IERR = IERR4
        END

        SUBROUTINE EXGECPP (IDEXO, ITYPE, IDELB, COUNTS, IERR)
C
C       READ ENTITY_COUNT-PER-POLYHEDRA
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDELB           ! (R)
        INTEGER ITYPE           ! (R)
        INTEGER COUNTS (*)      ! (NUMELB) ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IDELB4        ! (R)
        INTEGER*4 ITYPE4        ! (R)
        INTEGER*4 COUNTS4 (1)   ! (NUMELB) ! (W)
        POINTER (PCOUNTS4, COUNTS4)
        INTEGER*4 IERR4         ! (W)

        CHARACTER* (MXSTLN) NAMELB
        INTEGER*4 NUMELB4
        INTEGER*4 NUMLNK4
        INTEGER*4 NUMATR4

        IDEXO4 = IDEXO
        ITYPE4 = ITYPE
        IDELB4 = IDELB
        CALL EXGELB4 (IDEXO4, IDELB4, NAMELB, NUMELB4, NUMLNK4,
     &          NUMATR4, IERR4)
        CALL I4ALLOC (NUMELB4, PCOUNTS4)
        CALL EXGECPP4 (IDEXO4, ITYPE4, IDELB4, COUNTS4, IERR4)
        CALL I4I8 (NUMELB4, PCOUNTS4, COUNTS)
        IERR = IERR4
        END

        SUBROUTINE EXGENM (IDEXO, MAPEL, IERR)
C
C       READ ELEMENT NUMBER MAP
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER MAPEL (*)       ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 MAPEL4 (1)    ! (W)
        POINTER (PMAPEL4, MAPEL4)
        INTEGER*4 IERR4         ! (W)

        INTEGER*4 INFREQ4
        INTEGER*4 INTRET4
        REAL*4 RELRET4
        CHARACTER* (4) CHRRET

        INTEGER*4 NEXELEM

        IDEXO4 = IDEXO
        INFREQ4 = EXELEM
        CALL EXINQ4 (IDEXO4, INFREQ4, INTRET4, RELRET4, CHRRET, IERR4)
        NEXELEM = INTRET4
        CALL I4ALLOC (NEXELEM, PMAPEL4)
        CALL EXGENM4 (IDEXO4, MAPEL4, IERR4)
        CALL I4I8 (NEXELEM, PMAPEL4, MAPEL)
        IERR = IERR4
        END

        SUBROUTINE EXGEM (IDEXO, MAPID, MAPEL, IERR)
C
C       READ ELEMENT ORDER MAP
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER MAPID           ! (R)
        INTEGER MAPEL (*)       ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 MAPID4        ! (R)
        INTEGER*4 MAPEL4 (1)    ! (W)
        POINTER (PMAPEL4, MAPEL4)
        INTEGER*4 IERR4         ! (W)

        INTEGER*4 INFREQ4
        INTEGER*4 INTRET4
        REAL*4 RELRET4
        CHARACTER* (4) CHRRET

        INTEGER*4 NEXELEM

        IDEXO4 = IDEXO
        MAPID4 = MAPID
        INFREQ4 = EXELEM
        CALL EXINQ4 (IDEXO4, INFREQ4, INTRET4, RELRET4, CHRRET, IERR4)
        NEXELEM = INTRET4
        CALL I4ALLOC (NEXELEM, PMAPEL4)
        CALL EXGEM4 (IDEXO4, MAPID4, MAPEL4, IERR4)
        CALL I4I8 (NEXELEM, PMAPEL4, MAPEL)
        IERR = IERR4
        END

        SUBROUTINE EXGPEM (IDEXO, MAPID, ISTART, ICOUNT, MAPEL, IERR)
C
C       READ ELEMENT ORDER MAP
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER MAPID           ! (R)
        INTEGER ISTART          ! (R)
        INTEGER ICOUNT          ! (R)
        INTEGER MAPEL (*)       ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 MAPID4        ! (R)
        INTEGER*4 ISTART4       ! (R)
        INTEGER*4 ICOUNT4       ! (R)
        INTEGER*4 MAPEL4 (1)    ! (W)
        POINTER (PMAPEL4, MAPEL4)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        MAPID4 = MAPID
        ISTART4 = ISTART
        ICOUNT4 = ICOUNT
        CALL I4ALLOC (ICOUNT4, PMAPEL4)
        CALL EXGPEM4 (IDEXO4, MAPID4, ISTART4, ICOUNT4, MAPEL4, IERR4)
        CALL I4I8 (ICOUNT4, PMAPEL4, MAPEL)
        IERR = IERR4
        END

        SUBROUTINE EXGEV (IDEXO, ISTEP, IXELEV, IDELB, NUMELB,
     &          VALEV, IERR)
C
C       READ ELEMENT VARIABLE VALUES AT A TIME STEP
C
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
        INTEGER*4 IDELB4        ! (R)
        INTEGER*4 NUMELB4       ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        ISTEP4 = ISTEP
        IXELEV4 = IXELEV
        IDELB4 = IDELB
        NUMELB4 = NUMELB
        CALL EXGEV4 (IDEXO4, ISTEP4, IXELEV4, IDELB4, NUMELB4,
     &          VALEV, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGEVT (IDEXO, IXELEV, IELNUM, ISTPB, ISTPE,
     &          VALEV, IERR)
C
C       READ ELEMENT VARIABLE VALUES THROUGH TIME
C
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
        INTEGER*4 IELNUM4       ! (R)
        INTEGER*4 ISTPB4        ! (R)
        INTEGER*4 ISTPE4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4  = IDEXO
        IXELEV4 = IXELEV
        IELNUM4 = IELNUM
        ISTPB4  = ISTPB
        ISTPE4  = ISTPE

        CALL EXGEVT4 (IDEXO4, IXELEV4, IELNUM4, ISTPB4, ISTPE4,
     &          VALEV, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGNSV (IDEXO, ISTEP, IXV, ID, NUM, VAL, IERR)
C
C       READ NODESET VARIABLE VALUES AT A TIME STEP
C
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
        INTEGER*4 ID4   ! (R)
        INTEGER*4 NUM4  ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        ISTEP4 = ISTEP
        IXV4   = IXV
        ID4    = ID
        NUM4   = NUM
        CALL EXGNSV4 (IDEXO4, ISTEP4, IXV4, ID4, NUM4, VAL, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGSSV (IDEXO, ISTEP, IXV, ID, NUM, VAL, IERR)
C
C       READ SIDESET VARIABLE VALUES AT A TIME STEP
C
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
        INTEGER*4 ID4   ! (R)
        INTEGER*4 NUM4  ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        ISTEP4 = ISTEP
        IXV4   = IXV
        ID4    = ID
        NUM4   = NUM
        CALL EXGSSV4 (IDEXO4, ISTEP4, IXV4, ID4, NUM4, VAL, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGGV (IDEXO, ISTEP, NVARGL, VALGV, IERR)
C
C       READ GLOBAL VARIABLE VALUES AT A TIME STEP
C
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
C
C       READ GLOBAL VARIABLE VALUES THROUGH TIME
C
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
C
C       READ INFORMATION RECORDS
C
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
C
C       READ INITIALIZATION PARAMETERS
C
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
        INTEGER*4 NDIM4         ! (W)
        INTEGER*4 NUMNP4        ! (W)
        INTEGER*4 NUMEL4        ! (W)
        INTEGER*4 NELBLK4       ! (W)
        INTEGER*4 NUMNPS4       ! (W)
        INTEGER*4 NUMESS4       ! (W)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        CALL EXGINI4 (IDEXO4, TITLE, NDIM4, NUMNP4, NUMEL4, NELBLK4,
     &          NUMNPS4, NUMESS4, IERR4)
        NDIM = NDIM4
        NUMNP = NUMNP4
        NUMEL = NUMEL4
        NELBLK = NELBLK4
        NUMNPS = NUMNPS4
        NUMESS = NUMESS4
        IERR = IERR4
        END

        SUBROUTINE EXGMP (IDEXO, NNMAP, NEMAP, IERR)
C
C       GET MAP PARAMETERS
C
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
C
C       PUT MAP PARAMETERS
C
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
C
C       READ ELEMENT ORDER MAP
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER MAPEL (*)       ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 MAPEL4 (1)    ! (W)
        POINTER (PMAPEL4, MAPEL4)
        INTEGER*4 IERR4         ! (W)

        INTEGER*4 INFREQ4
        INTEGER*4 INTRET4
        REAL*4 RELRET4
        CHARACTER* (4) CHRRET

        INTEGER*4 NEXELEM

        IDEXO4 = IDEXO
        INFREQ4 = EXELEM
        CALL EXINQ4 (IDEXO4, INFREQ4, INTRET4, RELRET4, CHRRET, IERR4)
        NEXELEM = INTRET4
        CALL I4ALLOC (NEXELEM, PMAPEL4)
        CALL EXGMAP4 (IDEXO4, MAPEL4, IERR4)
        CALL I4I8 (NEXELEM, PMAPEL4, MAPEL)
        IERR = IERR4
        END

        SUBROUTINE EXGNNM (IDEXO, MAPNOD, IERR)
C
C       READ NODE NUMBER MAP
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER MAPNOD (*)      ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 MAPNOD4 (1)   ! (W)
        POINTER (PMAPNOD4, MAPNOD4)
        INTEGER*4 IERR4         ! (W)

        INTEGER*4 INFREQ4
        INTEGER*4 INTRET4
        REAL*4 RELRET4
        CHARACTER* (4) CHRRET

        INTEGER*4 NEXNODE

        IDEXO4 = IDEXO
        INFREQ4 = EXNODE
        CALL EXINQ4 (IDEXO4, INFREQ4, INTRET4, RELRET4, CHRRET, IERR4)
        NEXNODE = INTRET4
        CALL I4ALLOC (NEXNODE, PMAPNOD4)
        CALL EXGNNM4 (IDEXO4, MAPNOD4, IERR4)
        CALL I4I8 (NEXNODE, PMAPNOD4, MAPNOD)
        IERR = IERR4
        END

        SUBROUTINE EXGNM (IDEXO, MAPID, MAPNOD, IERR)
C
C       READ NODE ORDER MAP
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER MAPID           ! (R)
        INTEGER MAPNOD (*)      ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 MAPID4        ! (R)
        INTEGER*4 MAPNOD4 (1)   ! (W)
        POINTER (PMAPNOD4, MAPNOD4)
        INTEGER*4 IERR4         ! (W)

        INTEGER*4 INFREQ4
        INTEGER*4 INTRET4
        REAL*4 RELRET4
        CHARACTER* (4) CHRRET

        INTEGER*4 NEXNODE

        IDEXO4 = IDEXO
        MAPID4 = MAPID
        INFREQ4 = EXNODE
        CALL EXINQ4 (IDEXO4, INFREQ4, INTRET4, RELRET4, CHRRET, IERR4)
        NEXNODE = INTRET4
        CALL I4ALLOC (NEXNODE, PMAPNOD4)
        CALL EXGNM4 (IDEXO4, MAPID4, MAPNOD4, IERR4)
        CALL I4I8 (NEXNODE, PMAPNOD4, MAPNOD)
        IERR = IERR4
        END

        SUBROUTINE EXGNP (IDEXO, IDNPS, NNNPS, NDNPS, IERR)
C
C       READ NODE SET PARAMETERS
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDNPS           ! (R)
        INTEGER NNNPS           ! (W)
        INTEGER NDNPS           ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IDNPS4        ! (R)
        INTEGER*4 NNNPS4        ! (W)
        INTEGER*4 NDNPS4        ! (W)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        IDNPS4 = IDNPS
        CALL EXGNP4 (IDEXO4, IDNPS4, NNNPS4, NDNPS4, IERR4)
        NNNPS = NNNPS4
        NDNPS = NDNPS4
        IERR = IERR4
        END

        SUBROUTINE EXGNS (IDEXO, IDNPS, LTNNPS, IERR)
C
C       READ NODE SET
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDNPS           ! (R)
        INTEGER LTNNPS (*)      ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IDNPS4        ! (R)
        INTEGER*4 LTNNPS4 (1)   ! (W)
        POINTER (PLTNNPS4, LTNNPS4)
        INTEGER*4 IERR4         ! (W)

        INTEGER*4 NNNPS4
        INTEGER*4 NDNPS4

        IDEXO4 = IDEXO
        IDNPS4 = IDNPS
        CALL EXGNP4 (IDEXO4, IDNPS4, NNNPS4, NDNPS4, IERR4)
        CALL I4ALLOC (NNNPS4, PLTNNPS4)
        CALL EXGNS4 (IDEXO4, IDNPS4, LTNNPS4, IERR4)
        CALL I4I8 (NNNPS4, PLTNNPS4, LTNNPS)
        IERR = IERR4
        END

        SUBROUTINE EXGNSD (IDEXO, IDNPS, FACNPS, IERR)
C
C       READ NODE SET DISTRIBUTION FACTORS
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDNPS           ! (R)
        REAL FACNPS (*)         ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IDNPS4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        IDNPS4 = IDNPS

        CALL EXGNSD4 (IDEXO4, IDNPS4, FACNPS, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGNSI (IDEXO, IDNPSS, IERR)
C
C       READ NODE SET IDS
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDNPSS (*)      ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IDNPSS4 (1)   ! (W)
        POINTER (PIDNPSS4, IDNPSS4)
        INTEGER*4 IERR4         ! (W)

        INTEGER*4 INFREQ4
        INTEGER*4 INTRET4
        REAL*4 RELRET4
        CHARACTER* (4) CHRRET

        INTEGER*4 NEXNODS

        IDEXO4 = IDEXO
        INFREQ4 = EXNODS
        CALL EXINQ4 (IDEXO4, INFREQ4, INTRET4, RELRET4, CHRRET, IERR4)
        NEXNODS = INTRET4
        CALL I4ALLOC (NEXNODS, PIDNPSS4)
        CALL EXGNSI4 (IDEXO4, IDNPSS4, IERR4)
        CALL I4I8 (NEXNODS, PIDNPSS4, IDNPSS)
        IERR = IERR4
        END

        SUBROUTINE EXGNV (IDEXO, ISTEP, IXNODV, NUMNP, VALNV, IERR)
C
C       READ NODAL VARIABLE VALUES AT A TIME STEP
C
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
        INTEGER*4 NUMNP4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4  = IDEXO
        ISTEP4  = ISTEP
        IXNODV4 = IXNODV
        NUMNP4  = NUMNP
        CALL EXGNV4 (IDEXO4, ISTEP4, IXNODV4, NUMNP4, VALNV, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGNVT (IDEXO, IXNODV, NODNUM, ISTPB, ISTPE,
     &          VALNV, IERR)
C
C       READ NODAL VARIABLE VALUES THROUGH TIME
C
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
        INTEGER*4 NODNUM4       ! (R)
        INTEGER*4 ISTPB4        ! (R)
        INTEGER*4 ISTPE4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4  = IDEXO
        IXNODV4 = IXNODV
        NODNUM4 = NODNUM
        ISTPB4  = ISTPB
        ISTPE4  = ISTPE
        CALL EXGNVT4 (IDEXO4, IXNODV4, NODNUM4, ISTPB4, ISTPE4,
     &          VALNV, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGP (IDEXO, ITYPE, ID, NAMEPR, IVAL, IERR)
C
C       READ OBJECT PROPERTY
C
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
        INTEGER*4 ID4           ! (R)
        INTEGER*4 IVAL4         ! (W)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        ITYPE4 = ITYPE
        ID4 = ID
        CALL EXGP4 (IDEXO4, ITYPE4, ID4, NAMEPR, IVAL4, IERR4)
        IVAL = IVAL4
        IERR = IERR4
        END

        SUBROUTINE EXGPA (IDEXO, ITYPE, NAMEPR, IVAL, IERR)
C
C       READ OBJECT PROPERTY ARRAY
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER ITYPE           ! (R)
        CHARACTER* (*) NAMEPR   ! (R)
        INTEGER IVAL (*)        ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 ITYPE4        ! (R)
        INTEGER*4 IVAL4 (1)     ! (W)
        POINTER (PIVAL4, IVAL4)
        INTEGER*4 IERR4         ! (W)

        INTEGER*4 INFREQ4
        INTEGER*4 INTRET4
        REAL*4 RELRET4
        CHARACTER* (4) CHRRET

        INTEGER*4 NUMOBJ

        IDEXO4 = IDEXO
        ITYPE4 = ITYPE
        IF (ITYPE4 .EQ. EXEBLK) THEN
                INFREQ4 = EXELBL
        ELSE IF (ITYPE4 .EQ. EXNSET) THEN
                INFREQ4 = EXNODS
        ELSE IF (ITYPE4 .EQ. EXSSET) THEN
                INFREQ4 = EXSIDS
        ELSE IF (ITYPE4 .EQ. EXEMAP) THEN
                INFREQ4 = EXNEM
        ELSE IF (ITYPE4 .EQ. EXNMAP) THEN
                INFREQ4 = EXNNM
        END IF
        CALL EXINQ4 (IDEXO4, INFREQ4, INTRET4, RELRET4, CHRRET, IERR4)
        NUMOBJ = INTRET4
        CALL I4ALLOC (NUMOBJ, PIVAL4)
        CALL EXGPA4 (IDEXO4, ITYPE4, NAMEPR, IVAL4, IERR4)
        CALL I4I8 (NUMOBJ, PIVAL4, IVAL)
        IERR = IERR4
        END

        SUBROUTINE EXGPN (IDEXO, ITYPE, NAMEPR, IERR)
C
C       READ PROPERTY ARRAY NAMES
C
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
C
C       READ QA RECORDS
C
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
C
C       READ SIDE SET PARAMETERS
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDESS           ! (R)
        INTEGER NSESS           ! (W)
        INTEGER NDESS           ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IDESS4        ! (R)
        INTEGER*4 NSESS4        ! (W)
        INTEGER*4 NDESS4        ! (W)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        IDESS4 = IDESS
        CALL EXGSP4 (IDEXO4, IDESS4, NSESS4, NDESS4, IERR4)
        NSESS = NSESS4
        NDESS = NDESS4
        IERR = IERR4
        END

        SUBROUTINE EXGSNL (IDEXO, IDESS, NNESS, IERR)
C
C       GET SIDE SET NODELIST LENGTH
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDESS           ! (R)
        INTEGER NNESS           ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IDESS4        ! (R)
        INTEGER*4 NNESS4        ! (W)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        IDESS4 = IDESS
        CALL EXGSNL4 (IDEXO4, IDESS4, NNESS4, IERR4)
        NNESS = NNESS4
        IERR = IERR4
        END

        SUBROUTINE EXGSS (IDEXO, IDESS, LTEESS, LTSESS, IERR)
C
C       READ SIDE SET
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDESS           ! (R)
        INTEGER LTEESS (*)      ! (W)
        INTEGER LTSESS (*)      ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IDESS4        ! (R)
        INTEGER*4 LTEESS4 (1)   ! (W)
        INTEGER*4 LTSESS4 (1)   ! (W)
        POINTER (PLTEESS4, LTEESS4)
        POINTER (PLTSESS4, LTSESS4)
        INTEGER*4 IERR4         ! (W)

        INTEGER*4 NSESS4
        INTEGER*4 NDESS4

        IDEXO4 = IDEXO
        IDESS4 = IDESS
        CALL EXGSP4 (IDEXO4, IDESS4, NSESS4, NDESS4, IERR4)
        CALL I4ALLOC (NSESS4, PLTEESS4)
        CALL I4ALLOC (NSESS4, PLTSESS4)
        CALL EXGSS4 (IDEXO4, IDESS4, LTEESS4, LTSESS4, IERR4)
        CALL I4I8 (NSESS4, PLTEESS4, LTEESS)
        CALL I4I8 (NSESS4, PLTSESS4, LTSESS)
        IERR = IERR4
        END

        SUBROUTINE EXGSSD (IDEXO, IDESS, FACESS, IERR)
C
C       READ SIDE SET DISTRIBUTION FACTORS
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDESS           ! (R)
        REAL FACESS (*)         ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IDESS4        ! (R)
        INTEGER*4 IERR4         ! (W)


        IDEXO4 = IDEXO
        IDESS4 = IDESS
        CALL EXGSSD4 (IDEXO4, IDESS4, FACESS, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXGSSI (IDEXO, IDESSS, IERR)
C
C       READ SIDE SET IDS
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDESSS (*)      ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IDESSS4 (1)   ! (W)
        POINTER (PIDESSS4, IDESSS4)
        INTEGER*4 IERR4         ! (W)

        INTEGER*4 INFREQ4
        INTEGER*4 INTRET4
        REAL*4 RELRET4
        CHARACTER* (4) CHRRET
        INTEGER*4 NEXSIDS

        IDEXO4 = IDEXO
        INFREQ4 = EXSIDS
        CALL EXINQ4 (IDEXO4, INFREQ4, INTRET4, RELRET4, CHRRET, IERR4)
        NEXSIDS = INTRET4
        CALL I4ALLOC (NEXSIDS, PIDESSS4)
        CALL EXGSSI4 (IDEXO4, IDESSS4, IERR4)
        CALL I4I8 (NEXSIDS, PIDESSS4, IDESSS)
        IERR = IERR4
        END

        SUBROUTINE EXGSSN (IDEXO, IDESS, INCNT, LTNESS, IERR)
C
C       READ SIDE SET NODE LIST
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDESS           ! (R)
        INTEGER INCNT (*)       ! (W)
        INTEGER LTNESS (*)      ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IDESS4        ! (R)
        INTEGER*4 INCNT4 (1)    ! (W)
        INTEGER*4 LTNESS4 (1)   ! (W)
        INTEGER*4 IERR4         ! (W)

        INTEGER*4 TMP4(1)
        POINTER (PTMP4, TMP4)
        POINTER (PINCNT4, INCNT4)
        POINTER (PLTNESS4, LTNESS4)

        INTEGER*4 NSESS4
        INTEGER*4 NDESS4
        INTEGER*4 NLTNESS

        INTEGER*4 NINCNT
        integer i

        IDEXO4 = IDEXO
        IDESS4 = IDESS
C ... NOTE: It is not permissible to assume that the number of
C           distribution factors is the same as the number of
C           nodes associated with a sideset since some codes
C           are not maintaining that associativity...
        CALL EXGSP4 (IDEXO4, IDESS4, NSESS4, NDESS4, IERR4)
        NINCNT = NSESS4

C       MUST COMPUTE NODE LIST SIZE

        NINCNT = NSESS4
        CALL I4ALLOC (NINCNT, PTMP4)

        CALL EXGSSC4 (IDEXO4, IDESS4, TMP4, IERR4)
        nltness = 0
        do 11 i=1, nincnt
           nltness = nltness + tmp4(i)
 11     continue
        CALL I4FREE (PTMP4)
        CALL I4ALLOC (NINCNT, PINCNT4)
        CALL I4ALLOC (NLTNESS, PLTNESS4)

        CALL EXGSSN4 (IDEXO4, IDESS4, INCNT4, LTNESS4, IERR4)
        CALL I4I8 (NLTNESS, PLTNESS4, LTNESS)
        CALL I4I8 (NINCNT, PINCNT4, INCNT)
        IERR = IERR4
        END

        SUBROUTINE EXGSSC (IDEXO, IDESS, INCNT, IERR)
C
C       READ SIDE SET NODE LIST
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDESS           ! (R)
        INTEGER INCNT (*)       ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IDESS4        ! (R)
        INTEGER*4 INCNT4 (1)    ! (W)
        POINTER (PINCNT4, INCNT4)
        INTEGER*4 IERR4         ! (W)

        integer*4 nsess4
        integer*4 ndess4
        INTEGER*4 NINCNT

        IDEXO4 = IDEXO
        IDESS4 = IDESS
        CALL EXGSP4 (IDEXO4, IDESS4, NSESS4, NDESS4, IERR4)
        NINCNT = NSESS4
        CALL I4ALLOC (NINCNT, PINCNT4)

        CALL EXGSSC4 (IDEXO4, IDESS4, INCNT4, IERR4)

        CALL I4I8 (NINCNT, PINCNT4, INCNT)
        IERR = IERR4
        END

        SUBROUTINE EXGCSSC (IDEXO, INCNT, IERR)
C
C       READ CONCATENATED SIDE SET NODE COUNT LIST
C       (nodes/face for all faces)
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER INCNT (*)       ! (W)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 INCNT4 (1)    ! (W)
        POINTER (PINCNT4, INCNT4)
        INTEGER*4 IERR4         ! (W)

        integer*4 nsess4
        INTEGER*4 NINCNT

        INTEGER*4 INFREQ4
        REAL*4 RELRET4
        CHARACTER* (4) CHRRET

        IDEXO4 = IDEXO
        INFREQ4 = EXSSEL
        CALL EXINQ4 (IDEXO4, INFREQ4, NSESS4, RELRET4, CHRRET, IERR4)
        NINCNT = NSESS4
        CALL I4ALLOC (NINCNT, PINCNT4)

        CALL EXGCSSC4 (IDEXO4, INCNT4, IERR4)

        CALL I4I8 (NINCNT, PINCNT4, INCNT)
        IERR = IERR4
        END

        SUBROUTINE EXGTIM (IDEXO, NSTEP, TIME, IERR)
C
C       READ TIME VALUE FOR A TIME STEP
C
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
C
C       READ RESULTS VARIABLE NAMES
C
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
C
C       READ RESULTS VARIABLE NAME
C
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
C
C       WRITE RESULTS VARIABLE NAME
C
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
C
C       READ RESULTS VARIABLE PARAMETERS
C
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
C
C       READ ELEMENT VARIABLE TRUTH TABLE
C
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
        INTEGER*4 ISEVOK4(1)    ! (NVAREL,NELBLK) (W)
        POINTER (PISEVOK4, ISEVOK4)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        NELBLK4 = NELBLK
        NVAREL4 = NVAREL
        CALL I4ALLOC (NVAREL4*NELBLK4, PISEVOK4)
        CALL EXGVTT4 (IDEXO4, NELBLK4, NVAREL4, ISEVOK4, IERR4)
        CALL I4I8 (NVAREL4*NELBLK4, PISEVOK4, ISEVOK)
        IERR = IERR4
        END

        SUBROUTINE EXGNSTT (IDEXO, NBLK, NVAR, ISVOK, IERR)
C
C       READ NODESET VARIABLE TRUTH TABLE
C
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
        INTEGER*4 ISVOK4(1)     ! (NVAR,NBLK) (W)
        POINTER (PISVOK4, ISVOK4)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        NBLK4  = NBLK
        NVAR4  = NVAR
        CALL I4ALLOC (NVAR4*NBLK4, PISVOK4)
        CALL EXGNSTT4 (IDEXO4, NBLK4, NVAR4, ISVOK4, IERR4)
        CALL I4I8 (NVAR4*NBLK4, PISVOK4, ISVOK)
        IERR = IERR4
        END

        SUBROUTINE EXGSSTT (IDEXO, NBLK, NVAR, ISVOK, IERR)
C
C       READ SIDESET VARIABLE TRUTH TABLE
C
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
        INTEGER*4 ISVOK4(1)     ! (NVAR,NBLK) (W)
        POINTER (PISVOK4, ISVOK4)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        NBLK4  = NBLK
        NVAR4  = NVAR
        CALL I4ALLOC (NVAR4*NBLK4, PISVOK4)
        CALL EXGSSTT4 (IDEXO4, NBLK4, NVAR4, ISVOK4, IERR4)
        CALL I4I8 (NVAR4*NBLK4, PISVOK4, ISVOK)
        IERR = IERR4
        END

        SUBROUTINE EXINQ (IDEXO, INFREQ, INTRET, RELRET, CHRRET, IERR)
C
C       INQUIRE EXODUS PARAMETERS
C
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
        INTEGER*4 INTRET4       ! (W)
        REAL*4 RELRET4          ! (W)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        INFREQ4 = INFREQ
        CALL EXINQ4 (IDEXO4, INFREQ4, INTRET4, RELRET4, CHRRET, IERR4)
        INTRET = INTRET4
        RELRET = RELRET4
        IERR = IERR4
        END

        INTEGER FUNCTION EXOPEN (PATH, IMODE, ICOMPWS, IOWS, VERS, IERR)
C
C       OPEN EXODUS II FILE
C
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

        IMODE4 = IMODE

C ... This needs to be set explicitly since all wrapper routines
C     are passing REAL*8 and the exo_jack.c routines all take double*
        ICOMPWS4 = 8
        IOWS4 = IOWS
        IRET4 = EXOPEN4 (PATH, IMODE4, ICOMPWS4, IOWS4, VERS4, IERR4)
        ICOMPWS = ICOMPWS4
        IOWS = IOWS4
        VERS = VERS4
        IERR = IERR4
        IRET = IRET4
        EXOPEN = IRET
        END

        INTEGER FUNCTION EXLGMD (IDEXO)
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
C
C       SET ERROR REPORTING LEVEL
C
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
C
C       SET MAXIMUM NAME LENGTH
C
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


        SUBROUTINE EXPCNS (IDEXO, IDNPSS, NNNPS, NDNPS, IXNNPS,
     &          IXDNPS, LTNNPS, FACNPS, IERR)
C
C       WRITE CONCATENATED NODE SETS
C
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
        INTEGER*4 IDNPSS4 (1)   ! (R)
        INTEGER*4 NNNPS4 (1)    ! (R)
        INTEGER*4 NDNPS4 (1)    ! (R)
        INTEGER*4 IXNNPS4 (1)   ! (R)
        INTEGER*4 IXDNPS4 (1)   ! (R)
        INTEGER*4 LTNNPS4 (1)   ! (R)
        POINTER (PIDNPSS4, IDNPSS4)
        POINTER (PNNNPS4, NNNPS4)
        POINTER (PNDNPS4, NDNPS4)
        POINTER (PIXNNPS4, IXNNPS4)
        POINTER (PIXDNPS4, IXDNPS4)
        POINTER (PLTNNPS4, LTNNPS4)
        INTEGER*4 IERR4         ! (W)

        INTEGER*4 INFREQ4
        INTEGER*4 INTRET4
        REAL*4 RELRET4
        CHARACTER* (4) CHRRET

        INTEGER*4 NEXNODS
        INTEGER*4 NEXNSNL
        INTEGER*4 I

        IDEXO4 = IDEXO
        INFREQ4 = EXNODS
        CALL EXINQ4 (IDEXO4, INFREQ4, INTRET4, RELRET4, CHRRET, IERR4)
        NEXNODS = INTRET4
        NEXNSNL = 0
        DO I = 1, NEXNODS
                NEXNSNL = NEXNSNL + NNNPS(I)
        END DO
        CALL I8I4 (NEXNODS, IDNPSS, PIDNPSS4)
        CALL I8I4 (NEXNODS, NNNPS, PNNNPS4)
        CALL I8I4 (NEXNODS, NDNPS, PNDNPS4)
        CALL I8I4 (NEXNODS, IXNNPS, PIXNNPS4)
        CALL I8I4 (NEXNODS, IXDNPS, PIXDNPS4)
        CALL I8I4 (NEXNSNL, LTNNPS, PLTNNPS4)
        CALL EXPCNS4 (IDEXO4, IDNPSS4, NNNPS4, NDNPS4, IXNNPS4,
     &          IXDNPS4, LTNNPS4, FACNPS, IERR4)
        CALL I4FREE (PIDNPSS4)
        CALL I4FREE (PNNNPS4)
        CALL I4FREE (PNDNPS4)
        CALL I4FREE (PIXNNPS4)
        CALL I4FREE (PIXDNPS4)
        CALL I4FREE (PLTNNPS4)
        IERR = IERR4
        END

        SUBROUTINE EXPCON (IDEXO, NAMECO, IERR)
C
C       WRITE COORDINATE NAMES
C
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
C
C       WRITE NODAL COORDINATES
C
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
C
C       WRITE CONCATENATED SIDE SETS
C
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
        INTEGER*4 IDESSS4 (1)   ! (R)
        INTEGER*4 NSESS4 (1)    ! (R)
        INTEGER*4 NDESS4 (1)    ! (R)
        INTEGER*4 IXEESS4 (1)   ! (R)
        INTEGER*4 IXDESS4 (1)   ! (R)
        INTEGER*4 LTEESS4 (1)   ! (R)
        INTEGER*4 LTSESS4 (1)   ! (R)
        POINTER (PIDESSS4, IDESSS4)
        POINTER (PNSESS4, NSESS4)
        POINTER (PNDESS4, NDESS4)
        POINTER (PIXEESS4, IXEESS4)
        POINTER (PIXDESS4, IXDESS4)
        POINTER (PLTEESS4, LTEESS4)
        POINTER (PLTSESS4, LTSESS4)
        INTEGER*4 IERR4         ! (W)

        INTEGER*4 INFREQ4
        INTEGER*4 INTRET4
        REAL*4 RELRET4
        CHARACTER* (4) CHRRET

        INTEGER*4 NEXSIDS
        INTEGER*4 NEXSSEL
        INTEGER*4 I

        IDEXO4 = IDEXO
        INFREQ4 = EXSIDS
        CALL EXINQ4 (IDEXO4, INFREQ4, INTRET4, RELRET4, CHRRET, IERR4)
        NEXSIDS = INTRET4
C       INFREQ4 = EXSSEL
C       CALL EXINQ4 (IDEXO4, INFREQ4, INTRET4, RELRET4, CHRRET, IERR4)
C       NEXSSEL = INTRET4
C       INFREQ4 = EXSSDF
C       CALL EXINQ4 (IDEXO4, INFREQ4, INTRET4, RELRET4, CHRRET, IERR4)
C       NEXSSDF = INTRET4
        NEXSSEL = 0
        DO I = 1, NEXSIDS
                NEXSSEL = NEXSSEL + NSESS(I)
        END DO
        CALL I8I4 (NEXSIDS, IDESSS, PIDESSS4)
        CALL I8I4 (NEXSIDS, NSESS, PNSESS4)
        CALL I8I4 (NEXSIDS, NDESS, PNDESS4)
        CALL I8I4 (NEXSIDS, IXEESS, PIXEESS4)
        CALL I8I4 (NEXSIDS, IXDESS, PIXDESS4)
        CALL I8I4 (NEXSSEL, LTEESS, PLTEESS4)
        CALL I8I4 (NEXSSEL, LTSESS, PLTSESS4)
        CALL EXPCSS4 (IDEXO4, IDESSS4, NSESS4, NDESS4, IXEESS4,
     &          IXDESS4, LTEESS4, LTSESS4, FACESS, IERR4)
        CALL I4FREE (PIDESSS4)
        CALL I4FREE (PNSESS4)
        CALL I4FREE (PNDESS4)
        CALL I4FREE (PIXEESS4)
        CALL I4FREE (PIXDESS4)
        CALL I4FREE (PLTEESS4)
        CALL I4FREE (PLTSESS4)
        IERR = IERR4
        END

        SUBROUTINE EXPEAT (IDEXO, IDELB, ATRIB, IERR)
C
C       WRITE ELEMENT BLOCK ATTRIBUTES
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDELB           ! (R)
        REAL ATRIB(*)           ! (NUMATR,NUMELB) ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IDELB4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        IDELB4 = IDELB

        CALL EXPEAT4 (IDEXO4, IDELB4, ATRIB, IERR4)

        IERR = IERR4
        END

        SUBROUTINE EXPEAN (IDEXO, IDELB, NATTR, NAMES, IERR)
C
C       WRITE ELEMENT BLOCK ATTRIBUTE NAMES
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDELB           ! (R)
        INTEGER NATTR           ! (R)
        CHARACTER* (*) NAMES (*) ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IDELB4        ! (R)
        INTEGER*4 NATTR4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        IDELB4 = IDELB
        NATTR4 = NATTR

        CALL EXPEAN4 (IDEXO4, IDELB4, NATTR4, NAMES, IERR4)

        IERR = IERR4
        END

        SUBROUTINE EXPCLB (IDEXO, IDELB, NAMELB, NUMELB, NUMLNK,
     &          NUMATR, MKMAP, IERR)
C
C       WRITE ELEMENT BLOCK PARAMETERS
C
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
        INTEGER*4 IDELB4(1)     ! (R)
        INTEGER*4 NUMELB4(1)    ! (R)
        INTEGER*4 NUMLNK4(1)    ! (R)
        INTEGER*4 NUMATR4(1)    ! (R)
        INTEGER*4 MKMAP4        ! (R)
        POINTER (PIDELB4,  IDELB4)
        POINTER (PNUMELB4, NUMELB4)
        POINTER (PNUMLNK4, NUMLNK4)
        POINTER (PNUMATR4, NUMATR4)
        INTEGER*4 IERR4         ! (W)

        INTEGER*4 INFREQ4
        INTEGER*4 INTRET4
        REAL*4 RELRET4
        CHARACTER* (4) CHRRET

        INTEGER*4 NUMBLK

        IDEXO4 = IDEXO

        INFREQ4 = EXELBL
        CALL EXINQ4 (IDEXO4, INFREQ4, INTRET4, RELRET4, CHRRET, IERR4)
        NUMBLK = INTRET4

        CALL I8I4 (NUMBLK, IDELB,  PIDELB4)
        CALL I8I4 (NUMBLK, NUMELB, PNUMELB4)
        CALL I8I4 (NUMBLK, NUMLNK, PNUMLNK4)
        CALL I8I4 (NUMBLK, NUMATR, PNUMATR4)
        MKMAP4  = MKMAP

        CALL EXPCLB4 (IDEXO4, IDELB4, NAMELB, NUMELB4, NUMLNK4,
     &          NUMATR4, MKMAP4, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPELB (IDEXO, IDELB, NAMELB, NUMELB, NUMLNK,
     &          NUMATR, IERR)
C
C       WRITE ELEMENT BLOCK PARAMETERS
C
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
        INTEGER*4 IDELB4        ! (R)
        INTEGER*4 NUMELB4       ! (R)
        INTEGER*4 NUMLNK4       ! (R)
        INTEGER*4 NUMATR4       ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        IDELB4 = IDELB
        NUMELB4 = NUMELB
        NUMLNK4 = NUMLNK
        NUMATR4 = NUMATR
        CALL EXPELB4 (IDEXO4, IDELB4, NAMELB, NUMELB4, NUMLNK4,
     &          NUMATR4, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPELC (IDEXO, IDELB, LINK, IERR)
C
C       WRITE ELEMENT BLOCK CONNECTIVITY
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDELB           ! (R)
        INTEGER LINK (*)        ! (NUMLNK,NUMELB) ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IDELB4        ! (R)
        INTEGER*4 LINK4 (1)     ! (NUMLNK,NUMELB) ! (R)
        POINTER (PLINK4, LINK4)
        INTEGER*4 IERR4         ! (W)

        CHARACTER* (MXSTLN) NAMELB
        INTEGER*4 NUMELB4
        INTEGER*4 NUMLNK4
        INTEGER*4 NUMATR4

        IDEXO4 = IDEXO
        IDELB4 = IDELB
        CALL EXGELB4 (IDEXO4, IDELB4, NAMELB, NUMELB4, NUMLNK4,
     &          NUMATR4, IERR4)
        IF (NAMELB .EQ. 'NSIDED' .OR. NAMELB .EQ. 'nsided') THEN
          CALL I8I4 (NUMLNK4, LINK, PLINK4)
        ELSE
          CALL I8I4 (NUMLNK4*NUMELB4, LINK, PLINK4)
        ENDIF
        CALL EXPELC4 (IDEXO4, IDELB4, LINK4, IERR4)
        CALL I4FREE (PLINK4)
        IERR = IERR4
        END

        SUBROUTINE EXPECPP (IDEXO, ITYPE, IDELB, COUNTS, IERR)
C
C       WRITE ELEMENT BLOCK CONNECTIVITY
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER ITYPE           ! (R)
        INTEGER IDELB           ! (R)
        INTEGER COUNTS (*)      ! (NUMELB) ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 ITYPE4        ! (R)
        INTEGER*4 IDELB4        ! (R)
        INTEGER*4 COUNTS4 (1)   ! (NUMELB) ! (R)
        POINTER (PCOUNTS4, COUNTS4)
        INTEGER*4 IERR4         ! (W)

        CHARACTER* (MXSTLN) NAMELB
        INTEGER*4 NUMELB4
        INTEGER*4 NUMLNK4
        INTEGER*4 NUMATR4

        IDEXO4 = IDEXO
        ITYPE4 = ITYPE
        IDELB4 = IDELB
        CALL EXGELB4 (IDEXO4, IDELB4, NAMELB, NUMELB4, NUMLNK4,
     &          NUMATR4, IERR4)
        CALL I8I4 (NUMELB4, COUNTS, PCOUNTS4)
        CALL EXPECPP4 (IDEXO4, ITYPE4, IDELB4, COUNTS4, IERR4)
        CALL I4FREE (PCOUNTS4)
        IERR = IERR4
        END

        SUBROUTINE EXPENM (IDEXO, MAPEL, IERR)
C
C       WRITE ELEMENT NUMBER MAP
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER MAPEL (*)       ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 MAPEL4 (1)    ! (R)
        POINTER (PMAPEL4, MAPEL4)
        INTEGER*4 IERR4         ! (W)

        INTEGER*4 INFREQ4
        INTEGER*4 INTRET4
        REAL*4 RELRET4
        CHARACTER* (4) CHRRET

        INTEGER*4 NEXELEM

        IDEXO4 = IDEXO
        INFREQ4 = EXELEM
        CALL EXINQ4 (IDEXO4, INFREQ4, INTRET4, RELRET4, CHRRET, IERR4)
        NEXELEM = INTRET4
        CALL I8I4 (NEXELEM, MAPEL, PMAPEL4)
        CALL EXPENM4 (IDEXO4, MAPEL4, IERR4)
        CALL I4FREE (PMAPEL4)
        IERR = IERR4
        END

        SUBROUTINE EXPEM (IDEXO, MAPID, MAPEL, IERR)
C
C       WRITE ELEMENT NUMBER MAP
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER MAPID           ! (R)
        INTEGER MAPEL (*)       ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 MAPID4        ! (R)
        INTEGER*4 MAPEL4 (1)    ! (R)
        POINTER (PMAPEL4, MAPEL4)
        INTEGER*4 IERR4         ! (W)

        INTEGER*4 INFREQ4
        INTEGER*4 INTRET4
        REAL*4 RELRET4
        CHARACTER* (4) CHRRET

        INTEGER*4 NEXELEM

        IDEXO4 = IDEXO
        INFREQ4 = EXELEM
        MAPID4 = MAPID
        CALL EXINQ4 (IDEXO4, INFREQ4, INTRET4, RELRET4, CHRRET, IERR4)
        NEXELEM = INTRET4
        CALL I8I4 (NEXELEM, MAPEL, PMAPEL4)
        CALL EXPEM4 (IDEXO4, MAPID4, MAPEL4, IERR4)
        CALL I4FREE (PMAPEL4)
        IERR = IERR4
        END

        SUBROUTINE EXPPEM (IDEXO, MAPID, ISTART, ICOUNT, MAPEL, IERR)
C
C       WRITE ELEMENT NUMBER MAP
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER MAPID           ! (R)
        INTEGER ISTART          ! (R)
        INTEGER ICOUNT          ! (R)
        INTEGER MAPEL (*)       ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 MAPID4        ! (R)
        INTEGER*4 ISTART4       ! (R)
        INTEGER*4 ICOUNT4       ! (R)
        INTEGER*4 MAPEL4 (1)    ! (R)
        POINTER (PMAPEL4, MAPEL4)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        MAPID4 = MAPID
        ISTART4 = ISTART
        ICOUNT4 = ICOUNT
        CALL I8I4 (ICOUNT4, MAPEL, PMAPEL4)
        CALL EXPPEM4 (IDEXO4, MAPID4, ISTART4, ICOUNT4, MAPEL4, IERR4)
        CALL I4FREE (PMAPEL4)
        IERR = IERR4
        END

        SUBROUTINE EXPEV (IDEXO, ISTEP, IXELEV, IDELB, NUMELB,
     &          VALEV, IERR)
C
C       WRITE ELEMENT VARIABLE VALUES AT A TIME STEP
C
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
        INTEGER*4 IDELB4        ! (R)
        INTEGER*4 NUMELB4       ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4  = IDEXO
        ISTEP4  = ISTEP
        IXELEV4 = IXELEV
        IDELB4  = IDELB
        NUMELB4 = NUMELB
        CALL EXPEV4 (IDEXO4, ISTEP4, IXELEV4, IDELB4, NUMELB4,
     &          VALEV, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPNSV (IDEXO, ISTEP, IXV, ID, NUM, VAL, IERR)
C
C       WRITE NODESET VARIABLE VALUES AT A TIME STEP
C
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
        INTEGER*4 ID4   ! (R)
        INTEGER*4 NUM4  ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4  = IDEXO
        ISTEP4  = ISTEP
        IXV4 = IXV
        ID4  = ID
        NUM4 = NUM
        CALL EXPNSV4 (IDEXO4, ISTEP4, IXV4, ID4, NUM4, VAL, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPSSV (IDEXO, ISTEP, IXV, ID, NUM, VAL, IERR)
C
C       WRITE SIDESET VARIABLE VALUES AT A TIME STEP
C
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
        INTEGER*4 ID4   ! (R)
        INTEGER*4 NUM4  ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4  = IDEXO
        ISTEP4  = ISTEP
        IXV4 = IXV
        ID4  = ID
        NUM4 = NUM
        CALL EXPSSV4 (IDEXO4, ISTEP4, IXV4, ID4, NUM4, VAL, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPGV (IDEXO, ISTEP, NVARGL, VALGV, IERR)
C
C       WRITE GLOBAL VARIABLE VALUES AT A TIME STEP
C
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
C
C       WRITE INFORMATION RECORDS
C
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
C
C       WRITE INITIALIZATION PARAMETERS
C
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
        INTEGER*4 NDIM4         ! (R)
        INTEGER*4 NUMNP4        ! (R)
        INTEGER*4 NUMEL4        ! (R)
        INTEGER*4 NELBLK4       ! (R)
        INTEGER*4 NUMNPS4       ! (R)
        INTEGER*4 NUMESS4       ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        NDIM4 = NDIM
        NUMNP4 = NUMNP
        NUMEL4 = NUMEL
        NELBLK4 = NELBLK
        NUMNPS4 = NUMNPS
        NUMESS4 = NUMESS
        CALL EXPINI4 (IDEXO4, TITLE, NDIM4, NUMNP4, NUMEL4, NELBLK4,
     &          NUMNPS4, NUMESS4, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPMAP (IDEXO, MAPEL, IERR)
C
C       WRITE ELEMENT ORDER MAP
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER MAPEL (*)       ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 MAPEL4 (1)    ! (R)
        POINTER (PMAPEL4, MAPEL4)
        INTEGER*4 IERR4         ! (W)

        INTEGER*4 INFREQ4
        INTEGER*4 INTRET4
        REAL*4 RELRET4
        CHARACTER* (4) CHRRET

        INTEGER*4 NEXELEM

        IDEXO4 = IDEXO
        INFREQ4 = EXELEM
        CALL EXINQ4 (IDEXO4, INFREQ4, INTRET4, RELRET4, CHRRET, IERR4)
        NEXELEM = INTRET4
        CALL I8I4 (NEXELEM, MAPEL, PMAPEL4)
        CALL EXPMAP4 (IDEXO4, MAPEL4, IERR4)
        CALL I4FREE (PMAPEL4)
        IERR = IERR4
        END

        SUBROUTINE EXPNNM (IDEXO, MAPNOD, IERR)
C
C       WRITE NODE NUMBER MAP
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER MAPNOD (*)      ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 MAPNOD4 (1)   ! (R)
        POINTER (PMAPNOD4, MAPNOD4)
        INTEGER*4 IERR4         ! (W)

        INTEGER*4 INFREQ4
        INTEGER*4 INTRET4
        REAL*4 RELRET4
        CHARACTER* (4) CHRRET

        INTEGER*4 NEXNODE

        IDEXO4 = IDEXO
        INFREQ4 = EXNODE
        CALL EXINQ4 (IDEXO4, INFREQ4, INTRET4, RELRET4, CHRRET, IERR4)
        NEXNODE = INTRET4
        CALL I8I4 (NEXNODE, MAPNOD, PMAPNOD4)
        CALL EXPNNM4 (IDEXO4, MAPNOD4, IERR4)
        CALL I4FREE (PMAPNOD4)
        IERR = IERR4
        END

        SUBROUTINE EXPNM (IDEXO, MAPID, MAPNOD, IERR)
C
C       WRITE NODE NUMBER MAP
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER MAPID           ! (R)
        INTEGER MAPNOD (*)      ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 MAPID4        ! (R)
        INTEGER*4 MAPNOD4 (1)   ! (R)
        POINTER (PMAPNOD4, MAPNOD4)
        INTEGER*4 IERR4         ! (W)

        INTEGER*4 INFREQ4
        INTEGER*4 INTRET4
        REAL*4 RELRET4
        CHARACTER* (4) CHRRET

        INTEGER*4 NEXNODE

        IDEXO4 = IDEXO
        INFREQ4 = EXNODE
        MAPID4 = MAPID
        CALL EXINQ4 (IDEXO4, INFREQ4, INTRET4, RELRET4, CHRRET, IERR4)
        NEXNODE = INTRET4
        CALL I8I4 (NEXNODE, MAPNOD, PMAPNOD4)
        CALL EXPNM4 (IDEXO4, MAPID4, MAPNOD4, IERR4)
        CALL I4FREE (PMAPNOD4)
        IERR = IERR4
        END

        SUBROUTINE EXPNP (IDEXO, IDNPS, NNNPS, NDNPS, IERR)
C
C       WRITE NODE SET PARAMETERS
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDNPS           ! (R)
        INTEGER NNNPS           ! (R)
        INTEGER NDNPS           ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IDNPS4        ! (R)
        INTEGER*4 NNNPS4        ! (R)
        INTEGER*4 NDNPS4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        IDNPS4 = IDNPS
        NNNPS4 = NNNPS
        NDNPS4 = NDNPS
        CALL EXPNP4 (IDEXO4, IDNPS4, NNNPS4, NDNPS4, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPNS (IDEXO, IDNPS, LTNNPS, IERR)
C
C       WRITE NODE SET
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDNPS           ! (R)
        INTEGER LTNNPS (*)      ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IDNPS4        ! (R)
        INTEGER*4 LTNNPS4 (1)   ! (R)
        POINTER (PLTNNPS4, LTNNPS4)
        INTEGER*4 IERR4         ! (W)

        INTEGER*4 NNNPS4
        INTEGER*4 NDNPS4

        IDEXO4 = IDEXO
        IDNPS4 = IDNPS
        CALL EXGNP4 (IDEXO4, IDNPS4, NNNPS4, NDNPS4, IERR4)
        CALL I8I4 (NNNPS4, LTNNPS, PLTNNPS4)
        CALL EXPNS4 (IDEXO4, IDNPS4, LTNNPS4, IERR4)
        CALL I4FREE (PLTNNPS4)
        IERR = IERR4
        END

        SUBROUTINE EXPNSD (IDEXO, IDNPS, FACNPS, IERR)
C
C       WRITE NODE SET DISTRIBUTION FACTORS
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDNPS           ! (R)
        REAL FACNPS (*)         ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IDNPS4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        IDNPS4 = IDNPS
        CALL EXPNSD4 (IDEXO4, IDNPS4, FACNPS, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPNV (IDEXO, ISTEP, IXNODV, NUMNP, VALNV, IERR)
C
C       WRITE NODAL VARIABLE VALUES AT A TIME STEP
C
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
        INTEGER*4 NUMNP4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4  = IDEXO
        ISTEP4  = ISTEP
        IXNODV4 = IXNODV
        NUMNP4  = NUMNP
        CALL EXPNV4 (IDEXO4, ISTEP4, IXNODV4, NUMNP4, VALNV, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPP (IDEXO, ITYPE, ID, NAMEPR, IVAL, IERR)
C
C       WRITE OBJECT PROPERTY
C
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
        INTEGER*4 ID4           ! (R)
        INTEGER*4 IVAL4         ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        ITYPE4 = ITYPE
        ID4 = ID
        IVAL4 = IVAL
        CALL EXPP4 (IDEXO4, ITYPE4, ID4, NAMEPR, IVAL4, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPPA (IDEXO, ITYPE, NAMEPR, IVAL, IERR)
C
C       WRITE OBJECT PROPERTY ARRAY
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER ITYPE           ! (R)
        CHARACTER* (*) NAMEPR   ! (R)
        INTEGER IVAL (*)        ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 ITYPE4        ! (R)
        INTEGER*4 IVAL4 (1)     ! (R)
        POINTER (PIVAL4, IVAL4)
        INTEGER*4 IERR4         ! (W)

        INTEGER*4 INFREQ4
        INTEGER*4 INTRET4
        REAL*4 RELRET4
        CHARACTER* (4) CHRRET

        INTEGER*4 NUMOBJ

        IDEXO4 = IDEXO
        ITYPE4 = ITYPE
        IF (ITYPE4 .EQ. EXEBLK) THEN
                INFREQ4 = EXELBL
        ELSE IF (ITYPE4 .EQ. EXNSET) THEN
                INFREQ4 = EXNODS
        ELSE IF (ITYPE4 .EQ. EXSSET) THEN
                INFREQ4 = EXSIDS
        ELSE IF (ITYPE4 .EQ. EXEMAP) THEN
                INFREQ4 = EXNEM
        ELSE IF (ITYPE4 .EQ. EXNMAP) THEN
                INFREQ4 = EXNNM
        END IF
        CALL EXINQ4 (IDEXO4, INFREQ4, INTRET4, RELRET4, CHRRET, IERR4)
        NUMOBJ = INTRET4
        CALL I8I4 (NUMOBJ, IVAL, PIVAL4)
        CALL EXPPA4 (IDEXO4, ITYPE4, NAMEPR, IVAL4, IERR4)
        CALL I4FREE (PIVAL4)
        IERR = IERR4
        END

        SUBROUTINE EXPPN (IDEXO, ITYPE, NPROPS, NAMEPR, IERR)
C
C       WRITE PROPERTY ARRAY NAMES
C
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
C
C       WRITE QA RECORDS
C
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
C
C       WRITE SIDE SET PARAMETERS
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDESS           ! (R)
        INTEGER NSESS           ! (R)
        INTEGER NDESS           ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IDESS4        ! (R)
        INTEGER*4 NSESS4        ! (R)
        INTEGER*4 NDESS4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        IDESS4 = IDESS
        NSESS4 = NSESS
        NDESS4 = NDESS
        CALL EXPSP4 (IDEXO4, IDESS4, NSESS4, NDESS4, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPSS (IDEXO, IDESS, LTEESS, LTSESS, IERR)
C
C       WRITE SIDE SET
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDESS           ! (R)
        INTEGER LTEESS (*)      ! (R)
        INTEGER LTSESS (*)      ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IDESS4        ! (R)
        INTEGER*4 LTEESS4 (1)   ! (R)
        INTEGER*4 LTSESS4 (1)   ! (R)
        POINTER (PLTEESS4, LTEESS4)
        POINTER (PLTSESS4, LTSESS4)
        INTEGER*4 IERR4         ! (W)

        INTEGER*4 NSESS4
        INTEGER*4 NDESS4

        IDEXO4 = IDEXO
        IDESS4 = IDESS
        CALL EXGSP4 (IDEXO4, IDESS4, NSESS4, NDESS4, IERR4)
        CALL I8I4 (NSESS4, LTEESS, PLTEESS4)
        CALL I8I4 (NSESS4, LTSESS, PLTSESS4)
        CALL EXPSS4 (IDEXO4, IDESS4, LTEESS4, LTSESS4, IERR4)
        CALL I4FREE (PLTEESS4)
        CALL I4FREE (PLTSESS4)
        IERR = IERR4
        END

        SUBROUTINE EXPSSD (IDEXO, IDESS, FACESS, IERR)
C
C       WRITE SIDE SET DISTRIBUTION FACTORS
C
        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER IDEXO           ! (R)
        INTEGER IDESS           ! (R)
        REAL FACESS (*)         ! (R)
        INTEGER IERR            ! (W)

        INTEGER*4 IDEXO4        ! (R)
        INTEGER*4 IDESS4        ! (R)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        IDESS4 = IDESS
        CALL EXPSSD4 (IDEXO4, IDESS4, FACESS, IERR4)
        IERR = IERR4
        END

        SUBROUTINE EXPTIM (IDEXO, NSTEP, TIME, IERR)
C
C       WRITE TIME VALUE FOR A TIME STEP
C
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
C
C       WRITE RESULTS VARIABLE NAMES
C
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
C
C       WRITE RESULTS VARIABLE PARAMETERS
C
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
C
C       WRITE ELEMENT VARIABLE TRUTH TABLE
C
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
        INTEGER*4 ISEVOK4 (1)   ! (NVAREL,NELBLK) (R)
        POINTER (PISEVOK4, ISEVOK4)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        NELBLK4 = NELBLK
        NVAREL4 = NVAREL
        CALL I8I4 (NVAREL4*NELBLK4, ISEVOK, PISEVOK4)
        CALL EXPVTT4 (IDEXO4, NELBLK4, NVAREL4, ISEVOK4, IERR4)
        CALL I4FREE (PISEVOK4)
        IERR = IERR4
        END

        SUBROUTINE EXPNSTT (IDEXO, NBLK, NVAR, ISVOK, IERR)
C
C       WRITE NODESET VARIABLE TRUTH TABLE
C
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
        INTEGER*4 ISVOK4 (1)    ! (NVAR,NBLK) (R)
        POINTER (PISVOK4, ISVOK4)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        NBLK4 = NBLK
        NVAR4 = NVAR
        CALL I8I4 (NVAR4*NBLK4, ISVOK, PISVOK4)
        CALL EXPNSTT4 (IDEXO4, NBLK4, NVAR4, ISVOK4, IERR4)
        CALL I4FREE (PISVOK4)
        IERR = IERR4
        END

        SUBROUTINE EXPSSTT (IDEXO, NBLK, NVAR, ISVOK, IERR)
C
C       WRITE SIDESET VARIABLE TRUTH TABLE
C
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
        INTEGER*4 ISVOK4 (1)    ! (NVAR,NBLK) (R)
        POINTER (PISVOK4, ISVOK4)
        INTEGER*4 IERR4         ! (W)

        IDEXO4 = IDEXO
        NBLK4 = NBLK
        NVAR4 = NVAR
        CALL I8I4 (NVAR4*NBLK4, ISVOK, PISVOK4)
        CALL EXPSSTT4 (IDEXO4, NBLK4, NVAR4, ISVOK4, IERR4)
        CALL I4FREE (PISVOK4)
        IERR = IERR4
        END

        SUBROUTINE EXGNAMS (IDEXO, ITYPE, NVAR, NAMES, IERR)
C
C       READ OBJECT NAMES
C
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

        SUBROUTINE EXPNAMS (IDEXO, ITYPE, NVAR, NAMES, IERR)
C
C       WRITE OBJECT NAMES
C
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

        SUBROUTINE EXUPDA (IDEXO, IERR)
C
C       UPDATE EXODUS II FILE
C
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
C
C       COPY EXODUS II FILE
C
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


        SUBROUTINE I8I4 (N, I8, PI4)
C
C       CREATE I4 ARRAY AND COPY I8 ARRAY CONTENTS TO I4
C
        IMPLICIT NONE
        INTEGER*4 N
        INTEGER I8 (*)
        INTEGER*4 I4 (1)
        POINTER (PI4, I4)
        INTEGER*4 I

        CALL I4ALLOC (N, PI4)
        DO I = 1, N
                I4 (I) = I8 (I)
        END DO
        END

        SUBROUTINE I4I8 (N, PI4, I8)
C
C       COPY I4 ARRAY CONTENTS TO I8 AND FREE I4 ARRAY
C
        IMPLICIT NONE
        INTEGER*4 N
        INTEGER*4 I4 (1)
        POINTER (PI4, I4)
        INTEGER I8 (*)
        INTEGER*4 I

        DO I = 1, N
                I8 (I) = I4 (I)
        END DO
        CALL I4FREE (PI4)
        END

        SUBROUTINE I4ALLOC (N, PI4)
C
C       ALLOCATE DYNAMIC I4 ARRAY N ELEMENTS IN SIZE
C
        IMPLICIT NONE
        INTEGER*4 N
        INTEGER*4 I4 (1)
        POINTER (PI4, I4)
        INTEGER MALLOC
        INTEGER*4 NB

        PI4 = 0
        IF (N .EQ. 0) RETURN
        NB = N * 4
        PI4 = MALLOC (NB)
        END

        SUBROUTINE I4FREE (PI4)
C
C       FREE DYNAMIC MEMORY ASSOCIATED WITH I4 ARRAY REFERENCED BY PI4
C
        IMPLICIT NONE
        INTEGER*4 I4 (1)
        POINTER (PI4, I4)

        IF (PI4 .EQ. 0) RETURN
        CALL FREE (PI4)
        PI4 = 0
        END
#endif
