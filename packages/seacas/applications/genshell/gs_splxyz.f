C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C     -*- Mode: fortran -*-
C=======================================================================
      SUBROUTINE SPLXYZ (XN, YN, XN3, YN3, ZN3, ATRIB,
     &     NSPL1, RSA, ZSA, ZS2A, DISTA, SCRA,
     &     SLLFT, SLRGT, RDTHET, SWEEP, NOSCAL )
C=======================================================================

C     --*** SPLXYZ *** (GEN3D) Calculate 3D coordinates for Double-surface
C     Spline projection
C     --   Written by Greg Sjaardema - 01/08/90
C     --
C     --SPLXYZ calculates the coordinate array for the 3D spline translations.
C     --
C     --Parameters:
C     --   XN, YN - IN - the 2D coordinates, destroyed
C     --   XN3, YN3, ZN3 - OUT - the 3D coordinates
C     --   IXNP - IN - the new index for each node
C     --   NRNP - IN - the number of new nodes generated for each node
C     --   ZCORD - SCRATCH - size = NNREPL, holds z coordinate for transformations
C     --
C     --Common Variables:
C     --   Uses NDIM, NUMNP of /DBNUMS/
C     --   Uses NDIM3, NUMNP3 of /DBNUM3/
C     --   Uses DOTRAN, NNREPL, DIM3, NRTRAN, D3TRAN, ZGRAD,
C     --      CENTER, NUMCOL, NUMROW of /PARAMS/

      INCLUDE 'gs_dbnums.blk'
      INCLUDE 'gs_dbnum3.blk'
      INCLUDE 'gs_params.blk'

      PARAMETER (BINGO = 1.0E38)
      PARAMETER (TOLER = 1.0E-8)

      INTEGER XSWEEP, YSWEEP, SPHERI, SWEEP
      PARAMETER (XSWEEP = 10)
      PARAMETER (YSWEEP =  1)
      PARAMETER (SPHERI = 11)

      REAL XN(NUMNP), YN(NUMNP),
     &     XN3(NUMNP3), YN3(NUMNP3), ZN3(NUMNP3)
      REAL ATRIB(NUMEL)
      REAL RSA(NSPL1), ZSA(NSPL1), ZS2A(NSPL1), DISTA(NSPL1),
     &     SCRA(NSPL1)
C ... NOTE: Only sllft(1) and slrgt(1) are used, (2) left for comp. with gen3d
      REAL SLLFT(2), SLRGT(2)

      LOGICAL RDTHET, NOSCAL

      rad = 0.0
      xt = 0.0
      yt = 0.0

      PI = ATAN2(0.0, -1.0)

C ... Check for valid options...
      IF (RDTHET .AND. NOSCAL) THEN
        CALL PRTERR('ERROR', 'Cannot use NOSCALE with ANGULAR spline')
        RETURN
      END IF

C     ... Convert angles from degrees to radians
      IF (RDTHET) THEN
        DO 10 IPT = 1, NSPL1
            RSA(IPT) = RSA(IPT) * PI / 180.0
   10    CONTINUE
      END IF

      CALL SPLINE (RSA, ZSA, NSPL1, SLLFT(1), SLRGT(1), ZS2A, SCRA)

C     ... Calculate approximate distance along curve.  Distance is computed
C     as the straight line distance of the segments.  Accurate for
C     smooth curves; bad for non-smooth curves.
C     ... If 'NOSCALE' is set, then the spline must be linear and
C     we use the input RSA values for the distance.

      IF (NOSCAL) THEN
        DO 20 IPT = 1, NSPL1
          DISTA(IPT) = RSA(IPT)
 20     CONTINUE
      ELSE
        IF (RDTHET) THEN
          DISTA(1) = ZSA(1) * SIN(RSA(1))
          DO 30 IPT = 2, NSPL1
            DISTA(IPT) = DISTA(IPT-1)  + SQRT(
     &        (ZSA(IPT)   * SIN(RSA(IPT)) -
     &        ZSA(IPT-1) * SIN(RSA(IPT-1)))**2 +
     &        (ZSA(IPT)   * COS(RSA(IPT)) -
     &        ZSA(IPT-1) * COS(RSA(IPT-1)))**2)
 30       CONTINUE
        ELSE
          DISTA(1) = RSA(1)
          DO 50 IPT = 2, NSPL1
            DISTA(IPT) = SQRT( (RSA(IPT-1) - RSA(IPT))**2 +
     &        (ZSA(IPT-1) - ZSA(IPT))**2 ) + DISTA(IPT-1)
 50       CONTINUE
        END IF
      END IF
C...  Determine maximum radius of mesh

      RMAX = 0.0
      RMIN = BINGO
      IF (SWEEP .EQ. SPHERI) THEN
      DO 80 INP = 1, NUMNP
         DX = XN(INP)
         DY = YN(INP)
         RMIN = MIN (RMIN, (DX**2 + DY**2))
         RMAX = MAX (RMAX, (DX**2 + DY**2))
   80 CONTINUE
      RMIN = SQRT(RMIN)
      RMAX = SQRT(RMAX)
      ELSE IF (SWEEP .EQ. XSWEEP) THEN
         DO 82 INP = 1, NUMNP
            DY = YN(INP)
            RMIN = MIN (RMIN, DY)
            RMAX = MAX (RMAX, DY)
 82      CONTINUE
      ELSE IF (SWEEP .EQ. YSWEEP) THEN
         DO 84 INP = 1, NUMNP
            DX = XN(INP)
            RMIN = MIN (RMIN, DX)
            RMAX = MAX (RMAX, DX)
 84      CONTINUE
      END IF

C... If there is not a node at 0,0 - RMIN will not equal zero, therefore,
C     we assume that if the spline starts at zero, the mesh also starts
C     at zero.
      IF (DISTA(1) .EQ. 0.0) THEN
         RMIN = 0.0
      END IF

      IF (NOSCAL) THEN
        PROPA = 1.0
      ELSE
        PROPA = (DISTA(NSPL1) - DISTA(1)) / (RMAX - RMIN)
      END IF

C     ... Echo spline data

      WRITE (*,*) ' '
      WRITE (*,70) 'top spline', DISTA(NSPL1), PROPA
   70 FORMAT (' Total length of ',A,T32,' = ',1PE12.5,
     & ', Proportion = ',1PE12.5)

      KLOA = 1
      DO 100 INP = 1, NUMNP
         DX = XN(INP)
         DY = YN(INP)
         IF (SWEEP .EQ. SPHERI) THEN
            RAD = SQRT (DX**2 + DY**2)
         ELSE IF (SWEEP .EQ. XSWEEP) THEN
            RAD = DY
         ELSE IF (SWEEP .EQ. YSWEEP) THEN
            RAD = DX
         END IF
         RADA = DISTA(1) + PROPA * (RAD - RMIN)

C     ... Determine segment containing point and proportion within segment

         CALL HUNT (DISTA, NSPL1, RADA, KLOA)
         KLOA = MIN(NSPL1-1, MAX(1, KLOA))

C     ... We want to map the X-Y plane onto the spline surface without
C     stretching.  Therefore, calculate a new radius RNEW which is
C     the distance from the Z axis to the point on the spline such
C     that the arc length from R=0 to R=RAD is the distance RNEW.
C     The new X and Y coordinates corresponding to RNEW are calculated
C     as the proportion of RNEW / RAD.  FIXR is added to RAD in case
C     RAD is 0; FIXR = 1 iff RAD = 0, else FIXR = 0

C --- On IEEE Machines, at least on Ultrix 4.0 IEEE, the SIGN
C     trick does not work since it defines a signed 0.0

         IF (RAD .EQ. 0.0) THEN
            FIXR = 1.0
         ELSE
            FIXR = 0.0
         END IF
C         FIXR  = SIGN(0.5, RAD) + SIGN(0.5, -RAD)

         PROP  = (RADA - DISTA(KLOA)) / (DISTA(KLOA+1) - DISTA(KLOA))
         RNEWA = RSA(KLOA) + PROP * (RSA(KLOA+1) - RSA(KLOA))

         H  = RSA(KLOA+1) - RSA(KLOA)
         A  = (RSA(KLOA+1)-RNEWA) / H
         B  = (RNEWA-RSA(KLOA)) / H

         ZT = A * ZSA(KLOA) + B * ZSA(KLOA+1) +
     *     ((A**3-A) * ZS2A(KLOA)+(B**3-B) * ZS2A(KLOA+1)) * (H**2)/6.

         IF (RDTHET) THEN
            RSAV  = RNEWA
            ZSAV  = ZT
            ZT    = ZSAV * COS(RSAV)
            RNEWA = ZSAV * SIN(RSAV)
        END IF

         IF (SWEEP .EQ. SPHERI) THEN
C ... Spherical Sweep of Spline Surface
            XT = RNEWA * DX / (RAD + FIXR)
            YT = RNEWA * DY / (RAD + FIXR)

         ELSE IF (SWEEP .EQ. XSWEEP) THEN
C ... Sweep spline along the X axis
            XT = DX
            YT = RNEWA * DY / (RAD + FIXR)

         ELSE IF (SWEEP .EQ. YSWEEP) THEN
C ... Sweep spline along the X axis
            XT = RNEWA * DX / (RAD + FIXR)
            YT = DY

         END IF

C$$$$         write (*,*) inp, xn(inp), yn(inp), xt, yt, zt
         ZN3(INP) = ZT
         YN3(INP) = YT
         XN3(INP) = XT
  100 CONTINUE

C ... Now do the attributes for all of the elements

      DO 370 IEL = 1, NUMEL
         ATRIB(IEL) = DIM3
 370  CONTINUE

      RETURN
      END
