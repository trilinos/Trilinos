C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SPLXYZ (XN, YN, XN3, YN3, ZN3, IXNP, NRNP, ZCORD,
     &     NSPL1, NSPL2, RSA, ZSA, ZS2A, DISTA, SCRA,
     &     RSB, ZSB, ZS2B, DISTB, SCRB, SLLFT, SLRGT, RDTHET,
     &     SWEEP, NOSCAL)
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

      INCLUDE 'g3_dbnums.blk'
      INCLUDE 'g3_dbnum3.blk'
      INCLUDE 'g3_params.blk'

      PARAMETER (BINGO = 1.0E38)
      PARAMETER (TOLER = 1.0E-8)

      INTEGER XSWEEP, YSWEEP, SPHERI, SWEEP
      PARAMETER (XSWEEP = 10)
      PARAMETER (YSWEEP =  1)
      PARAMETER (SPHERI = 11)

      REAL XN(NUMNP), YN(NUMNP),
     &     XN3(NUMNP3), YN3(NUMNP3), ZN3(NUMNP3)
      INTEGER IXNP(*), NRNP(*)
      REAL ZCORD(NNREPL)
      REAL RSA(NSPL1), ZSA(NSPL1), ZS2A(NSPL1), DISTA(NSPL1),
     &     SCRA(NSPL1)
      REAL RSB(NSPL2), ZSB(NSPL2), ZS2B(NSPL2), DISTB(NSPL2),
     &     SCRB(NSPL2)
      REAL SLLFT(2), SLRGT(2)

      LOGICAL RDTHET, NOSCAL

      PI = ATAN2(0.0, -1.0)
      xb = 0.0
      xt = 0.0
      yb = 0.0
      yt = 0.0

C ... Check for valid options...
      IF (RDTHET .AND. NOSCAL) THEN
        CALL PRTERR('ERROR', 'Cannot use NOSCALE with ANGULAR spline')
        RETURN
      END IF

C     ... CALCULATE THE THICKNESS INCREMENT FOR EACH TRANSLATION
      IBLK = 0
      ZTOT = 0.0
    1 CONTINUE
      IBLK = IBLK + 1
      IF (NRTRAN(IBLK) .GT. 0) THEN
         ZTOT = ZTOT + D3TRAN(IBLK)
         IF (IBLK .LT. MAXINT) GO TO 1
      END IF

      NXTNR = 1
      IBLK = 0
      ZEND = 0.0

    2 CONTINUE
      IBLK = IBLK + 1
      IF (NRTRAN(IBLK) .GT. 0) THEN
         ZBEG = ZEND
         ZEND = ZBEG + D3TRAN(IBLK)
         CALL INIGRD (ZBEG/ZTOT, ZEND/ZTOT, ZGRAD(IBLK),
     *        NRTRAN(IBLK), NRTRAN(IBLK)+1, ZCORD(NXTNR) )
         NXTNR = NXTNR + NRTRAN(IBLK)
         IF (IBLK .LT. MAXINT) GO TO 2
      END IF

C     ... Convert angles from degrees to radians

      IF (RDTHET) THEN
         DO 10 IPT = 1, NSPL1
            RSA(IPT) = RSA(IPT) * PI / 180.0
   10    CONTINUE
         DO 20 IPT = 1, NSPL2
            RSB(IPT) = RSB(IPT) * PI / 180.0
   20    CONTINUE

      END IF

      CALL SPLINE (RSA, ZSA, NSPL1, SLLFT(1), SLRGT(1), ZS2A, SCRA)
      CALL SPLINE (RSB, ZSB, NSPL2, SLLFT(2), SLRGT(2), ZS2B, SCRB)

C     ... Calculate approximate distance along curve.  Distance is computed
C     as the straight line distance of the segments.  Accurate for
C     smooth curves; bad for non-smooth curves.
C     ... If 'NOSCALE' is set, then the spline must be linear and
C     we use the input RSA values for the distance.

      IF (NOSCAL) THEN
        DO 22 IPT = 1, NSPL1
          DISTA(IPT) = RSA(IPT)
 22     CONTINUE
        DO 24 IPT = 1, NSPL2
          DISTB(IPT) = RSB(IPT)
 24     CONTINUE
      ELSE
        IF (RDTHET) THEN
          DISTA(1) = ZSA(1) * SIN(RSA(1))
          DISTB(1) = ZSB(1) * SIN(RSB(1))
          DO 30 IPT = 2, NSPL1
            DISTA(IPT) = DISTA(IPT-1)  + SQRT(
     &        (ZSA(IPT)   * SIN(RSA(IPT)) -
     &        ZSA(IPT-1) * SIN(RSA(IPT-1)))**2 +
     &        (ZSA(IPT)   * COS(RSA(IPT)) -
     &        ZSA(IPT-1) * COS(RSA(IPT-1)))**2)
 30       CONTINUE
          DO 40 IPT = 2, NSPL2
            DISTB(IPT) = DISTB(IPT-1)  + SQRT(
     &        (ZSB(IPT)   * SIN(RSB(IPT)) -
     &        ZSB(IPT-1) * SIN(RSB(IPT-1)))**2 +
     &        (ZSB(IPT)   * COS(RSB(IPT)) -
     &        ZSB(IPT-1) * COS(RSB(IPT-1)))**2)
 40       CONTINUE
        ELSE
          DISTA(1) = RSA(1)
          DO 50 IPT = 2, NSPL1
            DISTA(IPT) = SQRT( (RSA(IPT-1) - RSA(IPT))**2 +
     &        (ZSA(IPT-1) - ZSA(IPT))**2 ) + DISTA(IPT-1)
 50       CONTINUE
          DISTB(1) = RSB(1)
          DO 60 IPT = 2, NSPL2
            DISTB(IPT) = SQRT( (RSB(IPT-1) - RSB(IPT))**2 +
     &        (ZSB(IPT-1) - ZSB(IPT))**2 ) + DISTB(IPT-1)
 60       CONTINUE
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
 80      CONTINUE
         RMAX = SQRT(RMAX)
         RMIN = SQRT(RMIN)
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
      IF (DISTA(1) .EQ. 0.0 .OR. DISTB(1) .EQ. 0.0) THEN
         RMIN = 0.0
      END IF

      IF (NOSCAL) THEN
        PROPA = 1.0
        PROPB = 1.0
      ELSE
        PROPA = (DISTA(NSPL1) - DISTA(1)) / (RMAX - RMIN)
        PROPB = (DISTB(NSPL2) - DISTB(1)) / (RMAX - RMIN)
      END IF

C     ... Echo spline data
      WRITE (*,*) ' '
      WRITE (*,70) 'top spline', DISTA(NSPL1), PROPA
      WRITE (*,70) 'bottom spline',DISTB(NSPL2), PROPB
   70 FORMAT (' Total length of ',A,T32,' = ',1PE12.5,
     &     ', Proportion = ',1PE12.5)

      KLOA = 1
      KLOB = 1
      DO 100 INP = 1, NUMNP
         JNP0 = IXNP(INP) - 1
         DX = XN(INP)
         DY = YN(INP)
         RAD = 0.0
         IF (SWEEP .EQ. SPHERI) THEN
            RAD = SQRT (DX**2 + DY**2)
         ELSE IF (SWEEP .EQ. XSWEEP) THEN
            RAD = DY
         ELSE IF (SWEEP .EQ. YSWEEP) THEN
            RAD = DX
         END IF

         RADA = DISTA(1) + PROPA * (RAD - RMIN)
         RADB = DISTB(1) + PROPB * (RAD - RMIN)

C     ... Determine segment containing point and proportion within segment
C     Assume last segment was close to this segment and start search
C     there.

         CALL HUNT (DISTA, NSPL1, RADA, KLOA)
         KLOA = MIN(NSPL1-1, MAX(1, KLOA))

         CALL HUNT (DISTB, NSPL2, RADB, KLOB)
         KLOB = MIN(NSPL2-1, MAX(1, KLOB))

C     ... We want to map the X-Y plane onto the spline surface without
C     stretching.  Therefore, calculate a new radius RNEW which is
C     the distance from the Z axis to the point on the spline such
C     that the arc length from R=0 to R=RAD is the distance RNEW.
C     The new X and Y coordinates corresponding to RNEW are calculated
C     as the proportion of RNEW / RAD.  FIXR is added to RAD in case
C     RAD is 0; FIXR = 1 iff RAD = 0, else FIXR = 0

         IF (RAD .EQ. 0.0) THEN
            FIXR = 1.0
         ELSE
            FIXR = 0.0
         END IF
C         FIXR  = SIGN(0.5, RAD) + SIGN(0.5, -RAD)

         PROP  = (RADA - DISTA(KLOA)) / (DISTA(KLOA+1) - DISTA(KLOA))
         RNEWA = RSA(KLOA) + PROP * (RSA(KLOA+1) - RSA(KLOA))
         PROP  = (RADB - DISTB(KLOB)) / (DISTB(KLOB+1) - DISTB(KLOB))
         RNEWB = RSB(KLOB) + PROP * (RSB(KLOB+1) - RSB(KLOB))

         H  = RSA(KLOA+1) - RSA(KLOA)
         A  = (RSA(KLOA+1)-RNEWA) / H
         B  = (RNEWA-RSA(KLOA)) / H

         ZT = A * ZSA(KLOA) + B * ZSA(KLOA+1) +
     *       ((A**3-A) * ZS2A(KLOA)+(B**3-B) * ZS2A(KLOA+1)) * (H**2)/6.

         H  = RSB(KLOB+1) - RSB(KLOB)
         A  = (RSB(KLOB+1)-RNEWB) / H
         B  = (RNEWB-RSB(KLOB)) / H

         ZB = A * ZSB(KLOB) + B * ZSB(KLOB+1) +
     *       ((A**3-A) * ZS2B(KLOB)+(B**3-B) * ZS2B(KLOB+1)) * (H**2)/6.

         IF (RDTHET) THEN
            RSAV  = RNEWA
            ZSAV  = ZT
            ZT    = ZSAV * COS(RSAV)
            RNEWA = ZSAV * SIN(RSAV)

            RSAV  = RNEWB
            ZSAV  = ZB
            ZB    = ZSAV * COS(RSAV)
            RNEWB = ZSAV * SIN(RSAV)
         END IF

         IF (SWEEP .EQ. SPHERI) THEN
C ... Spherical Sweep of Spline Surface
            XT = RNEWA * DX / (RAD + FIXR)
            YT = RNEWA * DY / (RAD + FIXR)
            XB = RNEWB * DX / (RAD + FIXR)
            YB = RNEWB * DY / (RAD + FIXR)

         ELSE IF (SWEEP .EQ. XSWEEP) THEN
C ... Sweep spline along the X axis
            XT = DX
            YT = RNEWA * DY / (RAD + FIXR)
            XB = DX
            YB = RNEWB * DY / (RAD + FIXR)

         ELSE IF (SWEEP .EQ. YSWEEP) THEN
C ... Sweep spline along the X axis
            XT = RNEWA * DX / (RAD + FIXR)
            YT = DY
            XB = RNEWB * DX / (RAD + FIXR)
            YB = DY

         END IF
         DELZ = ZT - ZB
         DELX = XT - XB
         DELY = YT - YB

         DO 90 NR = 1, NRNP(INP)
            ZN3(JNP0+NR) = ZT - DELZ * ZCORD(NR)
            YN3(JNP0+NR) = YT - DELY * ZCORD(NR)
            XN3(JNP0+NR) = XT - DELX * ZCORD(NR)
   90    CONTINUE
  100 CONTINUE

      RETURN
      END
