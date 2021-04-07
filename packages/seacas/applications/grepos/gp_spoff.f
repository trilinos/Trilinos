C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C     -*- Mode: fortran -*-
C=======================================================================
      SUBROUTINE SPOFF (XN, YN, ZN, NSPL, ZS, XS, XS2, YS, YS2,
     &     SLTOP, SLBOT, NUMNP, NDIM )
C=======================================================================

C     --*** SPLXYZ *** (GEN3D) Calculate 3D coordinates for Spline projection
C     --   Written by Amy Gilkey - revised 05/09/88
C     --   Modified by Greg Sjaardema - 02/06/89
C     --
C     --WRPXYZ calculates the coordinate array for the 3D warp translations.
C     --
C     --Parameters:
C     --   XN, YN, ZN - IN/OUT - the coordinates, modified
C     --   NSPL - Number of points on spline curve
C     --   ZS, XS, YS - IN - spline point coordinates
C     --   XS2, YS2 - spline slopes
C     --   SLTOP, SLBOT - slopes at top and bottom of spline
C     --   NUMNP - number of nodal points
C     --   NDIM - spatial dimension (2 or 3)
C     --

      REAL XN(*), YN(*), ZN(*)
      REAL ZS(NSPL), XS(NSPL), XS2(NSPL), YS(NSPL), YS2(NSPL)
      REAL SLTOP(2), SLBOT(2)

C$$$      CALL RDSPLN (.TRUE., NDB, SLTOP, SLBOT, NSPL, ZS, XS, YS)

      CALL SPLINE (ZS, XS, NSPL, SLTOP(1), SLBOT(1), XS2)
      CALL SPLINE (ZS, YS, NSPL, SLTOP(2), SLBOT(2), YS2)

      KLO = 1

      DO 80 INP = 1, NUMNP
         Z = ZN(INP)
         CALL HUNT (ZS, NSPL, Z, KLO)
         IF (KLO .EQ. NSPL) KLO = KLO - 1
         IF (KLO .EQ.    0) KLO = KLO + 1

         H  = ZS(KLO+1) - ZS(KLO)
         A  = (ZS(KLO+1)-Z)/H
         B  = (Z-ZS(KLO))/H

         XOFF = A * XS(KLO) + B * XS(KLO+1) +
     *        ((A**3-A) * XS2(KLO)+(B**3-B)*XS2(KLO+1)) * (H**2) / 6.

         YOFF = A * YS(KLO) + B * YS(KLO+1) +
     *        ((A**3-A) * YS2(KLO)+(B**3-B)*YS2(KLO+1)) * (H**2) / 6.

         XN(INP) = XN(INP) + XOFF
         YN(INP) = YN(INP) + YOFF
   80 CONTINUE

      RETURN
      END
