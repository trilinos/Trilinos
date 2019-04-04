C Copyright(C) 2011-2017 National Technology & Engineering Solutions of
C Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C
C * Redistributions of source code must retain the above copyright
C    notice, this list of conditions and the following disclaimer.
C
C * Redistributions in binary form must reproduce the above
C   copyright notice, this list of conditions and the following
C   disclaimer in the documentation and/or other materials provided
C   with the distribution.
C
C * Neither the name of NTESS nor the names of its
C   contributors may be used to endorse or promote products derived
C   from this software without specific prior written permission.
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
