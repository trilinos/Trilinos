C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SPTXYZ (XN, YN, XN3, YN3, ZN3, IXNP, NRNP, ZCORD,
     &     NSPL, ZS, XS, XS2, YS, YS2, SCR )
C=======================================================================

C     --*** SPLTYZ *** (GEN3D) Calculate 3D coordinates for Spline translation
C     --   Written by Greg Sjaardema - 02/06/89
C     --
C     --Parameters:
C     --   XN, YN - IN - the 2D coordinates, destroyed
C     --   XN3, YN3, ZN3 - OUT - the 3D coordinates
C     --   IXNP - IN - the new index for each node
C     --   NRNP - IN - the number of new nodes generated for each node
C     --   ZCORD - SCRATCH - size = NNREPL, holds z coordinate for transformations
C     --   NSPL - IN - number of spline points
C     --   ZS, XS, YS - IN - spline data points
C     --   XS2, YS2 - SCRATCH - temporary storage for spline data
C     --   SCR - SCRATCH - temporary storage for spoine data
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

      REAL XN(NUMNP), YN(NUMNP),
     &     XN3(NUMNP3), YN3(NUMNP3), ZN3(NUMNP3)
      INTEGER IXNP(*), NRNP(*)
      REAL ZCORD(NNREPL)
      REAL ZS(NSPL), XS(NSPL), XS2(NSPL), YS(NSPL), YS2(NSPL)
      REAL SCR(NSPL), SLLFT(2), SLRGT(2)

C     ... CALCULATE THE THICKNESS INCREMENT FOR EACH TRANSLATION

      IBLK = 0
      ZTOT = 0.0
   10 CONTINUE
      IBLK = IBLK + 1
      IF (NRTRAN(IBLK) .GT. 0) THEN
         ZTOT = ZTOT + D3TRAN(IBLK)
         IF (IBLK .LT. MAXINT) GO TO 10
      END IF

      NXTNR = 1
      IBLK = 0
      ZEND = 0.0

   20 CONTINUE
      IBLK = IBLK + 1
      IF (NRTRAN(IBLK) .GT. 0) THEN
         ZBEG = ZEND
         ZEND = ZBEG + D3TRAN(IBLK)
         CALL INIGRD (ZBEG/ZTOT, ZEND/ZTOT, ZGRAD(IBLK),
     *        NRTRAN(IBLK), NRTRAN(IBLK)+1, ZCORD(NXTNR) )
         NXTNR = NXTNR + NRTRAN(IBLK)
         IF (IBLK .LT. MAXINT) GO TO 20
      END IF

      CALL SPLINE (ZS, XS, NSPL, SLLFT(1), SLRGT(1), XS2, SCR)
      CALL SPLINE (ZS, YS, NSPL, SLLFT(2), SLRGT(2), YS2, SCR)

      ZMAX = ZS(NSPL)

      KLO = 1

      DO 80 NR = 1, NNREPL
         Z = ZCORD(NR) * ZMAX
         CALL HUNT (ZS, NSPL, Z, KLO)
         KLO = MIN(NSPL-1, MAX(1, KLO))

         H  = ZS(KLO+1) - ZS(KLO)
         A  = (ZS(KLO+1)-Z)/H
         B  = (Z-ZS(KLO))/H

         XOFF = A * XS(KLO) + B * XS(KLO+1) +
     *        ((A**3-A) * XS2(KLO)+(B**3-B)*XS2(KLO+1)) * (H**2) / 6.

         YOFF = A * YS(KLO) + B * YS(KLO+1) +
     *        ((A**3-A) * YS2(KLO)+(B**3-B)*YS2(KLO+1)) * (H**2) / 6.

         DO 70 INP = 1, NUMNP
            JNP = IXNP(INP) - 1
            XN3(JNP+NR) = XN(INP) + XOFF
            YN3(JNP+NR) = YN(INP) + YOFF
            ZN3(JNP+NR) = Z
   70    CONTINUE
   80 CONTINUE

      RETURN
      END

