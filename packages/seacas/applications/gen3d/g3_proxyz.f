C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PROXYZ (XN, YN, XN3, YN3, ZN3, IXNP, NRNP, ZCORD,
     *   SINANG, COSANG)
C=======================================================================

C   --*** PROXYZ *** (GEN3D) Calculate 3D coordinates for experimental
C   --   Modified by Greg Sjaardema - 02/06/89
C   --
C   --Parameters:
C   --   XN, YN - IN - the 2D coordinates, destroyed
C   --   XN3, YN3, ZN3 - OUT - the 3D coordinates
C   --   IXNP - IN - the new index for each node
C   --   NRNP - IN - the number of new nodes generated for each node
C   --   ZCORD - SCRATCH - size = NNREPL, holds z coordinate for transformations
C   --   SINANG, COSANG - SCRATCH - size = NNREPL, holds sin and cos of
C   --      angles for rotations
C   --
C   --Common Variables:
C   --   Uses NDIM, NUMNP of /DBNUMS/
C   --   Uses NDIM3, NUMNP3 of /DBNUM3/
C   --   Uses DOTRAN, NNREPL, DIM3, NRTRAN, D3TRAN, XXGRAD,
C   --      CENTER, NUMCOL, NUMROW of /PARAMS/
C   --   Uses XOFFS, YOFFS, ZOFFS of /XYZOFF/
C   --   Uses ROT3D, ROTMAT of /XYZROT/

      INCLUDE 'g3_dbnums.blk'
      INCLUDE 'g3_dbnum3.blk'
      INCLUDE 'g3_params.blk'
      INCLUDE 'g3_xxxxx.blk'

      REAL XN(NUMNP), YN(NUMNP),
     &   XN3(NUMNP3), YN3(NUMNP3), ZN3(NUMNP3)
      INTEGER IXNP(*), NRNP(*)
      REAL ZCORD(NNREPL)
      REAL SINANG(NNREPL), COSANG(NNREPL)

C      --Initialize the parametric interval distance

      ZTOT = 0.0
      DO 10 IBLK = 1, NBLK
         IF (NRTRAN(IBLK) .GT. 0) THEN
            ZTOT = ZTOT + D3TRAN(IBLK)
         ELSE
            CALL PRTERR ('PROGRAM',
     *         'Zero translations found')
            NBLK = IBLK
            GO TO 20
         END IF
   10 CONTINUE

      IF (ZTOT .EQ. 0.0) THEN
         CALL PRTERR ('CMDERR', 'Total translation distance is zero')
         STOP
      END IF

      NXTNR = 1
      ZEND = 0.0
   20 CONTINUE
      DO 30 IBLK = 1, NBLK
         ZBEG = ZEND
         ZEND = ZBEG + D3TRAN(IBLK)
         CALL INIGRD (ZBEG/ZTOT, ZEND/ZTOT, ZGRAD(IBLK),
     *      NRTRAN(IBLK), NRTRAN(IBLK)+1, ZCORD(NXTNR) )
         NXTNR = NXTNR + NRTRAN(IBLK)
   30 CONTINUE

C      --Project bottom surface onto a plane

      IF (ISXWRP .EQ. ISFLAT) THEN
         DO 50 INP = 1, NUMNP
            XB = XN(INP)
            YB = YN(INP)
            ZB = 0.0
            XT = (XN(INP) - XXSCL0) * XXSCAL + XXSCL0
            YT = (YN(INP) - XYSCL0) * XYSCAL + XYSCL0
            ZT =  - ZTOT + (XXA * XT + XXB * YT) / XXC
            XT =  XT + XXOFFS
            YT =  YT + XYOFFS

            JNP = IXNP(INP) - 1
            DO 40 NR = 1, NNREPL
               XN3(JNP+NR) = XB + (XT - XB) * ZCORD(NR)
               YN3(JNP+NR) = YB + (YT - YB) * ZCORD(NR)
               ZN3(JNP+NR) = ZB + (ZT - ZB) * ZCORD(NR)
   40       CONTINUE
   50    CONTINUE

C ... Warp of translated surface

      ELSE IF (ISXWRP .EQ. ISSPHE) THEN
         IF (CONVEX) THEN
            ZCEN = ZTOT - XWARP
            RMULT = 1.0
         ELSE
            ZCEN = ZTOT + XWARP
            RMULT = -1.0
         END IF
         DO 70 INP = 1, NUMNP
            XB = XN(INP)
            YB = YN(INP)
            ZB = 0.0
            XT = (XN(INP) - XXSCL0) * XXSCAL + XXSCL0
            YT = (YN(INP) - XYSCL0) * XYSCAL + XYSCL0

C ... Note: ZT must be negative

            ZT =  - (ZCEN + RMULT * SQRT(XWARP**2 - XB**2 - YB**2))
            XT =  XT + XXOFFS
            YT =  YT + XYOFFS

            JNP = IXNP(INP) - 1
            DO 60 NR = 1, NNREPL
               XN3(JNP+NR) = XB + (XT - XB) * ZCORD(NR)
               YN3(JNP+NR) = YB + (YT - YB) * ZCORD(NR)
               ZN3(JNP+NR) = ZB + (ZT - ZB) * ZCORD(NR)
   60       CONTINUE
   70    CONTINUE

C ... Toroidal Surface
      ELSE IF (ISXWRP .EQ. ISTORO) THEN
         IF (CONVEX) THEN
            ZCEN = ZTOT - XWARP - YWARP
            RMULT = 1.0
         ELSE
            ZCEN = ZTOT + XWARP + YWARP
            RMULT = -1.0
         END IF
         DO 90 INP = 1, NUMNP
            XB = XN(INP)
            YB = YN(INP)
            ZB = 0.0
            XT = (XN(INP) - XXSCL0) * XXSCAL + XXSCL0
            YT = (YN(INP) - XYSCL0) * XYSCAL + XYSCL0

C ... Note: ZT must be negative

            ZA =  ZCEN + RMULT * SQRT(XWARP**2 - YB**2)
            ZT =  - (ZA + RMULT * SQRT(YWARP**2 - XB**2))
            XT =  XT + XXOFFS
            YT =  YT + XYOFFS

            JNP = IXNP(INP) - 1
            DO 80 NR = 1, NNREPL
               XN3(JNP+NR) = XB + (XT - XB) * ZCORD(NR)
               YN3(JNP+NR) = YB + (YT - YB) * ZCORD(NR)
               ZN3(JNP+NR) = ZB + (ZT - ZB) * ZCORD(NR)
   80       CONTINUE
   90    CONTINUE

C ... Cylinder with Axis on X Axis
      ELSE IF (ISXWRP .EQ. ISXCYL) THEN
         IF (CONVEX) THEN
            ZCEN = ZTOT - XWARP
            RMULT = 1.0
         ELSE
            ZCEN = ZTOT + XWARP
            RMULT = -1.0
         END IF
         DO 110 INP = 1, NUMNP
            XB = XN(INP)
            YB = YN(INP)
            ZB = 0.0
            XT = (XN(INP) - XXSCL0) * XXSCAL + XXSCL0
            YT = (YN(INP) - XYSCL0) * XYSCAL + XYSCL0

C ... Note: ZT must be negative

            ZT =  - (ZCEN + RMULT * SQRT(XWARP**2 - YB**2))
            XT =  XT + XXOFFS
            YT =  YT + XYOFFS

            JNP = IXNP(INP) - 1
            DO 100 NR = 1, NNREPL
               XN3(JNP+NR) = XB + (XT - XB) * ZCORD(NR)
               YN3(JNP+NR) = YB + (YT - YB) * ZCORD(NR)
               ZN3(JNP+NR) = ZB + (ZT - ZB) * ZCORD(NR)
 100         CONTINUE
 110       CONTINUE

C ... Cylinder with Axis on Y Axis
      ELSE IF (ISXWRP .EQ. ISYCYL) THEN
         IF (CONVEX) THEN
            ZCEN = ZTOT - YWARP
            RMULT = 1.0
         ELSE
            ZCEN = ZTOT + YWARP
            RMULT = -1.0
         END IF
         DO 130 INP = 1, NUMNP
            XB = XN(INP)
            YB = YN(INP)
            ZB = 0.0
            XT = (XN(INP) - XXSCL0) * XXSCAL + XXSCL0
            YT = (YN(INP) - XYSCL0) * XYSCAL + XYSCL0

C ... Note: ZT must be negative

            ZT =  - (ZCEN + RMULT * SQRT(YWARP**2 - XB**2))
            XT =  XT + XXOFFS
            YT =  YT + XYOFFS

            JNP = IXNP(INP) - 1
            DO 120 NR = 1, NNREPL
               XN3(JNP+NR) = XB + (XT - XB) * ZCORD(NR)
               YN3(JNP+NR) = YB + (YT - YB) * ZCORD(NR)
               ZN3(JNP+NR) = ZB + (ZT - ZB) * ZCORD(NR)
 120         CONTINUE
 130       CONTINUE
      END IF

      RETURN
      END
