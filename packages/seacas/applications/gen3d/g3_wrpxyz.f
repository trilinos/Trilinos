C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE WRPXYZ (XN, YN, XN3, YN3, ZN3, IXNP, NRNP, ZCORD  )
C=======================================================================

C   --*** WRPXYZ *** (GEN3D) Calculate 3D coordinates
C   --   Written by Amy Gilkey - revised 05/09/88
C   --   Modified by Greg Sjaardema - 02/06/89
C   --       Added Warp Function
C   --       Added Gradient to Rotations (not for center blocks)
C   --
C   --WRPXYZ calculates the coordinate array for the 3D warp translations.
C   --
C   --Parameters:
C   --   XN, YN - IN - the 2D coordinates, destroyed
C   --   XN3, YN3, ZN3 - OUT - the 3D coordinates
C   --   IXNP - IN - the new index for each node
C   --   NRNP - IN - the number of new nodes generated for each node
C   --   ZCORD - SCRATCH - size = NNREPL, holds z coordinate for transformations
C   --
C   --Common Variables:
C   --   Uses NDIM, NUMNP of /DBNUMS/
C   --   Uses NDIM3, NUMNP3 of /DBNUM3/
C   --   Uses DOTRAN, NNREPL, DIM3, NRTRAN, D3TRAN, ZGRAD,
C   --      CENTER, NUMCOL, NUMROW of /PARAMS/

      INCLUDE 'g3_dbnums.blk'
      INCLUDE 'g3_dbnum3.blk'
      INCLUDE 'g3_params.blk'

      REAL XN(NUMNP), YN(NUMNP),
     &   XN3(NUMNP3), YN3(NUMNP3), ZN3(NUMNP3)
      INTEGER IXNP(*), NRNP(*)
      REAL ZCORD(NNREPL)

C ... Doing a Warp
C ... CALCULATE THE THICKNESS INCREMENT FOR EACH TRANSLATION

      IF (VEDGE) THEN
         CALL INIGRD (0.0, 1.0, ZGRAD(1), NRTRAN(1), NNREPL, ZCORD)
      ELSE
         CALL INIGRD (DWARP, DWARP+D3TRAN(1), ZGRAD(1), NRTRAN(1),
     *      NNREPL, ZCORD)
      END IF

C      NOGRAD = (ABS (ZGRAD(1) - 1.0) .LE. 0.001)
C      IF (NOGRAD) THEN
C         D3 = 1.0 / NRTRAN(1)
C      ELSE
C         D3 = (1.0 - ZGRAD(1)) / (1.0 - ZGRAD(1)**NRTRAN(1))
C      END IF

C      IF (.NOT. VEDGE) THEN
C         D3 = D3 * D3TRAN(1)
C      END IF

C      IF (NOGRAD) THEN
C         IF (VEDGE) THEN
C            DO 10 NR = 1, NNREPL
C               ZCORD(NR) = D3 * (NR-1)
C   10       CONTINUE
C         ELSE
C            DO 20 NR = 1, NNREPL
C               ZCORD(NR) = DWARP + D3 * (NR-1)
C   20       CONTINUE
C         END IF
C      ELSE
C         IF (VEDGE) THEN
C            ZCORD(1) = 0.0
C         ELSE
C            ZCORD(1) = DWARP
C         END IF
C         DO 30 NR = 2, NNREPL
C            ZCORD(NR) = ZCORD(NR-1) + D3 * ZGRAD(1)**(NR-2)
C   30    CONTINUE
C      END IF

      IF (IWARP .EQ. 1) THEN

C ... Warp type 1: Point Centered

         DO 60 INP = 1, NUMNP
            JNP0 = IXNP(INP) - 1
            DX = XN(INP)
            DY = YN(INP)
            ZT = DWARP - SQRT(DWARP**2 - DX**2 - DY**2)
            IF (VEDGE) THEN
               DWARPB = DWARP + D3TRAN(1)
               ZB = DWARP - SQRT(DWARPB**2 - DX**2 - DY**2)
               DO 40 NR = 1, NRNP(INP)
                  ZN3(JNP0+NR) = ZT + (ZB - ZT) * ZCORD(NR)
                  YN3(JNP0+NR) = YN(INP)
                  XN3(JNP0+NR) = XN(INP)
   40          CONTINUE
            ELSE
               C1 = SQRT (XN(INP)**2 + YN(INP)**2 + (ZT-DWARP)**2)
               DO 50 NR = 1, NRNP(INP)
                  TNR = ZCORD(NR) / C1
                  ZN3(JNP0+NR) = DWARP + (ZT - DWARP) * TNR
                  YN3(JNP0+NR) = YN(INP) * TNR
                  XN3(JNP0+NR) = XN(INP) * TNR
   50          CONTINUE
            END IF
   60    CONTINUE

         CONTINUE
      ELSE IF (IWARP .EQ. -1) THEN

C ... Warp type -1: X Axis Centered

         DO 90 INP = 1, NUMNP
            JNP0 = IXNP(INP) - 1
            THET = YN(INP) / DWARP
            C1 = SIN(THET)
            C2 = COS(THET)
            IF (VEDGE) THEN
               XT = XN(INP)
               YT = C1 * DWARP
               ZT = DWARP - C2 * DWARP

               DWARPB = DWARP + D3TRAN(1)
               ZB = DWARP - SQRT(ABS(DWARPB**2 - YT**2))
               XB = XT
               YB = YT

               DO 70 NR = 1, NRNP(INP)
                  XN3(JNP0+NR) = XT + (XB - XT) * ZCORD(NR)
                  YN3(JNP0+NR) = YT + (YB - YT) * ZCORD(NR)
                  ZN3(JNP0+NR) = ZT + (ZB - ZT) * ZCORD(NR)
   70          CONTINUE
            ELSE
               DO 80 NR = 1, NRNP(INP)
                  XN3(JNP0+NR) = XN(INP)
                  YN3(JNP0+NR) = C1 * ZCORD(NR)
                  ZN3(JNP0+NR) = DWARP - C2 * ZCORD(NR)
   80          CONTINUE
            END IF
   90    CONTINUE

      ELSE IF (IWARP .EQ. -2) THEN

C ... Warp type -2: Y Axis Centered

         DO 120 INP = 1, NUMNP
            JNP0 = IXNP(INP) - 1
            THET = XN(INP) / DWARP
            C1 = SIN(THET)
            C2 = COS(THET)
            IF (VEDGE) THEN
               XT = C1 * DWARP
               YT = YN(INP)
               ZT = DWARP - C2 * DWARP

               DWARPB = DWARP + D3TRAN(1)
               ZB = DWARP - SQRT(ABS(DWARPB**2 - XT**2))
               XB = XT
               YB = YT

               DO 100 NR = 1, NRNP(INP)
                  XN3(JNP0+NR) = XT + (XB - XT) * ZCORD(NR)
                  YN3(JNP0+NR) = YT + (YB - YT) * ZCORD(NR)
                  ZN3(JNP0+NR) = ZT + (ZB - ZT) * ZCORD(NR)
  100          CONTINUE
            ELSE
               DO 110 NR = 1, NRNP(INP)
                  XN3(JNP0+NR) = C1 * ZCORD(NR)
                  YN3(JNP0+NR) = YN(INP)
                  ZN3(JNP0+NR) = DWARP - C2 * ZCORD(NR)
  110          CONTINUE
            END IF
  120    CONTINUE

      ELSE IF (IWARP .EQ. 2) THEN

C ... Warp type 1: Point-Centered Ellipse

         DO 360 INP = 1, NUMNP
            JNP0 = IXNP(INP) - 1
            DX = XN(INP)
            DY = YN(INP)
            ZT = DWARP - DWARP / HRAD * SQRT(HRAD**2 - DX**2 - DY**2)
            IF (VEDGE) THEN
               DWARPB = DWARP + D3TRAN(1)
               ZB = DWARP - DWARPB/HRAD * SQRT(HRAD**2 - DX**2 - DY**2)
               DO 340 NR = 1, NRNP(INP)
                  ZN3(JNP0+NR) = ZT + (ZB - ZT) * ZCORD(NR)
                  YN3(JNP0+NR) = YN(INP)
                  XN3(JNP0+NR) = XN(INP)
  340          CONTINUE
            ELSE
               C1 = SQRT (XN(INP)**2 + YN(INP)**2 + (ZT-DWARP)**2)
               DO 350 NR = 1, NRNP(INP)
                  TNR = ZCORD(NR) / C1
                  ZN3(JNP0+NR) = DWARP + (ZT - DWARP) * TNR
                  YN3(JNP0+NR) = YN(INP) * TNR
                  XN3(JNP0+NR) = XN(INP) * TNR
  350          CONTINUE
            END IF
  360    CONTINUE

         CONTINUE
      END IF
      RETURN
      END
