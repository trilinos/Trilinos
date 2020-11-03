C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE WRPXYZ (XN, YN, XN3, YN3, ZN3, ATRIB)
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
C   --
C   --Common Variables:
C   --   Uses NDIM, NUMNP of /DBNUMS/
C   --   Uses NDIM3, NUMNP3 of /DBNUM3/
C   --   Uses DOTRAN, NNREPL, DIM3, NRTRAN, D3TRAN, ZGRAD,
C   --      CENTER, NUMCOL, NUMROW of /PARAMS/

      INCLUDE 'gs_dbnums.blk'
      INCLUDE 'gs_dbnum3.blk'
      INCLUDE 'gs_params.blk'

      REAL XN(NUMNP), YN(NUMNP),
     &   XN3(NUMNP3), YN3(NUMNP3), ZN3(NUMNP3)
      REAL ATRIB(NUMEL)

      IF (IWARP .EQ. 1) THEN

C ... Warp type 1: Point Centered

         DO 100 INP = 1, NUMNP
            XN3(INP) = XN(INP)
            YN3(INP) = YN(INP)
            ZN3(INP) = DWARP - SQRT(DWARP**2 - XN(INP)**2 - YN(INP)**2)
  100    CONTINUE

         CONTINUE
      ELSE IF (IWARP .EQ. -1) THEN

C ... Warp type -1: X Axis Centered

         DO 110 INP = 1, NUMNP
            THET = YN(INP) / DWARP
            XN3(INP) = XN(INP)
            YN3(INP) = SIN(THET) * DWARP
            ZN3(INP) = DWARP - COS(THET) * DWARP
  110    CONTINUE

      ELSE IF (IWARP .EQ. -2) THEN

C ... Warp type -2: Y Axis Centered

         DO 120 INP = 1, NUMNP
            THET = XN(INP) / DWARP
            XN3(INP) = SIN(THET) * DWARP
            YN3(INP) = YN(INP)
            ZN3(INP) = DWARP - COS(THET) * DWARP
  120    CONTINUE

      ELSE IF (IWARP .EQ. -3) THEN

C ... Warp type -3: X Axis Centered, Project straight up

         DO 130 INP = 1, NUMNP
            XN3(INP) = XN(INP)
            YN3(INP) = YN(INP)
            ZN3(INP) = DWARP - sqrt( dwarp**2 - yn(inp)**2 )
  130    CONTINUE

      ELSE IF (IWARP .EQ. -4) THEN

C ... Warp type -4: Y Axis Centered, Project straight up

         DO 140 INP = 1, NUMNP
            XN3(INP) = XN(INP)
            YN3(INP) = YN(INP)
            ZN3(INP) = DWARP - sqrt( dwarp**2 - xn(inp)**2 )
  140    CONTINUE

      ELSE IF (IWARP .EQ. 2) THEN

C ... Warp type 2: Point-Centered Ellipse

         DO 150 INP = 1, NUMNP
            DX = XN(INP)
            DY = YN(INP)
            ZT = DWARP * (1.0-SQRT(1.0-DX**2/XRAD**2-DY**2/YRAD**2))
            ZN3(INP) = ZT
            YN3(INP) = YN(INP)
            XN3(INP) = XN(INP)
  150    CONTINUE
      END IF

C ... Now do the attributes for all of the elements

      DO 160 IEL = 1, NUMEL
         ATRIB(IEL) = DIM3
  160 CONTINUE

      RETURN
      END
