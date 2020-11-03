C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE TWIXYZ (XN, YN, XN3, YN3, ZN3, IXNP, NRNP, ZCORD,
     *   SINANG, COSANG)
C=======================================================================

C   --*** TWIXYZ *** (GEN3D) Calculate 3D coordinates for TWISTS
C   --   Modified by Greg Sjaardema - 02/06/89
C   --
C   --TWIXYZ calculates the coordinate array for the 3D database.
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
C   --   Uses DOTRAN, NNREPL, DIM3, NRTRAN, D3TRAN, ZGRAD,
C   --      CENTER, NUMCOL, NUMROW of /PARAMS/
C   --   Uses XOFFS, YOFFS, ZOFFS of /XYZOFF/
C   --   Uses ROT3D, ROTMAT of /XYZROT/

      INCLUDE 'g3_dbnums.blk'
      INCLUDE 'g3_dbnum3.blk'
      INCLUDE 'g3_params.blk'
      INCLUDE 'g3_twist.blk'

      REAL XN(NUMNP), YN(NUMNP),
     &   XN3(NUMNP3), YN3(NUMNP3), ZN3(NUMNP3)
      INTEGER IXNP(*), NRNP(*)
      REAL ZCORD(NNREPL)
      REAL SINANG(NNREPL), COSANG(NNREPL)

      PI = ATAN2(0.0, -1.0)
      D2R = PI / 180.0

C      --For translations, repeat the X coordinate and add a Z coordinate

      IF (ITWTYP .EQ. 1) THEN

C         --Calculate the equal or gradient step for current translation

         IBLK = 0
         NXTNR = 1

         ZEND = 0.0
   10    CONTINUE
         IBLK = IBLK + 1
         IF (NRTRAN(IBLK) .GT. 0) THEN
            ZBEG = ZEND
            ZEND = ZBEG + D3TRAN(IBLK)
            CALL INIGRD (-ZBEG, -ZEND, ZGRAD(IBLK),
     *         NRTRAN(IBLK), NRTRAN(IBLK)+1, ZCORD(NXTNR) )
            NXTNR = NXTNR + NRTRAN(IBLK)
            IF (IBLK .LT. MAXINT) GO TO 10
         END IF

         D3 = - (TWANGL / NEREPL) * D2R
         DO 20 NR = 1, NNREPL
            ANG = D3 * (NR-1)
            SINANG(NR) = SIN (ANG)
            COSANG(NR) = COS (ANG)
   20    CONTINUE

C      --Repeat the X coordinate and add the calculated Z coordinate

         DO 40 INP = 1, NUMNP
            JNP = IXNP(INP) - 1
            DO 30 NR = 1, NNREPL
               XT = COSANG(NR) * (XN(INP) - TWXCEN) -
     *            SINANG(NR) * (YN(INP) - TWYCEN) + TWXCEN
               YT = SINANG(NR) * (XN(INP) - TWXCEN) +
     *            COSANG(NR) * (YN(INP) - TWYCEN) + TWYCEN
               XN3(JNP+NR) = XT
               YN3(JNP+NR) = YT
               ZN3(JNP+NR) = ZCORD(NR)
   30       CONTINUE
   40    CONTINUE
      ELSE IF (ITWTYP .EQ. 2) THEN

C      --Get the coordinates for the non-center nodes (evenly spaced)

         IBLK = 0
         NXTNR = 1

         ZEND = 0.0
   50    CONTINUE
         IBLK = IBLK + 1
         IF (NRTRAN(IBLK) .GT. 0) THEN
            ZBEG = ZEND
            ZEND = ZBEG + D3TRAN(IBLK) * D2R
            CALL INIGRD (-ZBEG, -ZEND, ZGRAD(IBLK),
     *         NRTRAN(IBLK), NRTRAN(IBLK)+1, COSANG(NXTNR) )
            NXTNR = NXTNR + NRTRAN(IBLK)
            IF (IBLK .LT. MAXINT) GO TO 50
         END IF

         DO 60 NR = 1, NNREPL
            ANG = COSANG(NR)
            SINANG(NR) = SIN (ANG)
            COSANG(NR) = COS (ANG)
   60    CONTINUE

         D3 = - (TWANGL / NEREPL) * D2R
         DO 80 INP = 1, NUMNP
            IF (NRNP(INP) .EQ. NNREPL) THEN
               JNP = IXNP(INP)
               DO 70 NR = 1, NNREPL
                  TWANG = D3 * (NR-1)
                  SINTW = SIN(TWANG)
                  COSTW = COS(TWANG)
                  XT = COSTW * (XN(INP) - TWXCEN) -
     *               SINTW * (YN(INP) - TWYCEN) + TWXCEN
                  YT = SINTW * (XN(INP) - TWXCEN) +
     *               COSTW * (YN(INP) - TWYCEN) + TWYCEN
                  XN3(JNP) = (XT - CENTER) * COSANG(NR)
                  YN3(JNP) = YT
                  ZN3(JNP) = (XT - CENTER) * SINANG(NR)
                  JNP = JNP + 1
   70          CONTINUE
            END IF
   80    CONTINUE

C      --Add center of rotation

         IF (CENTER .NE. 0.0) THEN
            DO 90 JNP = 1, NUMNP3
               XN3(JNP) = XN3(JNP) + CENTER
   90       CONTINUE
         END IF
      END IF
      RETURN
      END
