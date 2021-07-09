C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE NEWXYZ (XN, YN, XN3, YN3, ZN3, IXNP, NRNP, NPCEN,
     &   ZCORD, SINANG, COSANG, A)
C=======================================================================

C   --*** NEWXYZ *** (GEN3D) Calculate 3D coordinates
C   --   Written by Amy Gilkey - revised 05/09/88
C   --   Modified by Greg Sjaardema - 02/06/89
C   --       Added Warp Function
C   --       Added Gradient to Rotations (not for center blocks)
C   --       Split transformations into separate subroutines
C   --
C   --NEWXYZ calculates the coordinate array for the 3D database.
C   --
C   --Parameters:
C   --   XN, YN - IN - the 2D coordinates, destroyed
C   --   XN3, YN3, ZN3 - OUT - the 3D coordinates
C   --   IXNP - IN - the new index for each node
C   --   NRNP - IN - the number of new nodes generated for each node
C   --   NPCEN - IN - the node numbers of the center nodes by column and row
C   --   ZCORD - SCRATCH - size = NNREPL, holds z coordinate for transformations
C   --   SINANG, COSANG - SCRATCH - size = NNREPL, holds sin and cos of
C   --      angles for rotations
C   --
C   --Common Variables:
C   --   Uses NDIM, NUMNP of /DBNUMS/
C   --   Uses NDIM3, NUMNP3 of /DBNUM3/
C   --   Uses ITRANT, NNREPL, DIM3, NRTRAN, D3TRAN, ZGRAD,
C   --      CENTER, NUMCOL, NUMROW of /PARAMS/
C   --   Uses XOFFS, YOFFS, ZOFFS of /XYZOFF/
C   --   Uses ROT3D, ROTMAT of /XYZROT/

      INCLUDE 'g3_dbnums.blk'
      INCLUDE 'g3_dbnum3.blk'
      INCLUDE 'g3_params.blk'
      INCLUDE 'g3_xyzoff.blk'
      INCLUDE 'g3_xyzrot.blk'
      INCLUDE 'g3_xyzmir.blk'
      INCLUDE 'g3_xyzero.blk'
      INCLUDE 'g3_xyzscl.blk'
      INCLUDE 'g3_twist.blk'
      INCLUDE 'g3_splxyz.blk'

      REAL XN(NUMNP), YN(NUMNP),
     &   XN3(NUMNP3), YN3(NUMNP3), ZN3(NUMNP3)
      INTEGER IXNP(*), NRNP(*)
      INTEGER NPCEN(NUMCDM,*)
      REAL ZCORD(NNREPL)
      REAL SINANG(NNREPL), COSANG(NNREPL)
      REAL A(*)

      IF (ITRANT .EQ. 1) THEN
         CALL TRNXYZ (XN, YN, XN3, YN3, ZN3, IXNP, NRNP, ZCORD)
       ELSE IF (ITRANT .EQ. 2) THEN
         if (rotax .eq. 0) then
           CALL ARCXYZ (XN, YN, XN3, YN3, ZN3, IXNP, NRNP, NPCEN,
     &       SINANG, COSANG)
         else
           CALL ARCYXZ (XN, YN, XN3, YN3, ZN3, IXNP, NRNP, NPCEN,
     &       SINANG, COSANG)
         end if
      ELSE IF (ITRANT .EQ. 4) THEN
         CALL WRPXYZ (XN, YN, XN3, YN3, ZN3, IXNP, NRNP, ZCORD  )
      ELSE IF (ITRANT .EQ. 8) THEN
         CALL TWIXYZ (XN, YN, XN3, YN3, ZN3, IXNP, NRNP, ZCORD,
     *      SINANG, COSANG)
      ELSE IF (ITRANT .EQ. 16) THEN
         CALL PROXYZ (XN, YN, XN3, YN3, ZN3, IXNP, NRNP, ZCORD,
     *      SINANG, COSANG)
      ELSE IF (ITRANT .EQ. 32) THEN
         CALL EXPARC (XN, YN, XN3, YN3, ZN3, IXNP, NRNP, NPCEN,
     &      SINANG, COSANG)
      ELSE IF (ITRANT .EQ. 64) THEN
         CALL SPLXYZ (XN, YN, XN3, YN3, ZN3, IXNP, NRNP, ZCORD,
     $      NSPL(1), NSPL(2),
     &      A(KRSPLA), A(KZSPLA), A(KSPL2A), A(KDISTA), A(KSCRA),
     &      A(KRSPLB), A(KZSPLB), A(KSPL2B), A(KDISTB), A(KSCRB),
     $      SLLFT, SLRGT, RDTHET, SWEEP, NOSCAL)
      ELSE IF (ITRANT .EQ. 128) THEN
         CALL SPTXYZ (XN, YN, XN3, YN3, ZN3, IXNP, NRNP, ZCORD, NSPLT,
     $        A(KZSPL), A(KXSPL), A(KXSPL2), A(KYSPL), A(KYSPL2),
     $        A(KSCR))
      END IF

C   --Revolve 3D mesh, if needed

      IF (ROT3D) THEN
         DO 30 JNP = 1, NUMNP3
            X = XN3(JNP) - ROTCEN(1)
            Y = YN3(JNP) - ROTCEN(2)
            Z = ZN3(JNP) - ROTCEN(3)
            XN3(JNP) = X*ROTMAT(1,1) + Y*ROTMAT(2,1) + Z*ROTMAT(3,1)
     &         + ROTCEN(1)
            YN3(JNP) = X*ROTMAT(1,2) + Y*ROTMAT(2,2) + Z*ROTMAT(3,2)
     &         + ROTCEN(2)
            ZN3(JNP) = X*ROTMAT(1,3) + Y*ROTMAT(2,3) + Z*ROTMAT(3,3)
     &         + ROTCEN(3)
   30    CONTINUE
      END IF

C   --Add offset, if any, to coordinates

      IF (XOFFS .NE. 0.0) THEN
         DO 40 JNP = 1, NUMNP3
            XN3(JNP) = XN3(JNP) + XOFFS
   40    CONTINUE
      END IF
      IF (YOFFS .NE. 0.0) THEN
         DO 50 JNP = 1, NUMNP3
            YN3(JNP) = YN3(JNP) + YOFFS
   50    CONTINUE
      END IF
      IF (ZOFFS .NE. 0.0) THEN
         DO 60 JNP = 1, NUMNP3
            ZN3(JNP) = ZN3(JNP) + ZOFFS
   60    CONTINUE
      END IF

C   --Mirror coordinates if any specified

      IF (XMIRR .LT. 0.) THEN
         DO 70 JNP = 1, NUMNP3
            XN3(JNP) = -1.0 * XN3(JNP)
   70    CONTINUE
      END IF
      IF (YMIRR .LT. 0.) THEN
         DO 80 JNP = 1, NUMNP3
            YN3(JNP) = -1.0 * YN3(JNP)
   80    CONTINUE
      END IF
      IF (ZMIRR .LT. 0.) THEN
         DO 90 JNP = 1, NUMNP3
            ZN3(JNP) = -1.0 * ZN3(JNP)
   90    CONTINUE
      END IF

C --- Zero coordinates if *ZERO is not equal to zero

      IF (XZERO .NE. 0.) THEN
         DO 100 JNP = 1, NUMNP3
           if (ABS(XN3(JNP)) .LT. XZERO) XN3(JNP) = 0.0
  100    CONTINUE
      END IF
      IF (YZERO .NE. 0.) THEN
         DO 110 JNP = 1, NUMNP3
           if (ABS(YN3(JNP)) .LT. YZERO) YN3(JNP) = 0.0
  110    CONTINUE
      END IF
      IF (ZZERO .NE. 0.) THEN
         DO 120 JNP = 1, NUMNP3
           if (ABS(ZN3(JNP)) .LT. ZZERO) ZN3(JNP) = 0.0
  120    CONTINUE
      END IF

C --- Scale the coordinates if any Scaled

      IF (XSCAL .NE. 1.) THEN
         DO 130 JNP = 1, NUMNP3
            XN3(JNP) = XSCAL * XN3(JNP)
  130    CONTINUE
      END IF
      IF (YSCAL .NE. 1.) THEN
         DO 140 JNP = 1, NUMNP3
            YN3(JNP) = YSCAL * YN3(JNP)
  140    CONTINUE
      END IF
      IF (ZSCAL .NE. 1. .AND. NDIM .EQ. 3) THEN
         DO 150 JNP = 1, NUMNP3
            ZN3(JNP) = ZSCAL * ZN3(JNP)
  150    CONTINUE
      END IF

      CALL MINMAX (NUMNP3, XN3, XMIN, XMAX)
      CALL MINMAX (NUMNP3, YN3, YMIN, YMAX)
      CALL MINMAX (NUMNP3, ZN3, ZMIN, ZMAX)

      WRITE (*, 155) 'Output Mesh Limits:'
      WRITE (*, 160) 'X', XMIN, 'X', XMAX, XMAX-XMIN
      WRITE (*, 160) 'Y', YMIN, 'Y', YMAX, YMAX-YMIN
      WRITE (*, 160) 'Z', ZMIN, 'Z', ZMAX, ZMAX-ZMIN
 155  FORMAT(/' ',A)
 160  FORMAT( ' Minimum ',A1,' = ',1PE12.5,', Maximum ',A1,' = ',
     &     1PE12.5,', Range = ',1PE12.5)

      RETURN
      END

