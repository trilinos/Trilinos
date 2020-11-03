C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C -*- Mode: fortran -*-
C=======================================================================
      SUBROUTINE NEWXYZ (XN, YN, XN3, YN3, ZN3, ATRIB, A)
C=======================================================================

C   --*** NEWXYZ *** (GENSHELL) Calculate 3D coordinates
C   --
C   --NEWXYZ calculates the coordinate array for the 3D database.
C   --
C   --Parameters:
C   --   XN, YN - IN - the 2D coordinates, destroyed
C   --   XN3, YN3, ZN3 - OUT - the 3D coordinates
C   --
C   --Common Variables:
C   --   Uses NDIM, NUMNP of /DBNUMS/
C   --   Uses NDIM3, NUMNP3 of /DBNUM3/
C   --   Uses ITRANT, NNREPL, DIM3, NRTRAN, D3TRAN, ZGRAD,
C   --      CENTER, NUMCOL, NUMROW of /PARAMS/
C   --   Uses XOFFS, YOFFS, ZOFFS of /XYZOFF/
C   --   Uses ROT3D, ROTMAT of /XYZROT/

      INCLUDE 'gs_dbnums.blk'
      INCLUDE 'gs_dbnum3.blk'
      INCLUDE 'gs_params.blk'
      INCLUDE 'gs_xyzoff.blk'
      INCLUDE 'gs_xyzrot.blk'
      INCLUDE 'gs_xyzmir.blk'
      INCLUDE 'gs_xyzero.blk'
      INCLUDE 'gs_xyzscl.blk'
      INCLUDE 'gs_splxyz.blk'

      REAL XN(NUMNP), YN(NUMNP),
     &   XN3(NUMNP3), YN3(NUMNP3), ZN3(NUMNP3)
      REAL ATRIB(*)
      REAL A(*)

      IF (ITRANT .EQ. 1) THEN
         CALL TRNXYZ (XN, YN, XN3, YN3, ZN3, ATRIB)
C$$$      ELSE IF (ITRANT .EQ. 2) THEN
C$$$         CALL ARCXYZ (XN, YN, XN3, YN3, ZN3, ATRIB)
      ELSE IF (ITRANT .EQ. 4) THEN
         CALL WRPXYZ (XN, YN, XN3, YN3, ZN3, ATRIB)
C$$$      ELSE IF (ITRANT .EQ. 8) THEN
C$$$         CALL TWIXYZ (XN, YN, XN3, YN3, ZN3, ATRIB)
c$$$      ELSE IF (ITRANT .EQ. 16) THEN
c$$$         CALL PROXYZ (XN, YN, XN3, YN3, ZN3, ATRIB)
c$$$      ELSE IF (ITRANT .EQ. 32) THEN
c$$$         CALL EXPARC (XN, YN, XN3, YN3, ZN3, ATRIB)
      ELSE IF (ITRANT .EQ. 64) THEN
         CALL SPLXYZ (XN, YN, XN3, YN3, ZN3, ATRIB,
     $      NSPL(1), A(KRSPLA), A(KZSPLA), A(KSPL2A), A(KDISTA),
     $        A(KSCRA), SLLFT, SLRGT, RDTHET, SWEEP, NOSCAL)
c$$$      ELSE IF (ITRANT .EQ. 128) THEN
c$$$         CALL SPTXYZ (XN, YN, XN3, YN3, ZN3, ATRIB, NSPLT,
c$$$     $        A(KZSPL), A(KXSPL), A(KXSPL2), A(KYSPL), A(KYSPL2),
c$$$     $        A(KSCR))
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

C     --Randomize coordinates if any specified

      IDUM = 1
      IF (XRAND .NE. 0.) THEN
         DO 61 JNP = 1, NUMNP
            XN3(JNP) = (2.0*RAN1(IDUM)-1.0) * XRAND + XN3(JNP)
 61      CONTINUE
      END IF
      IF (YRAND .NE. 0.) THEN
         DO 71 JNP = 1, NUMNP
            YN3(JNP) = (2.0*RAN1(IDUM)-1.0) * YRAND + YN3(JNP)
 71      CONTINUE
      END IF
      IF (ZRAND .NE. 0.0) THEN
         DO 81 JNP = 1, NUMNP
            ZN3(JNP) = (2.0*RAN1(IDUM)-1.0) * ZRAND + ZN3(JNP)
 81      CONTINUE
      END IF

C --- Zero coordinates if *ZERO is not equal to zero

      IF (XZERO .NE. 0.) THEN
         DO 100 JNP = 1, NUMNP3
            XN3(JNP) = (0.5+SIGN( 0.5, ABS(XN3(JNP))-XZERO)) * XN3(JNP)
  100    CONTINUE
      END IF
      IF (YZERO .NE. 0.) THEN
         DO 110 JNP = 1, NUMNP3
            YN3(JNP) = (0.5+SIGN( 0.5, ABS(YN3(JNP))-YZERO)) * YN3(JNP)
  110    CONTINUE
      END IF
      IF (ZZERO .NE. 0.) THEN
         DO 120 JNP = 1, NUMNP3
            ZN3(JNP) = (0.5+SIGN( 0.5, ABS(ZN3(JNP))-ZZERO)) * ZN3(JNP)
  120    CONTINUE
      END IF

C --- Scale the coordinates if any Scaled

      IF (XSCAL .NE. 1.) THEN
         DO 130 JNP = 1, NUMNP
            XN3(JNP) = XSCAL * XN3(JNP)
  130    CONTINUE
      END IF
      IF (YSCAL .NE. 1.) THEN
         DO 140 JNP = 1, NUMNP
            YN3(JNP) = YSCAL * YN3(JNP)
  140    CONTINUE
      END IF
      IF (ZSCAL .NE. 1. .AND. NDIM .EQ. 3) THEN
         DO 150 JNP = 1, NUMNP
            ZN3(JNP) = ZSCAL * ZN3(JNP)
  150    CONTINUE
      END IF

      RETURN
      END
