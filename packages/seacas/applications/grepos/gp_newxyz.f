C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE NEWXYZ (XN, YN, ZN, NUMNP, NDIM, A)
C=======================================================================

C     --*** NEWXYZ *** (GEND) Modify coordinates
C     --   Written by Amy Gilkey - revised 05/09/88
C     --   Modified by Greg Sjaardema - 02/06/89
C     --
C     --NEWXYZ modifies the coordinate array for the database.
C     --
C     --Parameters:
C     --   XN, YN, ZN - OUT - the coordinates
C     --   NUMNP - IN - Number of nodes
C     --   NDIM  - IN - Number of spatial dimensions
C     --   A - IN - Base array for memory allocation
C     --
C     --Common Variables:
C     --   Uses XOFFS, YOFFS, ZOFFS, SPLOFF of /XYZOFF/
C     --   Uses ROT3D, ROTMAT of /XYZROT/

      include 'gp_xyzoff.blk'
      include 'gp_xyzrot.blk'
      include 'gp_xyzmir.blk'
      include 'gp_xyzero.blk'
      include 'gp_xyzscl.blk'
      include 'gp_splxyz.blk'

      REAL XN(numnp), YN(numnp), ZN(numnp)
      REAL A(*)

C     --Rotate 3D mesh, if needed

      IF (ROT3D .AND. NDIM .EQ. 3) THEN
         DO 10 JNP = 1, NUMNP
            X = XN(JNP) - ROTCEN(1)
            Y = YN(JNP) - ROTCEN(2)
            Z = ZN(JNP) - ROTCEN(3)
            XN(JNP) = X*ROTMAT(1,1) + Y*ROTMAT(2,1) + Z*ROTMAT(3,1)
     &           + ROTCEN(1)
            YN(JNP) = X*ROTMAT(1,2) + Y*ROTMAT(2,2) + Z*ROTMAT(3,2)
     &           + ROTCEN(2)
            ZN(JNP) = X*ROTMAT(1,3) + Y*ROTMAT(2,3) + Z*ROTMAT(3,3)
     &           + ROTCEN(3)
   10    CONTINUE
      ELSE IF (ROT3D .AND. NDIM .EQ. 2) THEN
         DO 20 JNP = 1, NUMNP
            X = XN(JNP) - ROTCEN(1)
            Y = YN(JNP) - ROTCEN(2)
            XN(JNP) = X*ROTMAT(1,1) + Y*ROTMAT(2,1) + ROTCEN(1)
            YN(JNP) = X*ROTMAT(1,2) + Y*ROTMAT(2,2) + ROTCEN(2)
   20    CONTINUE
      END IF

C     --Add offset, if any, to coordinates

      IF (SPLOFF) THEN
C     ... Finish processing in sub-routine

         CALL SPOFF (XN, YN, ZN, NSPL, A(KZSPL), A(KXSPL), A(KXSPL2),
     &        A(KYSPL), A(KYSPL2), SLTOP, SLBOT, NUMNP, NDIM)

      ELSE
         IF (XOFFS .NE. 0.0) THEN
            DO 30 JNP = 1, NUMNP
               XN(JNP) = XN(JNP) + XOFFS
   30       CONTINUE
         END IF
         IF (YOFFS .NE. 0.0) THEN
            DO 40 JNP = 1, NUMNP
               YN(JNP) = YN(JNP) + YOFFS
   40       CONTINUE
         END IF
         IF (ZOFFS .NE. 0.0 .AND. NDIM .EQ. 3) THEN
            DO 50 JNP = 1, NUMNP
               ZN(JNP) = ZN(JNP) + ZOFFS
   50       CONTINUE
         END IF
      END IF

C     --Mirror coordinates if any specified

      IF (XMIRR .LT. 0.) THEN
         DO 60 JNP = 1, NUMNP
            XN(JNP) = -XN(JNP)
   60    CONTINUE
      END IF
      IF (YMIRR .LT. 0.) THEN
         DO 70 JNP = 1, NUMNP
            YN(JNP) = -YN(JNP)
   70    CONTINUE
      END IF
      IF (ZMIRR .LT. 0. .AND. NDIM .EQ. 3) THEN
         DO 80 JNP = 1, NUMNP
            ZN(JNP) = -ZN(JNP)
   80    CONTINUE
      END IF

C     --Randomize coordinates if any specified

      IDUM = 1
      IF (XRAND .NE. 0.) THEN
         DO 61 JNP = 1, NUMNP
            XN(JNP) = (2.0*RAN1(IDUM)-1.0) * XRAND + XN(JNP)
 61      CONTINUE
      END IF
      IF (YRAND .NE. 0.) THEN
         DO 71 JNP = 1, NUMNP
            YN(JNP) = (2.0*RAN1(IDUM)-1.0) * YRAND + YN(JNP)
 71      CONTINUE
      END IF
      IF (ZRAND .NE. 0. .AND. NDIM .EQ. 3) THEN
         DO 81 JNP = 1, NUMNP
            ZN(JNP) = (2.0*RAN1(IDUM)-1.0) * ZRAND + ZN(JNP)
 81      CONTINUE
      END IF

C     --- Zero coordinates if *ZERO is not equal to zero

      IF (XZERO .NE. 0.) THEN
         DO 90 JNP = 1, NUMNP
            IF (ABS(XN(JNP)) .LE. XZERO) XN(JNP) = 0.0
   90    CONTINUE
      END IF
      IF (YZERO .NE. 0.) THEN
         DO 100 JNP = 1, NUMNP
            IF (ABS(YN(JNP)) .LE. YZERO) YN(JNP) = 0.0
  100    CONTINUE
      END IF
      IF (ZZERO .NE. 0. .AND. NDIM .EQ. 3) THEN
         DO 110 JNP = 1, NUMNP
            IF (ABS(ZN(JNP)) .LE. ZZERO) ZN(JNP) = 0.0
  110    CONTINUE
      END IF

C     --- Scale the coordinates if any Scaled

      IF (XSCAL .NE. 1.) THEN
         DO 120 JNP = 1, NUMNP
            XN(JNP) = XSCAL * XN(JNP)
  120    CONTINUE
      END IF
      IF (YSCAL .NE. 1.) THEN
         DO 130 JNP = 1, NUMNP
            YN(JNP) = YSCAL * YN(JNP)
  130    CONTINUE
      END IF
      IF (ZSCAL .NE. 1. .AND. NDIM .EQ. 3) THEN
         DO 140 JNP = 1, NUMNP
            ZN(JNP) = ZSCAL * ZN(JNP)
  140    CONTINUE
      END IF

      RETURN
      END
