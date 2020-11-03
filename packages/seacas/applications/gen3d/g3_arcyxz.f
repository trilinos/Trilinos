C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE ARCYXZ (XN, YN, XN3, YN3, ZN3, IXNP, NRNP, NPCEN,
     &   SINANG, COSANG)
C=======================================================================

C   --*** ARCYXZ *** (GEN3D) Calculate 3D coordinates - rotation about X
C   --   Written by Amy Gilkey - revised 05/09/88
C   --   Modified by Greg Sjaardema - 02/06/89
C   --       Added Warp Function
C   --       Added Gradient to Rotations (not for center blocks)
C   --
C   --NEWXYZ calculates the coordinate array for the 3D database.
C   --
C   --Parameters:
C   --   XN, YN - IN - the 2D coordinates, destroyed
C   --   XN3, YN3, ZN3 - OUT - the 3D coordinates
C   --   IXNP - IN - the new index for each node
C   --   NRNP - IN - the number of new nodes generated for each node
C   --   NPCEN - IN - the node numbers of the center nodes by column and row
C   --   SINANG, COSANG - SCRATCH - size = NNREPL, holds sin and cos of
C   --      angles for rotations
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
C .. Dimension of NPCEN is (NUMCOL,NUMROW), but NUMCOL may be 0
      INTEGER NPCEN(NUMCDM,*)
      REAL SINANG(NNREPL), COSANG(NNREPL)

      PI = ATAN2(0.0, -1.0)
C      --For rotations, change X and add Z so that they define pie-shaped
C      --regions around the Y axis

C   --Copy X coordinate from original

         DO 520 INP = 1, NUMNP
            JNP0 = IXNP(INP) - 1
            DO 510 NR = 1, NRNP(INP)
               XN3(JNP0+NR) = XN(INP)
  510       CONTINUE
  520    CONTINUE

C ... Check minimum Y coordinate to see if we will generate bad elements
         YMIN = YN(1)
         DO 5 INP = 2, NUMNP
           YMIN = MIN(YMIN, YN(INP))
 5       CONTINUE
         YMIN = YMIN - CENTER
         if (YMIN .LT. 0) THEN
           call prterr('WARNING',
     &       'Input mesh crosses over axis of rotation (X). Negative '//
     &       'elements may be generated. Adjust ROTCEN or input mesh.')
         end if

C      --If center block exists, calculate center of rotation

      IF (NUMCOL .GT. 0) THEN
         CENTER = 1E36
         DO 10 IROW = 1, NUMROW
            INP = NPCEN(1,IROW)
            IF (INP .LE. 0) GOTO 10
            CENTER = MIN (CENTER, YN(INP))
   10    CONTINUE
      END IF

C      --Subtract center, so rotation is around zero

      IF (CENTER .NE. 0.0) THEN
         DO 20 INP = 1, NUMNP
            YN(INP) = YN(INP) - CENTER
   20    CONTINUE
      END IF

C      --Get the coordinates for the non-center nodes

      IBLK = 0
      NXTNR = 1
      ZEND = 0.0

    1    CONTINUE
         IBLK = IBLK + 1
         IF (NRTRAN(IBLK) .GT. 0) THEN
            ZBEG = ZEND
            ZEND = ZBEG + D3TRAN(IBLK) * PI / 180.0

            CALL INIGRD (-ZBEG, -ZEND, ZGRAD(IBLK),
     *         NRTRAN(IBLK), MIN(NRTRAN(IBLK)+1,NNREPL-NXTNR+1),
     *        COSANG(NXTNR))
            NXTNR = NXTNR + NRTRAN(IBLK)
            IF (IBLK .LT. MAXINT) GO TO 1
         END IF

         DO 30 NR = 1, NNREPL
            ANG = COSANG(NR)
            SINANG(NR) = SIN (ANG)
            COSANG(NR) = COS (ANG)
   30    CONTINUE

      DO 60 INP = 1, NUMNP
         IF (NRNP(INP) .EQ. NNREPL) THEN
            JNP = IXNP(INP)
            Y = YN(INP)
            DO 50 NR = 1, NNREPL
               YN3(JNP) = Y * COSANG(NR)
               ZN3(JNP) = Y * SINANG(NR)
               JNP = JNP + 1
   50       CONTINUE
         END IF
   60 CONTINUE

      IF (NUMCOL .GT. 0) THEN

C         --For center rotations, process by column for each row;
C         --calculate the first 8th of the circle, then copy to other 8ths.

         RAD45 = 45.0 * (PI / 180.0)
         COSD45 = COS (RAD45)

         DO 70 IROW = 1, NUMROW
            INP = NPCEN(1,IROW)
            IF (INP .LE. 0) GOTO 70
            JNP = IXNP(INP)
            YN3(JNP) = 0.0
            ZN3(JNP) = 0.0
   70    CONTINUE

         DO 170 ICOL = 2, NUMCOL

            DO 160 IROW = 1, NUMROW
               INP = NPCEN(ICOL,IROW)
               IF (INP .LE. 0) GOTO 160

               JNP = IXNP(INP)
               JNP0 = JNP - 1
               NEND = JNP0 + NRNP(INP)
               YN3(JNP) = YN(INP)
               ZN3(JNP) = 0.0

               Y = YN(INP)
               YMID = Y * COSD45

               DO 80 IX = 2, ICOL-1
                  I = NPCEN(IX,IROW)
                  IF (I .NE. 0) THEN
                     Z = YN(ABS(I))
                  ELSE
                     Z = YN(INP) * (IX-1) / (ICOL-1)
                  END IF
                  IF (Y .NE. 0.0) THEN
                     R = Z / Y
                  ELSE
                     R = 0.0
                  END IF
                  JNP = JNP + 1
                  YN3(JNP) = SQRT (Y*Y - Z*Z/2)
                  ZN3(JNP) = - (R * YMID)
   80          CONTINUE

               JNP = JNP + 1
               YN3(JNP) = YMID
               ZN3(JNP) = - YMID

               DO 90 IX = ICOL-1, 1, -1
                  JNP = JNP + 1
                  YN3(JNP) = - ZN3(JNP0+IX)
                  ZN3(JNP) = - YN3(JNP0+IX)
   90          CONTINUE
               IF (JNP .LT. NEND) THEN
                  DO 100 IX = 2, ICOL
                     JNP = JNP + 1
                     YN3(JNP) = ZN3(JNP0+IX)
                     ZN3(JNP) = - YN3(JNP0+IX)
  100             CONTINUE
                  DO 110 IX = ICOL-1, 1, -1
                     JNP = JNP + 1
                     YN3(JNP) = - YN3(JNP0+IX)
                     ZN3(JNP) = ZN3(JNP0+IX)
  110             CONTINUE
               END IF
               IF (JNP .LT. NEND) THEN
                  DO 120 IX = 2, ICOL
                     JNP = JNP + 1
                     YN3(JNP) = - YN3(JNP0+IX)
                     ZN3(JNP) = - ZN3(JNP0+IX)
  120             CONTINUE
                  DO 130 IX = ICOL-1, 1, -1
                     JNP = JNP + 1
                     YN3(JNP) = ZN3(JNP0+IX)
                     ZN3(JNP) = YN3(JNP0+IX)
  130             CONTINUE
               END IF
               IF (JNP .LT. NEND) THEN
                  DO 140 IX = 2, ICOL
                     JNP = JNP + 1
                     YN3(JNP) = - ZN3(JNP0+IX)
                     ZN3(JNP) = YN3(JNP0+IX)
  140             CONTINUE
                  DO 150 IX = ICOL-1, 2, -1
                     JNP = JNP + 1
                     YN3(JNP) = YN3(JNP0+IX)
                     ZN3(JNP) = - ZN3(JNP0+IX)
  150             CONTINUE
               END IF
  160       CONTINUE
  170    CONTINUE
      END IF

C      --Add center of rotation

      IF (CENTER .NE. 0.0) THEN
         DO 180 JNP = 1, NUMNP3
            YN3(JNP) = YN3(JNP) + CENTER
  180    CONTINUE
      END IF
      RETURN
      END
