C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE BCROSS (MXND, MLN, XN, YN, ZN, LXK, KXL, NXL, LXN,
     &   LNODES, I1, I2, J1, J2, NLOOP, BOK, LLL, XMIN, XMAX, YMIN,
     &   YMAX, ZMIN, ZMAX, DEV1, KREG, ERR)
C***********************************************************************

C  SUBROUTINE BCROSS = CHECKS TO MAKE SURE THAT A BOUNDARY IS NOT
C                      BECOMING A PERMANENT CROSS

C***********************************************************************

      DIMENSION XN(MXND), YN(MXND), ZN(MXND)
      DIMENSION NXL(2, 3*MXND), LXN(4, MXND)
      DIMENSION LXK(4, MXND), KXL(2, 3*MXND)
      DIMENSION LNODES(MLN, MXND)

      LOGICAL BOK, LCROSS, ERR

      CHARACTER*3 DEV1

      BOK = .TRUE.
      ERR = .FALSE.

      J0 = LNODES(2, J1)
      J3 = LNODES(3, J2)

C  IF J0 TO I2, OR J3 TO I1 IS TO BECOME A BOUNDARY LINE,
C  THEN TEST THE CONNECTION TO SEE IF IT INTERSECTS ANY OTHER
C  BOUNDARY LINES

      KOUNT = 0

C  FIRST TEST THE J0 TO I2 LINE

      IF ((LXN(2, J0) .LT. 0) .AND. (LXN(2, I2) .LT. 0)) THEN
         NODE1 = I1
         NODE2 = I2
  100    CONTINUE
         NODE1 = NODE2
         NODE2 = LNODES(3, NODE2)
         IF ((LXN(2, NODE1) .LT. 0) .AND. (LXN(2, NODE2) .LT. 0)) THEN
            IF (NODE2 .EQ. J0) THEN
               CALL GETANG (MXND, MLN, XN, YN, LNODES, LXK, KXL,
     &            NXL, LXN, NODE1, J0, J1, ANGLE, ERR)
               IF (ANGLE .LT. 0) THEN
                  BOK = .FALSE.
                  GOTO 130
               ELSE
                  GOTO 110
               ENDIF

            ELSE
               CALL INTSCT (XN(NODE1), YN(NODE1), XN(NODE2), YN(NODE2),
     &            XN(J0), YN(J0), XN(I2), YN(I2), U, W, LCROSS)
               IF (LCROSS) THEN
                  BOK = .FALSE.
                  GOTO 130
               ENDIF
            ENDIF
         ENDIF
         KOUNT = KOUNT + 1
         IF (KOUNT .LT. NLOOP) THEN
            GOTO 100
         ELSE
            ERR = .TRUE.
            GOTO 130
         ENDIF
      ENDIF

  110 CONTINUE

C  NEXT TEST THE J3 TO I1 LINE

      KOUNT = 0
      IF ((LXN(2, J3) .LT. 0) .AND. (LXN(2, I1) .LT. 0)) THEN
         NODE1 = J3
         NODE2 = LNODES(3, J3)
  120    CONTINUE
         NODE1 = NODE2
         NODE2 = LNODES(3, NODE2)
         IF ((LXN(2, NODE1) .LT. 0) .AND. (LXN(2, NODE2) .LT. 0)) THEN
            IF (NODE2 .EQ. I1) THEN
               CALL GETANG (MXND, MLN, XN, YN, LNODES, LXK, KXL,
     &            NXL, LXN, NODE1, I1, I2, ANGLE, ERR)
               IF (ANGLE .LT. 0) THEN
                  BOK = .FALSE.
                  GOTO 130
               ELSE
                  GOTO 130
               ENDIF

            ELSE
               CALL INTSCT (XN(NODE1), YN(NODE1), XN(NODE2), YN(NODE2),
     &            XN(J3), YN(J3), XN(I1), YN(I1), U, W, LCROSS)
               IF (LCROSS) THEN
                  BOK = .FALSE.
                  GOTO 130
               ENDIF
            ENDIF
         ENDIF
         KOUNT = KOUNT + 1
         IF (KOUNT .LT. NLOOP) THEN
            GOTO 120
         ELSE
            ERR = .TRUE.
            GOTO 130
         ENDIF
      ENDIF

  130 CONTINUE

      RETURN

      END
