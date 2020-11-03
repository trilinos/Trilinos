C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SHWINT (ITRANT, NEREPL, DIM3, NRTRAN, D3TRAN, ZGRAD)
C=======================================================================

      INTEGER NRTRAN(*)
      REAL    D3TRAN(*), ZGRAD(*)

      CHARACTER*20 RSTR(9)
      CHARACTER*20 STRA, TYPE

      CALL INTSTR (1, 0, NEREPL, STRA, LSTRA)
      CALL NUMSTR (1, 4, DIM3, RSTR(1), LR)

      IF      (ITRANT .EQ.  0) THEN
         TYPE = 'Transform'
      ELSE IF (ITRANT .EQ.  1) THEN
         TYPE = 'Translate'
      ELSE IF (ITRANT .EQ.  2) THEN
         TYPE = 'Rotate'
      ELSE IF (ITRANT .EQ.  4) THEN
         TYPE = 'Warp'
      ELSE IF (ITRANT .EQ.  8) THEN
         TYPE = 'Twist'
      ELSE IF (ITRANT .EQ. 16) THEN
         TYPE = 'Project'
      ELSE IF (ITRANT .EQ. 32) THEN
         TYPE = 'ExpRotate'
      ELSE IF (ITRANT .EQ. 64) THEN
         TYPE = 'Spline'
      ELSE
         CALL PRTERR ('PROGRAM', 'Unknown transformation option')
         RETURN
      END IF

      LT = LENSTR(TYPE)

      IF (NEREPL .EQ. NRTRAN(1)) THEN
         IF (ABS (ZGRAD(1) - 1.0) .LE. 1.0E-6) THEN
            WRITE (*, 20) TYPE(:LT), ' mesh ', STRA(:LSTRA),
     &         ' times for a total of ', RSTR(1)(:LR)
         ELSE
            CALL NUMSTR (1, 3, ZGRAD(1), RSTR(2), LR2)
            WRITE (*, 20) TYPE(:LT), ' mesh ', STRA(:LSTRA),
     &         ' times for a total of ', RSTR(1)(:LR),
     &         ' with a gradient of ', RSTR(2)(:LR2)
         END IF
      ELSE
         IBLK = 0
         NR = 0
   10    CONTINUE
         IF (.TRUE.) THEN
            IBLK = IBLK + 1
            NR = NR + NRTRAN(IBLK)
            CALL INTSTR (1, 0, NRTRAN(IBLK), STRA, LSTRA)
            CALL NUMSTR (1, 4, D3TRAN(IBLK), RSTR(1), LR)
            IF (ABS (ZGRAD(IBLK) - 1.0) .LE. 0.001) THEN
               WRITE (*, 20) TYPE(:LT), ' mesh ',
     &            STRA(:LSTRA), ' times for a subtotal of ',
     &            RSTR(1)(:LR)
            ELSE
               CALL NUMSTR (1, 3, ZGRAD(IBLK), RSTR(2), LR2)
               WRITE (*, 20) TYPE(:LT), ' mesh ',
     &            STRA(:LSTRA), ' times for a subtotal of ',
     &            RSTR(1)(:LR), ' with a gradient of ',
     &            RSTR(2)(:LR2)
            END IF
            IF (NR .LT. NEREPL) GOTO 10
         END IF
      END IF
   20 FORMAT (1X, 10A)
      RETURN
      END
