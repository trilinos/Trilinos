C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
      SUBROUTINE SHFTC(IARAY,N,ISTRT,IEND,INT)

C     THIS SUBROUTINE SHIFTS THE ROWS IN A CHARACTER ARRAY. IF INT>0
C     THEN ALL ROWS ISTRT TO IEND ARE SHIFTED UP "INT" ROWS. IF INT<0
C     THEN ALL ROWS ISTRT TO IEND ARE SHIFTED DOWN "INT" ROWS.

C     MODIFIED JAN 1986 BY BILL MILLS-CURRAN
C     THE CHARACTER ARRAY IS NOW 1-D RATHER THAN 2-D.

      CHARACTER*(*) IARAY(N)

C  CALCULATE RANGE AND INCREMENT OF DO LOOP

      IF(INT.LT.0)THEN

C  SHIFT DOWN

         I1=IEND
         I2=ISTRT
         ID=-1
      ELSE IF(INT.GT.0)THEN

C  SHIFT UP

         I1=ISTRT
         I2=IEND
         ID=1
      ELSE
         RETURN
      END IF

C  PERFORM SHIFT

      DO I=I1,I2,ID
         IARAY(I-INT)=IARAY(I)
      END DO
      RETURN
      END
