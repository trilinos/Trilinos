C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      INTEGER FUNCTION MEMALL(LENGTH,MEMRY)
      INTEGER LENGTH
      INTEGER MEMRY(*)
      INTEGER LR
      INTEGER LMX
      INTEGER IP
      INTEGER IPN
      INTEGER LAV
      INTEGER LT
      INTEGER IB

      IF (LENGTH.LE.0) THEN
         WRITE (6,*) ' Cannot allocate a segment of length zero.'
      END IF

      LR = LENGTH + MOD(LENGTH,2)
      LMX = MEMRY(1)
      IB = 0
      IP = 3
 2650 IF (.NOT. (IB.EQ.0)) GO TO 2660
      IPN = MEMRY(IP)
      IF (IPN.LE.0) THEN
         LT = IP + LR + 1
         IF (LT.GE.LMX) THEN
            WRITE (6,*) ' Cannot allocate space in memall.'
            WRITE (6,*) ' lt >= lmx '
            WRITE (6,*) ' lt,lr,lmx,length: ',LT,LR,LMX,LENGTH
            MEMALL = (-1)
            RETURN

         END IF

         IF (IPN.EQ.0) THEN
            LAV = LMX - IP - 1

         ELSE IF (IPN.LT.0) THEN
            LAV = -IPN - IP - 2
         END IF

         IF (LAV.GT.LR) THEN
            MEMRY(IP) = LT + 1
            MEMRY(LT+1) = IPN
            IF (IPN.EQ.0) THEN
               MEMRY(2) = LT + 1
            END IF

            IB = IP + 2
            DO 2670 J = 1,LR
               MEMRY(IP+1+J) = 0
 2670       CONTINUE

         ELSE IF (LAV.EQ.LR) THEN
            MEMRY(IP) = -IPN
            IB = IP + 2
            DO 2690 J = 1,LR
               MEMRY(IP+1+J) = 0
 2690       CONTINUE
         END IF

      END IF

      IP = ABS(IPN)
      GO TO 2650

 2660 CONTINUE
      MEMALL = (IB)
      RETURN

      END
