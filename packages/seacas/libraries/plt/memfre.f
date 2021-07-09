C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION MEMFRE(IB,MEMRY)
      INTEGER IB
      INTEGER MEMRY(*)

      IPR = IB - 2
      MEMFRE = .FALSE.
      IP = 3
      IQ = 3
      IPN = MEMRY(IP)
 2710 IF (.NOT. (IPN.NE.0)) GO TO 2720
 2730 IF (.NOT. (IPN.LT.0)) GO TO 2740
      IPM = MEMRY(-IPN)
      IF (IPM.GT.0) THEN
         IQ = IP
         IP = -IPN
         IPN = IPM

      ELSE IF (IPM.EQ.0) THEN
         MEMRY(IP) = 0
         MEMRY(2) = IP
         RETURN

      ELSE IF (IPM.LT.0) THEN
         IPN = IPM
         MEMRY(IP) = IPN
      END IF

      GO TO 2730

 2740 CONTINUE
      IF (IP.EQ.IPR) THEN
         IPN = -IPN
         IF (MEMRY(IQ).LT.0) THEN
            IP = IQ
         END IF

         MEMRY(IP) = IPN
         MEMFRE = .TRUE.
         IB = 0

      ELSE
         IQ = IP
         IP = ABS(IPN)
         IPN = MEMRY(IP)
      END IF

      GO TO 2710

 2720 CONTINUE
      RETURN

      END
