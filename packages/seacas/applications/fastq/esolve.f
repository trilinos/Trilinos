C $Id: esolve.f,v 1.1 1990/11/30 11:06:51 gdsjaar Exp $
C $Log: esolve.f,v $
C Revision 1.1  1990/11/30 11:06:51  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]ESOLVE.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      REAL FUNCTION ESOLVE (A7, A8, A2, ANG1, ANG2)
C***********************************************************************
C
C  FUNCTION ESOLVE = FINDS A SOLUTION TO THE ELIPSE EQUATION
C                    GIVEN AN INTERVAL THAT CONTAINS THE SOLUTION
C
C***********************************************************************
C
      EPS = 1.E-6
C
      F1 = ELIPSE (A7, A8, A2, ANG1)
      IF (ABS(F1) .LT. EPS) THEN
         ESOLVE = ANG1
         GO TO 110
      END IF
      F2 = ELIPSE (A7, A8, A2, ANG2)
      IF (ABS(F2) .LT. EPS) THEN
         ESOLVE = ANG2
         GO TO 110
      END IF
C
  100 CONTINUE
      IF (ABS(ANG1 - ANG2) .LT. EPS) THEN
         ESOLVE = (ANG1 + ANG2)/2.0
      ELSE
         ANG3 = (ANG1 + ANG2)/2.0
         F3 = ELIPSE (A7, A8, A2, ANG3)
C
         IF (ABS(F3) .LT. EPS) THEN
            ESOLVE = ANG3
            GO TO 110
         ELSE IF (F1/F3 .LT. 0.0) THEN
            ANG2 = ANG3
            F2 = F3
         ELSE
            ANG1 = ANG3
            F1 = F3
         END IF
         GO TO 100
      END IF
C
  110 CONTINUE
C
      RETURN
C
      END
