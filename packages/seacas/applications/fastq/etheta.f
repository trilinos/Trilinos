C $Id: etheta.f,v 1.2 1998/07/14 18:18:44 gdsjaar Exp $
C $Log: etheta.f,v $
C Revision 1.2  1998/07/14 18:18:44  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.1.1.1  1990/11/30 11:06:56  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:06:53  gdsjaar
c Initial revision
c 
C
CC* FILE: [.MAIN]ETHETA.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE ETHETA (A7, A8, A2, THETA1, THETA2, THETA, ERR)
C***********************************************************************
C
C  SUBROUTINE ETHETA = ITERATIVELY SOLVES THE ELIPTICAL PROBLEM OF
C                      FINDING AN "A" DISTANCE GIVEN TWO POINTS ON
C                      THE ELIPSE AND A CENTER POINT
C
C***********************************************************************
C
      LOGICAL ERR
C
C  START WITH 20 INCREMENTS, EACH PASS INCREMENTS DECREASE TEN FOLD
C
      ANGINC = (THETA2 - THETA1) * .05
      ANG = THETA1
      F1 = ELIPSE (A7, A8, A2, ANG)
      IF (F1 .EQ. 0.0) THEN
         THETA = ANG
         GO TO 110
      END IF
      ANG2 = ANG + ANGINC
  100 CONTINUE
      F2 = ELIPSE (A7, A8, A2, ANG2)
      IF (F2 .EQ. 0.0) THEN
         THETA = ANG2
      ELSE IF (F1*F2 .LT. 0.0) THEN
         THETA = ESOLVE (A7, A8, A2, ANG, ANG2)
      ELSE
         ANG = ANG2
         ANG2 = ANG2 + ANGINC
         IF (ANG2 .LE. THETA2) GO TO 100
         ERR = .TRUE.
         GO TO 120
      END IF
C
  110 CONTINUE
C
  120 CONTINUE
C
C  FIND THE SECOND ROOT IF THE FIRST ONE HAS BEEN LOCATED
C
      IF(.NOT.ERR)THEN
         ANG=THETA+ANGINC
         F1 = ELIPSE (A7, A8, A2, ANG)
         F2 = ELIPSE (A7, A8, A2, THETA2)
         IF (F1 .EQ. 0.0) THEN
            THETA = ANG
         ELSEIF (F2 .EQ. 0.0) THEN
            THETA = THETA2
         ELSE IF (F1*F2 .LT. 0.0) THEN
            THETA = ESOLVE (A7, A8, A2, ANG, THETA2)
         ELSE
            GO TO 130
         END IF
      END IF
C
  130 CONTINUE
      RETURN
C
      END
