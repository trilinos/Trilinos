C $Id: node12.f,v 1.1 1990/11/30 11:12:44 gdsjaar Exp $
C $Log: node12.f,v $
C Revision 1.1  1990/11/30 11:12:44  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]NODE12.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE NODE12 (MXND, MLN, LNODES, I1, I2, NLOOP1, NLOOP2,
     &   NODE1, NODE2, NODE, ERR)
C***********************************************************************
C
C  SUBROUTINE NODE12 = FINDS THE CURRENT NODE IN BOTH NEW LOOPS, AND
C                      KEEPS IT A CONSTANT
C
C***********************************************************************
C
      DIMENSION LNODES (MLN, MXND)
C
      LOGICAL ERR
C
      ERR = .FALSE.
C
      KOUNT = 0
      NTEST = I1
  100 CONTINUE
      KOUNT = KOUNT + 1
      IF (KOUNT .GT. NLOOP2) GOTO 110
      IF (NTEST .EQ. NODE) THEN
         NODE1 = NODE
         NODE2 = I2
         GOTO 130
      ENDIF
      NTEST = LNODES (3, NTEST)
      GOTO 100
C
  110 CONTINUE
C
      KOUNT = 0
      NTEST = I2
  120 CONTINUE
      KOUNT = KOUNT + 1
      IF (KOUNT .GT. NLOOP1) THEN
         CALL MESAGE ('** PROBLEMS IN NODE12 FINDING NODE **')
         ERR = .TRUE.
         GOTO 130
      ENDIF
      IF (NTEST .EQ. NODE) THEN
         NODE1 = I1
         NODE2 = NODE
         GOTO 130
      ENDIF
      NTEST = LNODES (3, NTEST)
      GOTO 120
C
  130 CONTINUE
C
      RETURN
C
      END
