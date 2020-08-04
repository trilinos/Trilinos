C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE NODE12 (MXND, MLN, LNODES, I1, I2, NLOOP1, NLOOP2,
     &   NODE1, NODE2, NODE, ERR)
C***********************************************************************

C  SUBROUTINE NODE12 = FINDS THE CURRENT NODE IN BOTH NEW LOOPS, AND
C                      KEEPS IT A CONSTANT

C***********************************************************************

      DIMENSION LNODES (MLN, MXND)

      LOGICAL ERR

      ERR = .FALSE.

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

  110 CONTINUE

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

  130 CONTINUE

      RETURN

      END
