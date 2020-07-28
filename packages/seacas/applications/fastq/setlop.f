C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE SETLOP (MXND, MLN, NLOOP, LNODES, NODE, IVALUE, ERR)
C***********************************************************************

C  SUBROUTINE SETLOP = MARKS ALL THE NODES IN THE LOOP AS DESIGNATED

C***********************************************************************

      DIMENSION LNODES (MLN, MXND)

      LOGICAL ERR

      ERR = .FALSE.

      KOUNT = 0
      INOW = NODE

  100 CONTINUE
      INOW = LNODES (3, INOW)
      LNODES (1, INOW) = IVALUE

      IF (INOW .EQ. NODE) RETURN

      KOUNT = KOUNT + 1
      IF (KOUNT .GT. NLOOP) THEN
         CALL MESAGE('PROBLEMS IN SETLOP WITH LOOP NOT CLOSING')
         ERR = .TRUE.
         GOTO 110
      ENDIF
      GOTO 100

  110 CONTINUE

      RETURN

      END
