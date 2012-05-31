C $Id: setlop.f,v 1.1 1990/11/30 11:15:27 gdsjaar Exp $
C $Log: setlop.f,v $
C Revision 1.1  1990/11/30 11:15:27  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]SETLOP.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE SETLOP (MXND, MLN, NLOOP, LNODES, NODE, IVALUE, ERR)
C***********************************************************************
C
C  SUBROUTINE SETLOP = MARKS ALL THE NODES IN THE LOOP AS DESIGNATED
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
      INOW = NODE
C
  100 CONTINUE
      INOW = LNODES (3, INOW)
      LNODES (1, INOW) = IVALUE
C
      IF (INOW .EQ. NODE) RETURN
C
      KOUNT = KOUNT + 1
      IF (KOUNT .GT. NLOOP) THEN
         CALL MESAGE('PROBLEMS IN SETLOP WITH LOOP NOT CLOSING')
         ERR = .TRUE.
         GOTO 110
      ENDIF
      GOTO 100
C
  110 CONTINUE
C
      RETURN
C
      END
