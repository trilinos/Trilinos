C $Id: setcir.f,v 1.1 1990/11/30 11:15:24 gdsjaar Exp $
C $Log: setcir.f,v $
C Revision 1.1  1990/11/30 11:15:24  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]SETCIR.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE SETCIR (MXND, MLN, NLOOP, LNODES, NODE, ERR)
C***********************************************************************
C
C  SUBROUTINE SETCIR = MARKS ALL THE NODES IN THE CIRCULAR LOOP
C  AS SIDES EXCEPT FOR ROW CORNERS
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
      NEWNOD = NODE
C
  100 CONTINUE
      INOW = LNODES (3, INOW)
      IF (LNODES (1, INOW) .LE. 4) THEN
         LNODES (1, INOW) = 3
      ELSEIF (LNODES (1, INOW) .EQ. 6) THEN
         LNODES (1, INOW) = 5
      ENDIF
      IF ( (LNODES (1, NEWNOD) .EQ. 3) .AND.
     &   ( (LNODES (1, INOW) .EQ. 5) .OR. (LNODES (1, INOW) .EQ. 7)) )
     &   THEN
         NEWNOD = INOW
      ENDIF
C
      IF (INOW .EQ. NODE) THEN
         NODE = NEWNOD
         RETURN
      ENDIF
C
      KOUNT = KOUNT + 1
      IF (KOUNT .GT. NLOOP) THEN
         CALL MESAGE('PROBLEMS IN SETCIR WITH LOOP NOT CLOSING')
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
