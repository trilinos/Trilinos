C $Id: fndlin.f,v 1.1 1990/11/30 11:07:41 gdsjaar Exp $
C $Log: fndlin.f,v $
C Revision 1.1  1990/11/30 11:07:41  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]FNDLIN.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE FNDLIN (MXND, LXN, NODE1, NODE2, LINE, ERR)
C
C***********************************************************************
C
C  SUBROUTINE FNDLIN =  FINDS THE LINE WITH ENDS NODE1 & NODE2
C
C***********************************************************************
C
      DIMENSION LXN(4, MXND)
      DIMENSION LINES1(20), LINES2(20)
      LOGICAL ERR
C
      ERR = .FALSE.
C
      CALL GETLXN (MXND, LXN, NODE1, LINES1, NL1, ERR)
      CALL GETLXN (MXND, LXN, NODE2, LINES2, NL2, ERR)
C
      IF (.NOT.ERR) THEN
         ERR = .TRUE.
         DO 110 I = 1, NL1
            DO 100 J = 1, NL2
               IF (LINES1(I) .EQ. LINES2(J)) THEN
                  LINE = LINES1(I)
                  ERR = .FALSE.
                  RETURN
               END IF
  100       CONTINUE
  110    CONTINUE
      END IF
C
      RETURN
C
      END
