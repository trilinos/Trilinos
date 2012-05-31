C $Id: bubble.f,v 1.1 1990/11/30 11:04:13 gdsjaar Exp $
C $Log: bubble.f,v $
C Revision 1.1  1990/11/30 11:04:13  gdsjaar
C Initial revision
C
CC* FILE: [.QMESH]BUBBLE.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE BUBBLE (X, KARRY, NORD, N)
C***********************************************************************
C
C  SUBROUTINE BUBBLE=SORTS ALL VALUES X(I), KARRY(I) INTO DECREASING
C                      ORDER, ASSUMING THAT VALUES 1 TO NORD ARE SORTED
C
C***********************************************************************
C
      DIMENSION X (N), KARRY (N)
C
      IF (N .LE. 1) RETURN
C
      ISTART = MAX0 (NORD + 1, 2)
      IF (ISTART .GT. N) RETURN
      DO 120 J = ISTART, N
         XVAL = X (J)
         KVAL = KARRY (J)
         JM1 = J - 1
         I = J
         DO 100 II = 1, JM1
            IF  (XVAL .LE. X (I - 1)) GO TO 110
            X (I) = X (I - 1)
            KARRY (I) = KARRY (I - 1)
            I = I - 1
  100    CONTINUE
  110    CONTINUE
         X (I) = XVAL
         KARRY (I) = KVAL
  120 CONTINUE
C
      RETURN
C
      END
