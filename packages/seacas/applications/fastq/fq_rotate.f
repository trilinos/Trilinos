C $Id: rotate.f,v 1.1 1990/11/30 11:15:08 gdsjaar Exp $
C $Log: rotate.f,v $
C Revision 1.1  1990/11/30 11:15:08  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]ROTATE.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE FQ_ROTATE (N, X, Y, NID, NEWF)
C***********************************************************************
C
C  SUBROUTINE ROTATE = CIRCULARLY SHIFTS THE DATA IN X,  Y,  AND NID
C
C***********************************************************************
C
      DIMENSION X (N), Y (N), NID (N)
C
      IF ((NEWF .LE. 1) .OR. (NEWF .GT. N)) RETURN
C
C  BUBBLE UP THROUGH THE ARRAYS AS MANY TIMES AS NEEDED
C
      DO 110 I = 1, NEWF - 1
         XLAST = X (1)
         YLAST = Y (1)
         NLAST = NID (1)
         DO 100 J = 1, N - 1
            X(J) = X (J + 1)
            Y(J) = Y (J + 1)
            NID(J) = NID (J + 1)
  100    CONTINUE
         X(N)   = XLAST
         Y(N)   = YLAST
         NID(N) = NLAST
  110 CONTINUE
C
      RETURN
C
      END
