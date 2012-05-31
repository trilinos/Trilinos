C $Id: revers.f,v 1.1 1990/11/30 11:14:53 gdsjaar Exp $
C $Log: revers.f,v $
C Revision 1.1  1990/11/30 11:14:53  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]REVERS.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE REVERS (X, N)
C***********************************************************************
C
C  SUBROUTINE REVERS = REVERS THE REAL ARRAY OF X (I), I=1, N
C
C***********************************************************************
C
      DIMENSION X (N)
C
      IF (N .LE. 1) RETURN
C
      NUP = N + 1
      M = N / 2
      DO 100 I = 1, M
         NUP = NUP - 1
         XK = X (I)
         X (I) = X (NUP)
         X (NUP) = XK
  100 CONTINUE
C
      RETURN
C
      END
