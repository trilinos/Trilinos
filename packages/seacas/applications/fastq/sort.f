C $Id: sort.f,v 1.1 1990/11/30 11:16:08 gdsjaar Exp $
C $Log: sort.f,v $
C Revision 1.1  1990/11/30 11:16:08  gdsjaar
C Initial revision
C
C
CC* FILE: [.RENUM]SORT.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE SORT (N, IX, IY)
C***********************************************************************
C
C  SUBROUTINE SORT = SORT THE ARRAY IX,  CARRYING ALONG IY
C
C***********************************************************************
C
      DIMENSION IX (N), IY (N)
      NN = N
      M = NN
  100 CONTINUE
      M =  (9 * M) / 16
      IF (M .LE. 0) RETURN
      M1 = M + 1
      DO 120 J = M1, NN
         L = J
         I = J - M
  110    CONTINUE
         IF (IX (L) .LT. IX (I)) THEN
            KEEPX = IX (I)
            KEEPY = IY (I)
            IX (I) = IX (L)
            IY (I) = IY (L)
            IX (L) = KEEPX
            IY (L) = KEEPY
            L = I
            I = I - M
            IF (I .GE. 1)GOTO 110
         ENDIF
  120 CONTINUE
C
      GOTO 100
C
      END
