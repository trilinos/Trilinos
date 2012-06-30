C $Id: sortia.f,v 1.1 1990/11/30 11:16:11 gdsjaar Exp $
C $Log: sortia.f,v $
C Revision 1.1  1990/11/30 11:16:11  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]SORTIA.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE SORTIA(N,IAR,IRANGE,I)
C***********************************************************************
C
C  SUBROUTINE SORTIA = THIS SUBROUTINE SORTS AN INTEGER ARRAY IN
C                      ASCENDING ORDER
C
C***********************************************************************
C
C VARIABLES   IN : N ...   NUMBER OF ELEMENTS IN THE ARRAY
C                  IAR ... INTEGER ARRAY WITH DATA TO BE SORTED
C             OUT: I ...   INDEX ARRAY WITH ITS SORTED VALUES IN
C                          ASCENDING ORDER.
C
C WRITTEN BY:  HORACIO RECALDE                 DATE: FEB 25,  1988
C***********************************************************************
C
      INTEGER I(N),IAR(N)
C
C-- COPY ELEMENTS IN THE I ARRAY
C
      DO 100 J = 1,IRANGE
         I(J) = IAR(J)
  100 CONTINUE
C
C---  PERFORM AN EXCHANGE SORT ON THE FIRST IRANGE-1
C
      DO 120 K = 1,IRANGE - 1
         MIN = I(K)
C
C---  EXCHANGE THE K-TH ELEMENTS WITH THE MINIMUM ELEMENT REMAIN
C
         DO 110 J = K+1,IRANGE
            IF (I(J) .LT. MIN) THEN
               L = I(J)
               I(J) = I(K)
               I(K) = L
               MIN = I(K)
            ENDIF
  110    CONTINUE
  120 CONTINUE
C
      RETURN
      END
