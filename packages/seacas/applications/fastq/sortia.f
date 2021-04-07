C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE SORTIA(N,IAR,IRANGE,I)
C***********************************************************************

C  SUBROUTINE SORTIA = THIS SUBROUTINE SORTS AN INTEGER ARRAY IN
C                      ASCENDING ORDER

C***********************************************************************

C VARIABLES   IN : N ...   NUMBER OF ELEMENTS IN THE ARRAY
C                  IAR ... INTEGER ARRAY WITH DATA TO BE SORTED
C             OUT: I ...   INDEX ARRAY WITH ITS SORTED VALUES IN
C                          ASCENDING ORDER.

C WRITTEN BY:  HORACIO RECALDE                 DATE: FEB 25,  1988
C***********************************************************************

      INTEGER I(N),IAR(N)

C-- COPY ELEMENTS IN THE I ARRAY

      DO 100 J = 1,IRANGE
         I(J) = IAR(J)
  100 CONTINUE

C---  PERFORM AN EXCHANGE SORT ON THE FIRST IRANGE-1

      DO 120 K = 1,IRANGE - 1
         MIN = I(K)

C---  EXCHANGE THE K-TH ELEMENTS WITH THE MINIMUM ELEMENT REMAIN

         DO 110 J = K+1,IRANGE
            IF (I(J) .LT. MIN) THEN
               L = I(J)
               I(J) = I(K)
               I(K) = L
               MIN = I(K)
            ENDIF
  110    CONTINUE
  120 CONTINUE

      RETURN
      END
