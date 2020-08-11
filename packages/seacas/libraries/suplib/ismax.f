C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C-------------------------------------------------------------     ************
C                                                                     ISMAX
C                                                                  ************
      INTEGER FUNCTION ISMAX(N,SX,INCX)

C        FINDS THE INDEX OF ELEMENT WITH MAX. VALUE.
C        JACK DONGARRA, LINPACK, 3/11/78.

      INTEGER I,INCX,IX,N
      REAL SX(*),SMAX

      ISMAX = 0
      IF (N .LT. 0) STOP
      IF (N .EQ. 0) RETURN
      ISMAX = 1
      IF (N .EQ. 1) RETURN
      IF (INCX .EQ. 1) GO TO 20

C        CODE FOR INCREMENT NOT EQUAL TO 1

      IX = 1
      IF (INCX .LE. 0) IX = (-N+1)*INCX + 1
      SMAX = SX(IX)
      IX = IX + INCX
      DO 10 I = 2,N
        IF (SX(IX) .LE. SMAX) GO TO 5
        ISMAX = I
        SMAX = SX(IX)
    5   IX = IX + INCX
   10 CONTINUE
      RETURN

C        CODE FOR INCREMENT EQUAL TO 1

   20 SMAX = SX(1)
      DO 30 I = 2,N,1
        IF (SX(I) .LE. SMAX) GO TO 30
        ISMAX = I
        SMAX = SX(I)
   30 CONTINUE
      RETURN
      END

