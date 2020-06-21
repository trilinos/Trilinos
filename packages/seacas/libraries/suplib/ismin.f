C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C-------------------------------------------------------------     ************
C                                                                     ISMIN
C                                                                  ************
      INTEGER FUNCTION ISMIN(N,SX,INCX)
C
C        FINDS THE INDEX OF ELEMENT WITH MIN. VALUE.
C        JACK DONGARRA, LINPACK, 3/11/78.
C
      INTEGER I,INCX,IX,N
      REAL SX(*),SMIN
C
      ISMIN = 0
      IF (N .LT. 0) STOP
      IF (N .EQ. 0) RETURN
      ISMIN = 1
      IF (N .EQ. 1) RETURN
      IF (INCX .EQ. 1) GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      IF (INCX .LE. 0) IX = (-N+1)*INCX + 1
      SMIN = SX(IX)
      IX = IX + INCX
      DO 10 I = 2,N,1
        IF (SX(IX) .GE. SMIN) GO TO 5
        ISMIN = I
        SMIN = SX(IX)
    5   IX = IX + INCX
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   20 SMIN = SX(1)
      DO 30 I = 2,N,1
        IF (SX(I) .GE. SMIN) GO TO 30
        ISMIN = I
        SMIN = SX(I)
   30 CONTINUE
      RETURN
      END

