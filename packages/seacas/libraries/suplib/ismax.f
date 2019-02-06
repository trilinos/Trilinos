C Copyright(C) 2009-2017 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.
C
C     * Neither the name of NTESS nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

C-------------------------------------------------------------     ************
C                                                                     ISMAX
C                                                                  ************
      INTEGER FUNCTION ISMAX(N,SX,INCX)
C
C        FINDS THE INDEX OF ELEMENT WITH MAX. VALUE.
C        JACK DONGARRA, LINPACK, 3/11/78.
C
      INTEGER I,INCX,IX,N
      REAL SX(1),SMAX
C
      ISMAX = 0
      IF (N .LT. 0) STOP
      IF (N .EQ. 0) RETURN
      ISMAX = 1
      IF (N .EQ. 1) RETURN
      IF (INCX .EQ. 1) GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
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
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   20 SMAX = SX(1)
      DO 30 I = 2,N,1
        IF (SX(I) .LE. SMAX) GO TO 30
        ISMAX = I
        SMAX = SX(I)
   30 CONTINUE
      RETURN
      END

