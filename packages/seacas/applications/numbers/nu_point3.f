C    Copyright(C) 1988-2017 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C    * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C              
C    * Redistributions in binary form must reproduce the above
C      copyright notice, this list of conditions and the following
C      disclaimer in the documentation and/or other materials provided
C      with the distribution.
C                            
C    * Neither the name of NTESS nor the names of its
C      contributors may be used to endorse or promote products derived
C      from this software without specific prior written permission.
C                                                    
C    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

      SUBROUTINE POINT3 (COORD, NUMNP, DIST, NDIM, P1, TOLER,
     *   NODEL, SORTYP, MAP, ANGLE, SORUP, INUM, OPT, SELECT)
      DIMENSION COORD (NUMNP,*), DIST(*), P1(*), TOLER(*),
     *   MAP(*), ANGLE(2,*)
      CHARACTER*(*) NODEL, SORTYP, OPT
      LOGICAL SORUP, SELECT(*), ISABRT
      include 'nu_io.blk'
      PI = ATAN2(0.0, -1.0)
C
      CALL LOCOUT ('POINT', NDIM, NODEL, TOLER, SORTYP, P1, P1, ' ')

      TEMP = TOLER(1)
      TOLER(1) = MAX(0.0, TEMP - TOLER(2))
      TOLER(2) = MAX(0.0, TEMP + TOLER(2))
C
      X1 = P1(1)
      Y1 = P1(2)
      Z1 = P1(3)
C
      DO 10 I=1, NUMNP
         IF (SELECT(I)) THEN
         X0 = COORD(I,1)
         Y0 = COORD(I,2)
         Z0 = COORD(I,3)
C
         DIST(I) = (X1 - X0)**2 + (Y1 - Y0)**2 + (Z1 - Z0)**2
         END IF
   10 CONTINUE

      INUM = 0
      DISMIN = 1.0E38
      DO 20 I=1, NUMNP
         IF (SELECT(I)) THEN
         DISMIN = MIN(DIST(I), ABS(DISMIN-TEMP))
         IF (DIST(I) .GE. TOLER(1)**2 .AND. DIST(I) .LE. TOLER(2)**2)
     *      THEN
            DX = COORD(I,1) - P1(1)
            DY = COORD(I,2) - P1(2)
            DZ = COORD(I,3) - P1(3)
            DS = SQRT(DX**2 + DZ**2)
C
            FIX1 = SIGN(0.5,ABS(DX+DZ)) + SIGN(0.5,-ABS(DX+DZ))
            FIX2 = SIGN(0.5,ABS(DS+DY)) + SIGN(0.5,-ABS(DS+DY))
            ANGLE(1,I) = ATAN2(DZ,DX+FIX1) * 180.0 / PI
            ANGLE(2,I) = ATAN2(DY,DS+FIX2) * 180.0 / PI
C
            INUM = INUM + 1
            MAP(INUM) = I
         END IF
         END IF
   20 CONTINUE
C
      IF      (SORTYP .EQ. 'X') THEN
         CALL INDEXX (COORD(1,1), MAP, INUM, .FALSE.)
      ELSE IF (SORTYP .EQ. 'Y') THEN
         CALL INDEXX (COORD(1,2), MAP, INUM, .FALSE.)
      ELSE IF (SORTYP .EQ. 'Z') THEN
         CALL INDEXX (COORD(1,3), MAP, INUM, .FALSE.)
      ELSE IF (SORTYP .EQ. 'THETA') THEN
         CALL INDEXN (ANGLE,      2, 1, MAP, INUM, .FALSE.)
      ELSE IF (SORTYP .EQ. 'PHI') THEN
         CALL INDEXN (ANGLE,      2, 2, MAP, INUM, .FALSE.)
      ELSE IF (SORTYP .EQ. 'DISTANCE') THEN
         CALL INDEXX (DIST,       MAP, INUM, .FALSE.)
      END IF

      IF (SORUP) THEN
         IBEG = 1
         IEND = INUM
         IINC = 1
      ELSE
         IBEG = INUM
         IEND = 1
         IINC = -1
      END IF

      IF (OPT .EQ. '*' .OR. INDEX(OPT, 'P') .GT. 0) THEN
      DO 30 IO=IOMIN, IOMAX
         WRITE (IO, 40) NODEL
   30 CONTINUE
   40 FORMAT (/,2X,A8,T16,'X',T26,'Y',T36,'Z',T44,'DISTANCE',
     *   T58,'THETA',T71,'PHI')
      DO 60 IN = IBEG, IEND, IINC
         IF (ISABRT()) RETURN
         I = MAP(IN)
         DO 50 IO=IOMIN, IOMAX
            WRITE (IO, 90) I, (COORD(I,J),J=1,3), SQRT(DIST(I)),
     *         ANGLE(1,I), ANGLE(2,I)
   50    CONTINUE
   60 CONTINUE

      IF (INUM .EQ. 0) THEN
         DO 70 IO=IOMIN, IOMAX
            WRITE (IO, 80) SQRT(DISMIN)
   70    CONTINUE
      END IF
      END IF
   80 FORMAT (/' None found within range, minimum distance = ',
     *   1PE12.3,/)
   90 FORMAT (I10, 3(F10.4), 3(1PE12.3))
      RETURN
      END
