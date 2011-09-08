C    Copyright(C) 2008 Sandia Corporation.  Under the terms of Contract
C    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C    certain rights in this software
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
C    * Neither the name of Sandia Corporation nor the names of its
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
C    
      SUBROUTINE VTABLE (NEWLOC, NEWLEN, VOID, LVOID, NVOIDS, CHRCOL,
     *   ERR)
      IMPLICIT INTEGER (A-Z)
C
C     THIS SUBROUTINE INSERTS NEW VOIDS IN THE VOID TABLE AND
C     THEN CHECKS FOR CONTIGUOUS VOIDS WHICH ARE THEN JOINED.
C
C     ERROR CODES
C
C     ERROR VECTOR AND FLAGS.
C     THE ERROR PARAMETERS BELONG IN MDINIT ALSO.
C
      INCLUDE 'params.inc'
C
C     VFULL  = NO ROOM IN VOID TABLE
C     BDVOID = OVERLAPPING VOIDS
C
      DIMENSION VOID(LVOID,CHRCOL,2)
C
      IF (NEWLEN .GT. 0) THEN
C
         IF (NVOIDS .GE. LVOID) THEN
            ERR = VFULL
            RETURN
         END IF
C
C        FIND LOCATION FOR NEW ENTRY.
C
         CALL SRCHI(VOID,1,NVOIDS,NEWLOC,ERR,ROW)
         IF (ERR .NE. 0) THEN
            ERR = BDVOID
            RETURN
         END IF
C
C        NEW ENTRY IN TABLE.
C
         IF (ROW .LE. NVOIDS) THEN
C
C           MAKE ROOM FOR NEW ENTRY.
C
            CALL SHFTI (VOID, LVOID*CHRCOL, 2, ROW, NVOIDS, -1)
C
         END IF
C
         VOID(ROW,1,1) = NEWLOC
         VOID(ROW,1,2) = NEWLEN
         NVOIDS = NVOIDS + 1
C
      END IF
C
C     CHECK TABLE TO SEE IF ANY VOIDS HAVE JOINED OR ARE ZERO LENGTH.
C
C     NOTE THAT A STANDARD DO LOOP CANNOT BE USED BECAUSE THE UPPER
C     LIMIT OF THE LOOP CAN CHANGE INSIDE THE LOOP.
C
      I = 1
  100 IF (I .GE. NVOIDS) GO TO 110
         IF (VOID(I,1,1)+VOID(I,1,2) .EQ. VOID(I+1,1,1)) THEN
C
C           THESE TWO VOIDS SHOULD BE JOINED.
C
            VOID(I,1,2) = VOID(I,1,2) + VOID(I+1,1,2)
            CALL SHFTI (VOID, LVOID*CHRCOL, 2, I+2, NVOIDS, 1)
            NVOIDS = NVOIDS - 1
            GO TO 100
C
         ELSE IF (VOID(I,1,2) .EQ. 0) THEN
C
C           THIS VOID IS ZERO LENGTH.
C
            CALL SHFTI (VOID, LVOID*CHRCOL, 2, I+1, NVOIDS, 1)
            NVOIDS = NVOIDS - 1
C
         ELSE IF (VOID(I,1,1)+VOID(I,1,2) .GT. VOID(I+1,1,1)) THEN
C
C           OVERLAPPING VOIDS
C
            ERR = BDVOID
            RETURN
C
         END IF
C
         I = I + 1
         GO TO 100
C
  110 CONTINUE
C
C     CHECK LAST VOID
C
      IF (NVOIDS .GE. 1) THEN
         IF (VOID(NVOIDS,1,2) .EQ. 0) NVOIDS = NVOIDS - 1
      END IF
C
      ERR = SUCESS
      RETURN
      END
