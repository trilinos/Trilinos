C    Copyright(C) 2014-2017 National Technology & Engineering Solutions of
C    Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
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
C

CC* FILE: [.QMESH]CROSS.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE CCROSS (J1, J2, I1, I2, JXI, IXJ, ISTART, ICLEAR,
     &   NOROOM, ERR)
C***********************************************************************
C
C  SUBROUTINE CROSS = CREATE OR ADD TO THE CROSS - REFERENCE ARRAY FOR
C                     JXI (J1, J2) IN IXJ (I1, I2)
C
C***********************************************************************
C
C  NOTE:
C     THE NEW ITEMS MUST BEGIN AT J1=1,  J2=ISTART.
C     THE CROSS REFERENCE ARRAY WILL BE CLEARED FROM I1=1,  I2=ICLEAR
C     TO THE END OF THE ARRAY.
C
C***********************************************************************
C
      DIMENSION JXI (J1, J2), IXJ (I1, I2)
C
      LOGICAL ERR, NOROOM
C
C  CLEAR
C
      ERR = .TRUE.
      NOROOM = .FALSE.
      DO 110 J = ICLEAR, I2
         DO 100 I = 1, I1
            IXJ (I, J)  =  0
  100    CONTINUE
  110 CONTINUE
C
C  REFILE EACH ITEM
C
      DO 150 J = ISTART, J2
         DO 140 I = 1, J1
            L = IABS (JXI (I, J))
            IF (L .NE. 0) THEN
               IF (L .GT. I2) THEN
                  WRITE ( * , 10000)L, I2
                  RETURN
               ENDIF
C
C  FIND EMPTY SPOT FOR THIS ITEM
C
               DO 120 K = 1, I1
                  KK = K
                  IF (IXJ (K, L) .EQ. 0)GO TO 130
  120          CONTINUE
               CALL MESAGE ('NO ROOM FOR REFERENCE - ERROR IN CROSS')
               NOROOM = .TRUE.
               RETURN
  130          CONTINUE
C
C  FILE THIS ITEM
C
               IXJ (KK, L)  =  J
            ENDIF
  140    CONTINUE
  150 CONTINUE
      ERR = .FALSE.
      RETURN
C
10000 FORMAT (' OUT-OF-BOUNDS REFERENCE IN CROSS (INDEX = ', I5,
     &   ', MAX = ', I5, ')')
      END
