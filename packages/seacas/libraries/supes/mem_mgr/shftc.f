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
      SUBROUTINE SHFTC(IARAY,N,ISTRT,IEND,INT)
C
C     THIS SUBROUTINE SHIFTS THE ROWS IN A CHARACTER ARRAY. IF INT>0
C     THEN ALL ROWS ISTRT TO IEND ARE SHIFTED UP "INT" ROWS. IF INT<0
C     THEN ALL ROWS ISTRT TO IEND ARE SHIFTED DOWN "INT" ROWS.
C
C     MODIFIED JAN 1986 BY BILL MILLS-CURRAN
C     THE CHARACTER ARRAY IS NOW 1-D RATHER THAN 2-D.
C
      CHARACTER*(*) IARAY(N)
C
C  CALCULATE RANGE AND INCREMENT OF DO LOOP
C
      IF(INT.LT.0)THEN
C
C  SHIFT DOWN
C
         I1=IEND
         I2=ISTRT
         ID=-1
      ELSE IF(INT.GT.0)THEN
C
C  SHIFT UP
C
         I1=ISTRT
         I2=IEND
         ID=1
      ELSE
         RETURN
      END IF
C
C  PERFORM SHIFT
C
      DO 100 I=I1,I2,ID
  100 IARAY(I-INT)=IARAY(I)
      RETURN
      END
