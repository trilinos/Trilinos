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
      SUBROUTINE SRCHC(IVEC,ILBIN,IUBIN,NUM,ICOD,LOC)
C
C     THIS SUBROUTINE SEARCHES ROW ILB THOURGH ROW IUB OF A
C     NUMERICALLY ORDERED CHARACTER COLUMN VECTOR FOR THE OCCURANCE
C     THE VALUE NUM.
C     IF NUM IS FOUND ICOD IS SET TO UNITY AND LOC IS THE ROW NUMBER
C     WHERE NUM RESIDES. IF NUM IS NOT FOUND ICOD IS ZERO AND LOC IS
C     THE ROW NUMBER WHERE NUM WOULD RESIDE IF IT WERE IN THE NUMER-
C     ICALLY ORDER LIST.   BOB LUST
C
C     THIS SUBROUTINE HAS BEEN CHANGED FROM BOB LUST'S VERSION
C     AND NOW ASSUMES THAT THERE IS NO MORE THAN ONE MATCH
C     IN THE ORDERED LIST 'IVEC'.  BILL MILLS-CURRAN  JAN. 1, 1983
C
C     IVEC    ORDERED CHARACTER LIST (SINGLE COLUMN)
C
C     ILBIN   LOW NUMBERED ROW OF SEARCH RANGE
C
C     IUBIN   HIGH NUMBERED ROW OF SEARCH RANGE
C
C     NUM     VALUE TO BE LOCATED IN IVEC
C
C     ICOD    RETURN CODE  0 = NO MATCH   1 = MATCH
C
C     LOC     LOCATION IN IVEC FOR NUM
C
      CHARACTER*(*) IVEC(1),NUM
C
      ILB = ILBIN
      IUB = IUBIN
      ICOD=0
      IF (IUB .LT. ILB) THEN
         LOC = 1
         RETURN
      END IF
C
C     CHECK TO SEE IF NUM IS AT EITHER END OF LIST
C
      IF(IVEC(ILB).GT.NUM)THEN
         LOC=ILB
         RETURN
      ELSE IF(IVEC(IUB).LT.NUM) THEN
         LOC=IUB+1
         RETURN
      END IF
C
C     NUM IS INTERNAL TO IVEC
C
  100 MID=(ILB+IUB)/2
      IF(MID.LE.ILB)GO TO 110
C
C     SEARCH RANGE IS MORE THAN 2
C
      IF(IVEC(MID).LT.NUM) THEN
C
C        UPPER PART OF LIST
C
         ILB=MID
         GO TO 100
      ELSE IF(IVEC(MID).GT.NUM) THEN
C
C        LOWER PART OF LIST
C
         IUB=MID
         GO TO 100
      ELSE
C
C        MATCH HAS OCCURED AT "MID"
C
         ICOD=1
         LOC=MID
         RETURN
      END IF
  110 CONTINUE
C
C     SEARCH RANGE IS 2 OR LESS
C
      IF(NUM.EQ.IVEC(ILB)) THEN
C
C        MATCH AT "ILB"
C
         ICOD=1
         LOC=ILB
         RETURN
      ELSE IF(NUM.EQ.IVEC(IUB)) THEN
C
C        MATCH AT "IUB"
C
         ICOD=1
         LOC=IUB
         RETURN
      ELSE
C
C        NO MATCH IN LIST.
C        LOCATION FOR NEW ENTRY IS "IUB".
C
         LOC=IUB
         RETURN
      END IF
      END
