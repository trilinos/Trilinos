C Copyright (c) 2007 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
C retains certain rights in this software.
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
C     * Neither the name of Sandia Corporation nor the names of its
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
C 

C=======================================================================
*DECK,VOL
      SUBROUTINE VOL (ITYPE,XX,YY,ZZ,VOLUME)
C
C     ******************************************************************
C
C     SUBROUTINE TO FIND THE VOLUME OF AN ELEMENT
C
C     Called by ELTON3, ELTOND
C
C     ******************************************************************
C
      DIMENSION XX(*), YY(*), ZZ(*)
C
C     ******************************************************************
C
      IF (ITYPE .EQ. 3) THEN
C
C 4-node quad
C
        A1 = XX(1) - XX(3)
        B1 = YY(1) - YY(3)
        A2 = XX(2) - XX(4)
        B2 = YY(2) - YY(4)
        CCR = (A1 * B2) - (A2 * B1)
        VOLUME = 0.5 * ABS(CCR)
C
        RETURN
C
      ELSE IF (ITYPE .EQ. 10) THEN
C
C 8-node hex
C
        GR1 = ( YY(2)*(ZZ(6)-ZZ(3)-ZZ(4)+ZZ(5)) + YY(3)*(ZZ(2)-ZZ(4))
     &        + YY(4)*(ZZ(3)-ZZ(8)-ZZ(5)+ZZ(2)) 
     &        + YY(5)*(ZZ(8)-ZZ(6)-ZZ(2)+ZZ(4))
     &        + YY(6)*(ZZ(5)-ZZ(2)) + YY(8)*(ZZ(4)-ZZ(5)) ) / 12.
c
        GR2 = ( YY(3)*(ZZ(7)-ZZ(4)-ZZ(1)+ZZ(6)) + YY(4)*(ZZ(3)-ZZ(1))
     &        + YY(1)*(ZZ(4)-ZZ(5)-ZZ(6)+ZZ(3)) 
     &        + YY(6)*(ZZ(5)-ZZ(7)-ZZ(3)+ZZ(1))
     &        + YY(7)*(ZZ(6)-ZZ(3)) + YY(5)*(ZZ(1)-ZZ(6)) ) / 12.
c
        GR3 = ( YY(4)*(ZZ(8)-ZZ(1)-ZZ(2)+ZZ(7)) + YY(1)*(ZZ(4)-ZZ(2))
     &        + YY(2)*(ZZ(1)-ZZ(6)-ZZ(7)+ZZ(4)) 
     &        + YY(7)*(ZZ(6)-ZZ(8)-ZZ(4)+ZZ(2))
     &        + YY(8)*(ZZ(7)-ZZ(4)) + YY(6)*(ZZ(2)-ZZ(7)) ) / 12.
c
        GR4 = ( YY(1)*(ZZ(5)-ZZ(2)-ZZ(3)+ZZ(8)) + YY(2)*(ZZ(1)-ZZ(3))
     &        + YY(3)*(ZZ(2)-ZZ(7)-ZZ(8)+ZZ(1)) 
     &        + YY(8)*(ZZ(7)-ZZ(5)-ZZ(1)+ZZ(3))
     &        + YY(5)*(ZZ(8)-ZZ(1)) + YY(7)*(ZZ(3)-ZZ(8)) ) / 12.
c
        GR5 = ( YY(8)*(ZZ(4)-ZZ(7)-ZZ(6)+ZZ(1)) + YY(7)*(ZZ(8)-ZZ(6))
     &        + YY(6)*(ZZ(7)-ZZ(2)-ZZ(1)+ZZ(8)) 
     &        + YY(1)*(ZZ(2)-ZZ(4)-ZZ(8)+ZZ(6))
     &        + YY(4)*(ZZ(1)-ZZ(8)) + YY(2)*(ZZ(6)-ZZ(1)) ) / 12.
c
        GR6 = ( YY(5)*(ZZ(1)-ZZ(8)-ZZ(7)+ZZ(2)) + YY(8)*(ZZ(5)-ZZ(7))
     &        + YY(7)*(ZZ(8)-ZZ(3)-ZZ(2)+ZZ(5)) 
     &        + YY(2)*(ZZ(3)-ZZ(1)-ZZ(5)+ZZ(7))
     &        + YY(1)*(ZZ(2)-ZZ(5)) + YY(3)*(ZZ(7)-ZZ(2)) ) / 12.
c
        GR7 = ( YY(6)*(ZZ(2)-ZZ(5)-ZZ(8)+ZZ(3)) + YY(5)*(ZZ(6)-ZZ(8))
     &        + YY(8)*(ZZ(5)-ZZ(4)-ZZ(3)+ZZ(6)) 
     &        + YY(3)*(ZZ(4)-ZZ(2)-ZZ(6)+ZZ(8))
     &        + YY(2)*(ZZ(3)-ZZ(6)) + YY(4)*(ZZ(8)-ZZ(3)) ) / 12.
c
        GR8 = ( YY(7)*(ZZ(3)-ZZ(6)-ZZ(5)+ZZ(4)) + YY(6)*(ZZ(7)-ZZ(5))
     &        + YY(5)*(ZZ(6)-ZZ(1)-ZZ(4)+ZZ(7)) 
     &        + YY(4)*(ZZ(1)-ZZ(3)-ZZ(7)+ZZ(5))
     &        + YY(3)*(ZZ(4)-ZZ(7)) + YY(1)*(ZZ(5)-ZZ(4)) ) / 12.
c
        VOLUME = XX(1) * GR1 + XX(2) * GR2 + XX(3) * GR3 + XX(4) * GR4
     &         + XX(5) * GR5 + XX(6) * GR6 + XX(7) * GR7 + XX(8) * GR8
C
        RETURN
C
      ELSE IF (ITYPE .EQ. 13) THEN
C
C 4-node shell, ignore curvature
C
        A1 = XX(1) - XX(3)
        B1 = YY(1) - YY(3)
        C1 = ZZ(1) - ZZ(3)
        A2 = XX(2) - XX(4)
        B2 = YY(2) - YY(4)
        C2 = ZZ(2) - ZZ(4)
C
        ACR = (B1 * C2) - (B2 * C1)
        BCR = (A1 * C2) - (A2 * C1)
        CCR = (A1 * B2) - (A2 * B1)
C
        VOLUME = 0.5 * SQRT(ACR*ACR + BCR*BCR + CCR*CCR)
C
        RETURN
      ELSE IF (ITYPE .EQ. 6) THEN
C
C 4-node tet 
C (process Key's 8-node tet the same way using first 4 nodes)
C
        Y12 = YY(1) - YY(2)
        Y13 = YY(1) - YY(3)
        Y14 = YY(1) - YY(4)
        Y24 = YY(2) - YY(4)
        Y34 = YY(3) - YY(4)
        Z12 = ZZ(1) - ZZ(2)
        Z13 = ZZ(1) - ZZ(3)
        Z14 = ZZ(1) - ZZ(4)
        Z24 = ZZ(2) - ZZ(4)
        Z34 = ZZ(3) - ZZ(4)
C
        BX1 = 0.16666666667 * (Y34*Z24 - Y24*Z34)
        BX2 = 0.16666666667 * (Y13*Z14 - Y14*Z13)
        BX3 = 0.16666666667 * (Y14*Z12 - Y12*Z14)
        BX4 = 0.16666666667 * (Y12*Z13 - Y13*Z12)
C
        VOLUME = BX1*XX(1) + BX2*XX(2) + BX3*XX(3) + BX4*XX(4)
      ELSE
        CALL ERROR('VOL','ELEMENT TYPE',' ',ITYPE,
     &             'NOT YET IMPLEMENTED',0,' ',' ',1)
      END IF
      RETURN
      END
