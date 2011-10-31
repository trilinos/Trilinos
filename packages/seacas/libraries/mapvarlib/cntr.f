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
*DECK,CNTR
      SUBROUTINE CNTR (ITYPE,XX,YY,ZZ,XCNTR,YCNTR,ZCNTR)
C
C     ******************************************************************
C
C     SUBROUTINE TO FIND THE CENTER COORDINATES OF AN ELEMENT
C
C     Called by ELTON3, LDCENB, SRCH2D, SRCH3D
C
C     ******************************************************************
C
      DIMENSION XX(*), YY(*), ZZ(*)
C
C     ******************************************************************
C
C     SELECT ELEMENT TYPE
C
      GO TO (100, 110, 120, 130, 140, 150, 160, 170, 180, 190,
     &    200, 210, 150), ITYPE
C
C     3-NODE TRIANGLE
C
  100 CONTINUE
      XCNTR=(XX(1)+XX(2)+XX(3))/3.
      YCNTR=(YY(1)+YY(2)+YY(3))/3.
      ZCNTR=0.
      RETURN
C
C     6-NODE TRIANGLE
C
  110 CONTINUE
      XCNTR=-(XX(1)+XX(2)+XX(3))/9.+4.*(XX(4)+XX(5)+XX(6))/9.
      YCNTR=-(YY(1)+YY(2)+YY(3))/9.+4.*(YY(4)+YY(5)+YY(6))/9.
      ZCNTR=0.
      RETURN
C
C     4-NODE QUADRILATERAL
C
  120 CONTINUE
      XCNTR=(XX(1)+XX(2)+XX(3)+XX(4))/4.
      YCNTR=(YY(1)+YY(2)+YY(3)+YY(4))/4.
      ZCNTR=0.
      RETURN
C
C     8-NODE QUADRILATERAL
C
  130 CONTINUE
      XCNTR=-(XX(1)+XX(2)+XX(3)+XX(4))/4.+(XX(5)+XX(6)+XX(7)+XX(8))/2.
      YCNTR=-(YY(1)+YY(2)+YY(3)+YY(4))/4.+(YY(5)+YY(6)+YY(7)+YY(8))/2.
      ZCNTR=0.
      RETURN
C
C     9-NODE QUADRILATERAL
C
  140 CONTINUE
      XCNTR=XX(9)
      YCNTR=YY(9)
      ZCNTR=0.
      RETURN
C
C     4-NODE TETRAHEDRON OR SHELL
C
  150 CONTINUE
      XCNTR=(XX(1)+XX(2)+XX(3)+XX(4))/4.
      YCNTR=(YY(1)+YY(2)+YY(3)+YY(4))/4.
      ZCNTR=(ZZ(1)+ZZ(2)+ZZ(3)+ZZ(4))/4.
      RETURN
C
C     10-NODE TETRAHEDRON
C
  160 CONTINUE
      XCNTR=-(XX(1)+XX(2)+XX(3)+XX(4))/8.+(XX(5)+XX(6)+XX(7)+XX(8)+
     1XX(9)+XX(10))/4.
      YCNTR=-(YY(1)+YY(2)+YY(3)+YY(4))/8.+(YY(5)+YY(6)+YY(7)+YY(8)+
     1YY(9)+YY(10))/4.
      ZCNTR=-(ZZ(1)+ZZ(2)+ZZ(3)+ZZ(4))/8.+(ZZ(5)+ZZ(6)+ZZ(7)+ZZ(8)+
     1ZZ(9)+ZZ(10))/4.
      RETURN
C
C     6-NODE PRISM
C
  170 CONTINUE
      XCNTR=(XX(1)+XX(2)+XX(3)+XX(4)+XX(5)+XX(6))/6.
      YCNTR=(YY(1)+YY(2)+YY(3)+YY(4)+YY(5)+YY(6))/6.
      ZCNTR=(ZZ(1)+ZZ(2)+ZZ(3)+ZZ(4)+ZZ(5)+ZZ(6))/6.
      RETURN
C
C     15-NODE PRISM
C
  180 CONTINUE
      XCNTR=-2.*(XX(1)+XX(2)+XX(3)+XX(4)+XX(5)+XX(6))/9.+
     12.*(XX(7)+XX(8)+XX(9)+XX(13)+XX(14)+XX(15))/9.+
     2(XX(10)+XX(11)+XX(12))/3.
      YCNTR=-2.*(YY(1)+YY(2)+YY(3)+YY(4)+YY(5)+YY(6))/9.+
     12.*(YY(7)+YY(8)+YY(9)+YY(13)+YY(14)+YY(15))/9.+
     2(YY(10)+YY(11)+YY(12))/3.
      ZCNTR=-2.*(ZZ(1)+ZZ(2)+ZZ(3)+ZZ(4)+ZZ(5)+ZZ(6))/9.+
     12.*(ZZ(7)+ZZ(8)+ZZ(9)+ZZ(13)+ZZ(14)+ZZ(15))/9.+
     2(ZZ(10)+ZZ(11)+ZZ(12))/3.
      RETURN
C
C     8-NODE HEX
C
  190 CONTINUE
      XCNTR=(XX(1)+XX(2)+XX(3)+XX(4)+XX(5)+XX(6)+XX(7)+XX(8))/8.
      YCNTR=(YY(1)+YY(2)+YY(3)+YY(4)+YY(5)+YY(6)+YY(7)+YY(8))/8.
      ZCNTR=(ZZ(1)+ZZ(2)+ZZ(3)+ZZ(4)+ZZ(5)+ZZ(6)+ZZ(7)+ZZ(8))/8.
      RETURN
C
C     20-NODE HEX
C
  200 CONTINUE
      XCNTR=-(XX(1)+XX(2)+XX(3)+XX(4)+XX(5)+XX(6)+XX(7)+XX(8))/4.+
     1(XX(9)+XX(10)+XX(11)+XX(12)+XX(13)+XX(14)+XX(15)+XX(16)+XX(17)+
     2XX(18)+XX(19)+XX(20))/4.
      YCNTR=-(YY(1)+YY(2)+YY(3)+YY(4)+YY(5)+YY(6)+YY(7)+YY(8))/4.+
     1(YY(9)+YY(10)+YY(11)+YY(12)+YY(13)+YY(14)+YY(15)+YY(16)+YY(17)+
     2YY(18)+YY(19)+YY(20))/4.
      ZCNTR=-(ZZ(1)+ZZ(2)+ZZ(3)+ZZ(4)+ZZ(5)+ZZ(6)+ZZ(7)+ZZ(8))/4.+
     1(ZZ(9)+ZZ(10)+ZZ(11)+ZZ(12)+ZZ(13)+ZZ(14)+ZZ(15)+ZZ(16)+ZZ(17)+
     2ZZ(18)+ZZ(19)+ZZ(20))/4.
      RETURN
C
C     27-NODE HEX
C
  210 CONTINUE
      XCNTR=XX(27)
      YCNTR=YY(27)
      ZCNTR=ZZ(27)
      RETURN
C
      END
