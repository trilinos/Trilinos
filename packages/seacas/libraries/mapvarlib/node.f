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
*DECK,NODE
      SUBROUTINE NODE (ITYPE,INODE,S,T,R)
C
C     ******************************************************************
C
C     SUBROUTINE TO SET THE NORMALIZED COORDINATES FOR A POINT
C     COINCIDENT WITH THE NODE OF AN ELEMENT
C
C     Called by SRCH2D & SRCH3D
C
C     ******************************************************************
C
      DIMENSION ST3(3), TT3(3), ST6(6), TT6(6)
      DIMENSION SQ4(4), TQ4(4), SQ8(8), TQ8(8), SQ9(9), TQ9(9)
      DIMENSION STT4(4), TTT4(4), RTT4(4)
      DIMENSION STT10(10), TTT10(10), RTT10(10)
      DIMENSION SP6(6), TP6(6), RP6(6)
      DIMENSION SP15(15), TP15(15), RP15(15)
      DIMENSION SB8(8), TB8(8), RB8(8)
      DIMENSION SB20(20), TB20(20), RB20(20)
      DIMENSION SB27(27), TB27(27), RB27(27)
C
      DATA (ST3(I),I=1,3)/1.,0.,0./
      DATA (TT3(I),I=1,3)/0.,1.,0./
C
      DATA (ST6(I),I=1,6)/1.,0.,0.,.5,0.,.5/
      DATA (TT6(I),I=1,6)/0.,1.,0.,.5,.5,0./
C
      DATA (SQ4(I),I=1,4)/-1.,1.,1.,-1./
      DATA (TQ4(I),I=1,4)/-1.,-1.,1.,1./
C
      DATA (SQ8(I),I=1,8)/-1.,1.,1.,-1.,0.,1.,0.,-1./
      DATA (TQ8(I),I=1,8)/-1.,-1.,1.,1.,-1.,0.,1.,0./
C
      DATA (SQ9(I),I=1,9)/-1.,1.,1.,-1.,0.,1.,0.,-1.,0./
      DATA (TQ9(I),I=1,9)/-1.,-1.,1.,1.,-1.,0.,1.,0.,0./
C
      DATA (STT4(I),I=1,4)/1.,0.,0.,0./
      DATA (TTT4(I),I=1,4)/0.,1.,0.,0./
      DATA (RTT4(I),I=1,4)/0.,0.,0.,1./
C
      DATA (STT10(I),I=1,10)/1.,0.,0.,0.,.25,0.,.25,.25,0.,0./
      DATA (TTT10(I),I=1,10)/0.,1.,0.,0.,.25,.25,0.,0.,.25,0./
      DATA (RTT10(I),I=1,10)/0.,0.,0.,1.,0.,0.,0.,.25,.25,.25/
C
      DATA (SP6(I),I=1,6)/1.,0.,0.,1.,0.,0./
      DATA (TP6(I),I=1,6)/0.,1.,0.,0.,1.,0./
      DATA (RP6(I),I=1,6)/-1.,-1.,-1.,1.,1.,1./
C
      DATA (SP15(I),I=1,15)/1.,0.,0.,1.,0.,0.,.5,0.,.5,1.,0.,           
     1     0.,.5,0.,.5/
      DATA (TP15(I),I=1,15)/0.,1.,0.,0.,1.,0.,.5,.5,0.,0.,1.,           
     1     0.,.5,.5,0./
      DATA (RP15(I),I=1,15)/-1.,-1.,-1.,1.,1.,1.,-1.,-1.,-1.,           
     1     0.,0.,0.,1.,1.,1./
C
      DATA (SB8(I),I=1,8)/-1.,1.,1.,-1.,-1.,1.,1.,-1./
      DATA (TB8(I),I=1,8)/-1.,-1.,1.,1.,-1.,-1.,1.,1./
      DATA (RB8(I),I=1,8)/-1.,-1.,-1.,-1.,1.,1.,1.,1./
C
      DATA (SB20(I),I=1,20)/-1.,1.,1.,-1.,-1.,1.,1.,-1.,0.,1.,0.,       
     1     -1.,-1.,1.,1.,-1.,0.,1.,0.,-1./
      DATA (TB20(I),I=1,20)/-1.,-1.,1.,1.,-1.,-1.,1.,1.,-1.,0.,1.,      
     1     0.,-1.,-1.,1.,1.,-1.,0.,1.,0./
      DATA (RB20(I),I=1,20)/-1.,-1.,-1.,-1.,1.,1.,1.,1.,-1.,-1.,        
     1     -1.,-1.,0.,0.,0.,0.,1.,1.,1.,1./
C
      DATA (SB27(I),I=1,27)/-1.,1.,1.,-1.,-1.,1.,1.,-1.,0.,1.,0.,       
     1     -1.,-1.,1.,1.,-1.,0.,1.,0.,-1.,0.,1.,0.,-1.,0.,0.,0./
      DATA (TB27(I),I=1,27)/-1.,-1.,1.,1.,-1.,-1.,1.,1.,-1.,0.,1.,      
     1     0.,-1.,-1.,1.,1.,-1.,0.,1.,0.,-1.,0.,1.,0.,0.,0.,0./
      DATA (RB27(I),I=1,27)/-1.,-1.,-1.,-1.,1.,1.,1.,1.,-1.,-1.,        
     1     -1.,-1.,0.,0.,0.,0.,1.,1.,1.,1.,0.,0.,0.,0.,-1.,1.,0./
C
C     ******************************************************************
C
C     SELECT ELEMENT TYPE
C
      GO TO (10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120), ITYPE
C
C     3-NODE TRIANGLE
C
   10 CONTINUE
      S=ST3(INODE)
      T=TT3(INODE)
      R=0.
      RETURN
C
C     6-NODE TRIANGLE
C
   20 CONTINUE
      S=ST6(INODE)
      T=TT6(INODE)
      R=0.
      RETURN
C
C     4-NODE QUADRILATERAL
C
   30 CONTINUE
      S=SQ4(INODE)
      T=TQ4(INODE)
      R=0.
      RETURN
C
C     8-NODE QUADRILATERAL
C
   40 CONTINUE
      S=SQ8(INODE)
      T=TQ8(INODE)
      R=0.
      RETURN
C
C     9-NODE QUADRILATERAL
C
   50 CONTINUE
      S=SQ9(INODE)
      T=TQ9(INODE)
      R=0.
      RETURN
C
C     4-NODE TETRAHEDRON
C
   60 CONTINUE
      S=STT4(INODE)
      T=TTT4(INODE)
      R=RTT4(INODE)
      RETURN
C
C     10-NODE TETRAHEDRON
C
   70 CONTINUE
      S=STT10(INODE)
      T=TTT10(INODE)
      R=RTT10(INODE)
      RETURN
C
C     6-NODE PRISM
C
   80 CONTINUE
      S=SP6(INODE)
      T=TP6(INODE)
      R=RP6(INODE)
      RETURN
C
C     15-NODE PRISM
C
   90 CONTINUE
      S=SP15(INODE)
      T=TP15(INODE)
      R=RP15(INODE)
      RETURN
C
C     8-NODE HEX
C
  100 CONTINUE
      S=SB8(INODE)
      T=TB8(INODE)
      R=RB8(INODE)
      RETURN
C
C     20-NODE HEX
C
  110 CONTINUE
      S=SB20(INODE)
      T=TB20(INODE)
      R=RB20(INODE)
      RETURN
C
C     27-NODE HEX
C
  120 CONTINUE
      S=SB27(INODE)
      T=TB27(INODE)
      R=RB27(INODE)
      RETURN
C
      END
