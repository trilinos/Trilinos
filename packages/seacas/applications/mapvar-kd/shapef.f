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

C $Id: shapef.f,v 1.3 2007/10/17 18:43:50 gdsjaar Exp $ 
      SUBROUTINE SHAPEF (ITYPE,SP,TP,RP,SOLN,BVALUE)
C
C     ******************************************************************
C
C     SUBROUTINE TO EVALUATE THE SOLUTION AT A POINT (SP,TP,RP)
C     WITHIN AN ELEMENT
C
C     Called by INTRPE & INTRPN
C
C     ******************************************************************
C
C  SP     INT  The KSI isoparametric coord of the point in the element
C  TP     INT  The ETA isoparametric coord of the point in the element
C  RP     INT  The PHI isoparametric coord of the point in the element
C  SOLN   REAL The variables at the nodes of the element
C  BVALUE REAL The value of the variable at the point inside the element
C
C     ******************************************************************
C
      DIMENSION SOLN(*)
C
C     ******************************************************************
C
C     SELECT TYPE OF ELEMENT
C
      GO TO (10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120), ITYPE
C
C     3-NODE TRIANGLE
C
   10 CONTINUE
      PP=1.-SP-TP
      BVALUE=SOLN(1)*SP+SOLN(2)*TP+SOLN(3)*PP
      RETURN
C
C     6-NODE TRIANGLE
C
   20 CONTINUE
      PP=1.-SP-TP
      BVALUE=SOLN(1)*SP*(2.*SP-1.)+SOLN(2)*TP*(2.*TP-1.)+               
     1       SOLN(3)*PP*(2.*PP-1.)+SOLN(4)*4.*SP*TP+                    
     2       SOLN(5)*4.*TP*PP+SOLN(6)*4.*SP*PP
      RETURN
C
C     4-NODE QUADRILATERAL
C
   30 CONTINUE
      BVALUE=.25*(SOLN(1)*(1.-SP)*(1.-TP)+SOLN(2)*(1.+SP)*(1.-TP)+      
     1            SOLN(3)*(1.+SP)*(1.+TP)+SOLN(4)*(1.-SP)*(1.+TP))
      RETURN
C
C     8-NODE QUADRILATERAL
C
   40 CONTINUE
      BVALUE=.25*(SOLN(1)*(1.-SP)*(1.-TP)*(-SP-TP-1.)+                  
     1            SOLN(2)*(1.+SP)*(1.-TP)*(SP-TP-1.)+                   
     2            SOLN(3)*(1.+SP)*(1.+TP)*(SP+TP-1.)+                   
     3            SOLN(4)*(1.-SP)*(1.+TP)*(-SP+TP-1.))+                 
     4       .50*(SOLN(5)*(1.-SP*SP)*(1.-TP)+                           
     5            SOLN(6)*(1.+SP)*(1.-TP*TP)+                           
     6            SOLN(7)*(1.-SP*SP)*(1.+TP)+                           
     7            SOLN(8)*(1.-SP)*(1.-TP*TP))
      RETURN
C
C     9-NODE QUADRILATERAL
C
   50 CONTINUE
      BVALUE=.25*(SOLN(1)*(1.-SP)*(1.-TP)*(SP*TP)+                      
     1            SOLN(2)*(1.+SP)*(1.-TP)*(-SP*TP)+                     
     2            SOLN(3)*(1.+SP)*(1.+TP)*(SP*TP)+                      
     3            SOLN(4)*(1.-SP)*(1.+TP)*(-SP*TP))+                    
     4       .50*(SOLN(5)*(1.-SP*SP)*(1.-TP)*(-TP)+                     
     5            SOLN(6)*(1.+SP)*(1.-TP*TP)*SP+                        
     6            SOLN(7)*(1.-SP*SP)*(1.+TP)*TP+                        
     7            SOLN(8)*(1.-SP)*(1.-TP*TP)*(-SP))+                    
     8            SOLN(9)*(1.-SP*SP)*(1.-TP*TP)
      RETURN
C
C     4-NODE TETRAHEDRON
C
   60 CONTINUE
      PP=1.-SP-TP-RP
C      BVALUE=SOLN(1)*SP+SOLN(2)*TP+SOLN(3)*RP+SOLN(4)*PP
C fix gww 8/24/00 - probably same fix for 10-node tet but not done yet
C     merlin had opposite nodal order (left hand rule for + volume)
C
      BVALUE=SOLN(1)*SP+SOLN(2)*TP+SOLN(3)*PP+SOLN(4)*RP
      RETURN
C
C     10-NODE TETRAHEDRON
C
   70 CONTINUE
      PP=1.-SP-TP-RP
      BVALUE=  SOLN(1)*SP*(2.*SP-1.)+SOLN(2)*TP*(2.*TP-1.)+             
     1         SOLN(3)*RP*(2.*RP-1.)+SOLN(4)*PP*(2.*PP-1.)+             
     2     4.*(SOLN(5)*SP*TP+SOLN(6)*TP*RP+SOLN(7)*RP*SP+               
     3         SOLN(8)*SP*PP+SOLN(9)*TP*PP+SOLN(10)*RP*PP)
      RETURN
C
C     6-NODE PRISM
C
   80 CONTINUE
      PP=1.-SP-TP
      BVALUE=(SOLN(1)*SP+SOLN(2)*TP+SOLN(3)*PP)*.5*(1.-RP)+             
     1       (SOLN(4)*SP+SOLN(5)*TP+SOLN(6)*PP)*.5*(1.+RP)
      RETURN
C
C     15-NODE PRISM
C
   90 CONTINUE
      PP=1.-SP-TP
      AA=1.-RP
      BB=1.+RP
      CC=1.-RP**2
      BVALUE=.5*(SOLN(1)*SP*((2.*SP-1.)*AA-CC)+                         
     1           SOLN(2)*TP*((2.*TP-1.)*AA-CC)+                         
     2           SOLN(3)*PP*((2.*PP-1.)*AA-CC)+                         
     3           SOLN(4)*SP*((2.*SP-1.)*BB-CC)+                         
     4           SOLN(5)*TP*((2.*TP-1.)*BB-CC)+                         
     5           SOLN(6)*PP*((2.*PP-1.)*BB-CC))+                        
     6    2.*AA*(SOLN(7)*SP*TP+SOLN(8)*TP*PP+SOLN(9)*PP*SP)+            
     7          (SOLN(10)*SP+SOLN(11)*TP+SOLN(12)*PP)*CC+               
     8    2.*BB*(SOLN(13)*SP*TP+SOLN(14)*TP*PP+SOLN(15)*PP*SP)
      RETURN
C
C     8-NODE HEX
C
  100 CONTINUE
      BVALUE=.125*(SOLN(1)*(1.-SP)*(1.-TP)*(1.-RP)+                     
     1             SOLN(2)*(1.+SP)*(1.-TP)*(1.-RP)+                     
     2             SOLN(3)*(1.+SP)*(1.+TP)*(1.-RP)+                     
     3             SOLN(4)*(1.-SP)*(1.+TP)*(1.-RP)+                     
     4             SOLN(5)*(1.-SP)*(1.-TP)*(1.+RP)+                     
     5             SOLN(6)*(1.+SP)*(1.-TP)*(1.+RP)+                     
     6             SOLN(7)*(1.+SP)*(1.+TP)*(1.+RP)+                     
     7             SOLN(8)*(1.-SP)*(1.+TP)*(1.+RP))
      RETURN
C
C     20-NODE HEX
C
  110 CONTINUE
      AA=.125*(SOLN(1)*(1.-SP)*(1.-TP)*(1.-RP)*(-SP-TP-RP-2.)+          
     1         SOLN(2)*(1.+SP)*(1.-TP)*(1.-RP)*(SP-TP-RP-2.)+           
     2         SOLN(3)*(1.+SP)*(1.+TP)*(1.-RP)*(SP+TP-RP-2.)+           
     3         SOLN(4)*(1.-SP)*(1.+TP)*(1.-RP)*(-SP+TP-RP-2.)+          
     4         SOLN(5)*(1.-SP)*(1.-TP)*(1.+RP)*(-SP-TP+RP-2.)+          
     5         SOLN(6)*(1.+SP)*(1.-TP)*(1.+RP)*(SP-TP+RP-2.)+           
     6         SOLN(7)*(1.+SP)*(1.+TP)*(1.+RP)*(SP+TP+RP-2.)+           
     7         SOLN(8)*(1.-SP)*(1.+TP)*(1.+RP)*(-SP+TP+RP-2.))
      BB=.250*(SOLN(9)*(1.-SP**2)*(1.-TP)*(1.-RP)+                      
     1         SOLN(10)*(1.+SP)*(1.-TP**2)*(1.-RP)+                     
     2         SOLN(11)*(1.-SP**2)*(1.+TP)*(1.-RP)+                     
     3         SOLN(12)*(1.-SP)*(1.-TP**2)*(1.-RP)+                     
     4         SOLN(13)*(1.-SP)*(1.-TP)*(1.-RP**2)+                     
     5         SOLN(14)*(1.+SP)*(1.-TP)*(1.-RP**2)+                     
     6         SOLN(15)*(1.+SP)*(1.+TP)*(1.-RP**2)+                     
     7         SOLN(16)*(1.-SP)*(1.+TP)*(1.-RP**2)+                     
     8         SOLN(17)*(1.-SP**2)*(1.-TP)*(1.+RP)+                     
     9         SOLN(18)*(1.+SP)*(1.-TP**2)*(1.+RP)+                     
     $         SOLN(19)*(1.-SP**2)*(1.+TP)*(1.+RP)+                     
     $         SOLN(20)*(1.-SP)*(1.-TP**2)*(1.+RP))
      BVALUE=AA+BB
      RETURN
C
C     27-NODE HEX
C
  120 CONTINUE
      AA=.125*(-SOLN(1)*SP*TP*RP*(1.-SP)*(1.-TP)*(1.-RP)                
     1         +SOLN(2)*SP*TP*RP*(1.+SP)*(1.-TP)*(1.-RP)                
     2         -SOLN(3)*SP*TP*RP*(1.+SP)*(1.+TP)*(1.-RP)                
     3         +SOLN(4)*SP*TP*RP*(1.-SP)*(1.+TP)*(1.-RP)                
     4         +SOLN(5)*SP*TP*RP*(1.-SP)*(1.-TP)*(1.+RP)                
     5         -SOLN(6)*SP*TP*RP*(1.+SP)*(1.-TP)*(1.+RP)                
     6         +SOLN(7)*SP*TP*RP*(1.+SP)*(1.+TP)*(1.+RP)                
     7         -SOLN(8)*SP*TP*RP*(1.-SP)*(1.+TP)*(1.+RP))
      BB=.250*(+SOLN (9)*TP*RP*(1.-SP**2)*(1.-TP)*(1.-RP)               
     1         -SOLN(10)*SP*RP*(1.+SP)*(1.-TP**2)*(1.-RP)               
     2         -SOLN(11)*TP*RP*(1.-SP**2)*(1.+TP)*(1.-RP)               
     3         +SOLN(12)*SP*RP*(1.-SP)*(1.-TP**2)*(1.-RP)               
     4         +SOLN(13)*SP*TP*(1.-SP)*(1.-TP)*(1.-RP**2)               
     5         -SOLN(14)*SP*TP*(1.+SP)*(1.-TP)*(1.-RP**2)               
     6         +SOLN(15)*SP*TP*(1.+SP)*(1.+TP)*(1.-RP**2)               
     7         -SOLN(16)*SP*TP*(1.-SP)*(1.+TP)*(1.-RP**2))
      CC=.250*(-SOLN(17)*TP*RP*(1.-SP**2)*(1.-TP)*(1.+RP)               
     1         +SOLN(18)*SP*RP*(1.+SP)*(1.-TP**2)*(1.+RP)               
     2         +SOLN(19)*TP*RP*(1.-SP**2)*(1.+TP)*(1.+RP)               
     3         -SOLN(20)*SP*RP*(1.-SP)*(1.-TP**2)*(1.+RP))+.500*(-SOLN(2
     41)*TP*(1.-SP**2)*(1.-TP)*(1.-RP**2)                        +SOLN(2
     52)*SP*(1.+SP)*(1.-TP**2)*(1.-RP**2)                        +SOLN(2
     63)*TP*(1.-SP**2)*(1.+TP)*(1.-RP**2)                        -SOLN(2
     74)*SP*(1.-SP)*(1.-TP**2)*(1.-RP**2)                        -SOLN(2
     85)*RP*(1.-SP**2)*(1.-TP**2)*(1.-RP)                        +SOLN(2
     96)*RP*(1.-SP**2)*(1.-TP**2)*(1.+RP))+SOLN(27)*(1.-SP**2)*(1.-TP**2
     $)*(1.-RP**2)
      BVALUE=AA+BB+CC
      RETURN
C
      END
