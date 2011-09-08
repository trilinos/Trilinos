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

C $Id: jacobn.f,v 1.3 2007/10/17 18:40:35 gdsjaar Exp $
C $Log: jacobn.f,v $
C Revision 1.3  2007/10/17 18:40:35  gdsjaar
C Added copyright notice to all files.
C
C Mapvar is licensed under the BSD license
C
C Revision 1.2  2003/09/11 12:47:52  gwwellm
C tet element processing and 3 years of bug fixes I thought I had already committed
C
C Revision 1.1  1998/03/13 18:12:23  gdsjaar
C New code -- mapvar. Interpolates results form an exodusII results file
C to a differently mesh geometry.  Written by Gerry Wellman,
C 9117. Loosely based on MERLIN. Provides a superset of merlin
C functionality.
C
c Revision 1.3  1993/02/12  20:47:36  gdsjaar
c Formatting mods to make it easier to figure out the code
c
c Revision 1.2  1991/12/02  19:56:16  gdsjaar
c New version from Gartling
c
C=======================================================================
*DECK,JACOBN
      SUBROUTINE JACOBN (ITYPE,XX,YY,ZZ,SP,TP,RP,A11,A12,A13,A21,A22,
     1A23,A31,A32,A33,F1,F2,F3)
C
C     ******************************************************************
C
C     SUBROUTINE TO EVALUATE ELEMENT SHAPE FUNCTIONS AND DERIVATIVES
C     (JACOBIAN) AT A SPECIFIED POINT (SP,TP,RP)
C
C     Called by SRCH2D & SRCH3D
C
C     ******************************************************************
C
      DIMENSION XX(*), YY(*), ZZ(*)
C
C     ******************************************************************
C
      F1=0.
      F2=0.
      F3=0.
      A11=0.
      A12=0.
      A13=0.
      A21=0.
      A22=0.
      A23=0.
      A31=0.
      A32=0.
      A33=0.
C
C     SELECT ELEMENT
C
      GO TO (100, 110, 120, 130, 140, 150, 160, 170, 180, 190,
     &    200, 210), ITYPE
C
C     3-NODE TRIANGLE
C
  100 CONTINUE
      PP=1.-SP-TP
      F1=XX(1)*SP+XX(2)*TP+XX(3)*PP
      F2=YY(1)*SP+YY(2)*TP+YY(3)*PP
C
      A11=XX(1)-XX(3)
      A12=XX(2)-XX(3)
      A21=YY(1)-YY(3)
      A22=YY(2)-YY(3)
      RETURN
C
C     6-NODE TRIANGLE
C
  110 CONTINUE
      PP=1.-SP-TP
      F1=XX(1)*SP*(2.*SP-1.)+XX(2)*TP*(2.*TP-1.)+
     1   XX(3)*PP*(2.*PP-1.)+XX(4)*4.*SP*TP+
     2   XX(5)*4.*TP*PP+XX(6)*4.*SP*PP
      F2=YY(1)*SP*(2.*SP-1.)+YY(2)*TP*(2.*TP-1.)+
     1   YY(3)*PP*(2.*PP-1.)+YY(4)*4.*SP*TP+
     2   YY(5)*4.*TP*PP+YY(6)*4.*SP*PP
C
      A11=XX(1)*(4.*SP-1.)-XX(3)*(4.*PP-1.)+(XX(4)-XX(5))*4.*TP+
     1    XX(6)*4.*(PP-SP)
      A12=XX(2)*(4.*TP-1.)-XX(3)*(4.*PP-1.)+(XX(4)-XX(6))*4.*SP+
     1    XX(5)*4.*(PP-TP)
      A21=YY(1)*(4.*SP-1.)-YY(3)*(4.*PP-1.)+(YY(4)-YY(5))*4.*TP+
     1    YY(6)*4.*(PP-SP)
      A22=YY(2)*(4.*TP-1.)-YY(3)*(4.*PP-1.)+(YY(4)-YY(6))*4.*SP+
     1    YY(5)*4.*(PP-TP)
      RETURN
C
C     4-NODE QUADRILATERAL
C
  120 CONTINUE
      F1=.25*(XX(1)*(1.-SP)*(1.-TP)+XX(2)*(1.+SP)*(1.-TP)+
     1        XX(3)*(1.+SP)*(1.+TP)+XX(4)*(1.-SP)*(1.+TP))
      F2=.25*(YY(1)*(1.-SP)*(1.-TP)+YY(2)*(1.+SP)*(1.-TP)+
     1        YY(3)*(1.+SP)*(1.+TP)+YY(4)*(1.-SP)*(1.+TP))
C
      A11=.25*((XX(2)-XX(1))*(1.-TP)+(XX(3)-XX(4))*(1.+TP))
      A12=.25*((XX(4)-XX(1))*(1.-SP)+(XX(3)-XX(2))*(1.+SP))
      A21=.25*((YY(2)-YY(1))*(1.-TP)+(YY(3)-YY(4))*(1.+TP))
      A22=.25*((YY(4)-YY(1))*(1.-SP)+(YY(3)-YY(2))*(1.+SP))
      RETURN
C
C     8-NODE QUADRILATERAL
C
  130 CONTINUE
      F1=.25*(XX(1)*(1.-SP)*(1.-TP)*(-SP-TP-1.)+
     1        XX(2)*(1.+SP)*(1.-TP)*(SP-TP-1.)+
     2        XX(3)*(1.+SP)*(1.+TP)*(SP+TP-1.)+
     3        XX(4)*(1.-SP)*(1.+TP)*(-SP+TP-1.))+
     4   .50*(XX(5)*(1.-SP*SP)*(1.-TP)+XX(6)*(1.+SP)*(1.-TP*TP)+
     5        XX(7)*(1.-SP*SP)*(1.+TP)+XX(8)*(1.-SP)*(1.-TP*TP))
      F2=.25*(YY(1)*(1.-SP)*(1.-TP)*(-SP-TP-1.)+
     1        YY(2)*(1.+SP)*(1.-TP)*(SP-TP-1.)+
     2        YY(3)*(1.+SP)*(1.+TP)*(SP+TP-1.)+
     3        YY(4)*(1.-SP)*(1.+TP)*(-SP+TP-1.))+
     4   .50*(YY(5)*(1.-SP*SP)*(1.-TP)+YY(6)*(1.+SP)*(1.-TP*TP)+
     5        YY(7)*(1.-SP*SP)*(1.+TP)+YY(8)*(1.-SP)*(1.-TP*TP))
C
      A11=.25*(XX(1)*(1.-TP)*(TP+2.*SP)+XX(2)*(1.-TP)*(-TP+2.*SP)+
     1         XX(3)*(1.+TP)*(TP+2.*SP)+XX(4)*(1.+TP)*(-TP+2.*SP))+
     2    .50*(XX(5)*(1.-TP)*(-2.*SP)+(XX(6)-XX(8))*(1.-TP*TP)+
     3         XX(7)*(1.+TP)*(-2.*SP))
      A12=.25*(XX(1)*(1.-SP)*(SP+2.*TP)+XX(2)*(1.+SP)*(-SP+2.*TP)+
     1         XX(3)*(1.+SP)*(SP+2.*TP)+XX(4)*(1.-SP)*(-SP+2.*TP))+
     2    .50*(XX(6)*(1.+SP)*(-2.*TP)+(XX(7)-XX(5))*(1.-SP*SP)+
     3         XX(8)*(1.-SP)*(-2.*TP))
      A21=.25*(YY(1)*(1.-TP)*(TP+2.*SP)+YY(2)*(1.-TP)*(-TP+2.*SP)+
     1         YY(3)*(1.+TP)*(TP+2.*SP)+YY(4)*(1.+TP)*(-TP+2.*SP))+
     2    .50*(YY(5)*(1.-TP)*(-2.*SP)+(YY(6)-YY(8))*(1.-TP*TP)+
     3         YY(7)*(1.+TP)*(-2.*SP))
      A22=.25*(YY(1)*(1.-SP)*(SP+2.*TP)+YY(2)*(1.+SP)*(-SP+2.*TP)+
     1         YY(3)*(1.+SP)*(SP+2.*TP)+YY(4)*(1.-SP)*(-SP+2.*TP))+
     2    .50*(YY(6)*(1.+SP)*(-2.*TP)+(YY(7)-YY(5))*(1.-SP*SP)+
     3         YY(8)*(1.-SP)*(-2.*TP))
      RETURN
C
C     9-NODE QUADRILATERAL
C
  140 CONTINUE
      F1=.25*(XX(1)*(1.-SP)*(1.-TP)*(SP*TP)+
     1        XX(2)*(1.+SP)*(1.-TP)*(-SP*TP)+
     2        XX(3)*(1.+SP)*(1.+TP)*(SP*TP)+
     3        XX(4)*(1.-SP)*(1.+TP)*(-SP*TP))+
     4   .50*(XX(5)*(1.-SP*SP)*(1.-TP)*(-TP)+
     5        XX(6)*(1.+SP)*(1.-TP*TP)*SP+
     6        XX(7)*(1.-SP*SP)*(1.+TP)*TP+
     7        XX(8)*(1.-SP)*(1.-TP*TP)*(-SP))+
     8        XX(9)*(1.-SP*SP)*(1.-TP*TP)
      F2=.25*(YY(1)*(1.-SP)*(1.-TP)*(SP*TP)+
     1        YY(2)*(1.+SP)*(1.-TP)*(-SP*TP)+
     2        YY(3)*(1.+SP)*(1.+TP)*(SP*TP)+
     3        YY(4)*(1.-SP)*(1.+TP)*(-SP*TP))+
     4   .50*(YY(5)*(1.-SP*SP)*(1.-TP)*(-TP)+
     5        YY(6)*(1.+SP)*(1.-TP*TP)*SP+
     6        YY(7)*(1.-SP*SP)*(1.+TP)*TP+
     7        YY(8)*(1.-SP)*(1.-TP*TP)*(-SP))+
     8        YY(9)*(1.-SP*SP)*(1.-TP*TP)
C
      A11=.25*(XX(1)*(1.-TP)*(TP-2.*SP*TP)+
     1         XX(2)*(1.-TP)*(-TP-2.*SP*TP)+
     2         XX(3)*(1.+TP)*(TP+2.*SP*TP)+
     3         XX(4)*(1.+TP)*(-TP+2.*SP*TP))+
     4    .50*(XX(5)*(-TP+TP*TP)*(-2.*SP)+XX(6)*(1.-TP*TP)*(1.+2*SP)+
     5         XX(7)*(TP+TP*TP)*(-2.*SP)+XX(8)*(1.-TP*TP)*(-1.+2.*SP))+
     6         XX(9)*(1.-TP*TP)*(-2.*SP)
      A12=.25*(XX(1)*(1.-SP)*(SP-2.*SP*TP)+
     1         XX(2)*(1.+SP)*(-SP+2.*SP*TP)+
     2         XX(3)*(1.+SP)*(SP+2.*SP*TP)+
     3         XX(4)*(1.-SP)*(-SP-2.*SP*TP))+
     4    .50*(XX(5)*(1.-SP*SP)*(-1.+2.*TP)+XX(6)*(SP+SP*SP)*(-2.*TP)+
     5         XX(7)*(1.-SP*SP)*(1.+2.*TP)+XX(8)*(-SP+SP*SP)*(-2.*TP))+
     6         XX(9)*(1.-SP*SP)*(-2.*TP)
      A21=.25*(YY(1)*(1.-TP)*(TP-2.*SP*TP)+
     1         YY(2)*(1.-TP)*(-TP-2.*SP*TP)+
     2         YY(3)*(1.+TP)*(TP+2.*SP*TP)+
     3         YY(4)*(1.+TP)*(-TP+2.*SP*TP))+
     4    .50*(YY(5)*(-TP+TP*TP)*(-2.*SP)+YY(6)*(1.-TP*TP)*(1.+2*SP)+
     5         YY(7)*(TP+TP*TP)*(-2.*SP)+YY(8)*(1.-TP*TP)*(-1.+2.*SP))+
     6         YY(9)*(1.-TP*TP)*(-2.*SP)
      A22=.25*(YY(1)*(1.-SP)*(SP-2.*SP*TP)+
     1         YY(2)*(1.+SP)*(-SP+2.*SP*TP)+
     2         YY(3)*(1.+SP)*(SP+2.*SP*TP)+
     3         YY(4)*(1.-SP)*(-SP-2.*SP*TP))+
     4    .50*(YY(5)*(1.-SP*SP)*(-1.+2.*TP)+YY(6)*(SP+SP*SP)*(-2.*TP)+
     5         YY(7)*(1.-SP*SP)*(1.+2.*TP)+YY(8)*(-SP+SP*SP)*(-2.*TP))+
     6         YY(9)*(1.-SP*SP)*(-2.*TP)
      RETURN
C
C     4-NODE TETRAHEDRON
C
  150 CONTINUE
      PP=1.-SP-TP-RP
C      F1=XX(1)*SP+XX(2)*TP+XX(3)*RP+XX(4)*PP
C      F2=YY(1)*SP+YY(2)*TP+YY(3)*RP+YY(4)*PP
C      F3=ZZ(1)*SP+ZZ(2)*TP+ZZ(3)*RP+ZZ(4)*PP
C
C      A11=XX(1)-XX(4)
C      A12=XX(2)-XX(4)
C      A13=XX(3)-XX(4)
C      A21=YY(1)-YY(4)
C      A22=YY(2)-YY(4)
C      A23=YY(3)-YY(4)
C      A31=ZZ(1)-ZZ(4)
C      A32=ZZ(2)-ZZ(4)
C      A33=ZZ(3)-ZZ(4)
C fix gww 8/24/00 - need same fix for 10-node tet?? 
C     merlin had opposite nodal order (left hand rule for + volume)
C
      F1=XX(1)*SP+XX(2)*TP+XX(3)*PP+XX(4)*RP
      F2=YY(1)*SP+YY(2)*TP+YY(3)*PP+YY(4)*RP
      F3=ZZ(1)*SP+ZZ(2)*TP+ZZ(3)*PP+ZZ(4)*RP
C
      A11=XX(1)-XX(3)
      A12=XX(2)-XX(3)
      A13=XX(4)-XX(3)
      A21=YY(1)-YY(3)
      A22=YY(2)-YY(3)
      A23=YY(4)-YY(3)
      A31=ZZ(1)-ZZ(3)
      A32=ZZ(2)-ZZ(3)
      A33=ZZ(4)-ZZ(3)
      RETURN
C
C     10-NODE TETRAHEDRON
C
  160 CONTINUE
      PP=1.-SP-TP-RP
      F1=  XX(1)*SP*(2.*SP-1.)+XX(2)*TP*(2.*TP-1.)+
     1     XX(3)*RP*(2.*RP-1.)+XX(4)*PP*(2.*PP-1.)+
     2 4.*(XX(5)*SP*TP+XX(6)*TP*RP+XX(7)*RP*SP+
     3     XX(8)*SP*PP+XX(9)*TP*PP+XX(10)*RP*PP)
      F2=  YY(1)*SP*(2.*SP-1.)+YY(2)*TP*(2.*TP-1.)+
     1     YY(3)*RP*(2.*RP-1.)+YY(4)*PP*(2.*PP-1.)+
     2 4.*(YY(5)*SP*TP+YY(6)*TP*RP+YY(7)*RP*SP+
     3     YY(8)*SP*PP+YY(9)*TP*PP+YY(10)*RP*PP)
      F3=  ZZ(1)*SP*(2.*SP-1.)+ZZ(2)*TP*(2.*TP-1.)+
     1     ZZ(3)*RP*(2.*RP-1.)+ZZ(4)*PP*(2.*PP-1.)+
     2 4.*(ZZ(5)*SP*TP+ZZ(6)*TP*RP+ZZ(7)*RP*SP+
     3     ZZ(8)*SP*PP+ZZ(9)*TP*PP+ZZ(10)*RP*PP)
C
      A11=XX(1)*(4.*SP-1.)-XX(4)*(4.*PP-1.)+XX(5)*4.*TP+XX(7)*4.*RP+
     1    XX(8)*4.*(PP-SP)-XX(9)*4.*TP-XX(10)*4.*(RP+TP)
      A12=XX(2)*(4.*TP-1.)-XX(4)*(4.*PP-1.)+XX(5)*4.*SP+XX(6)*4.*RP-
     1    XX(8)*4.*SP+XX(9)*4.*(PP-TP)-XX(10)*4.*(RP+SP)
      A13=XX(3)*(4.*RP-1.)-XX(4)*(4.*PP-1.)+XX(6)*4.*TP+XX(7)*4.*SP-
     1    XX(8)*4.*SP-XX(9)*4.*TP+XX(10)*(4.*PP-RP)
      A21=YY(1)*(4.*SP-1.)-YY(4)*(4.*PP-1.)+YY(5)*4.*TP+YY(7)*4.*RP+
     1    YY(8)*4.*(PP-SP)-YY(9)*4.*TP-YY(10)*4.*(RP+TP)
      A22=YY(2)*(4.*TP-1.)-YY(4)*(4.*PP-1.)+YY(5)*4.*SP+YY(6)*4.*RP-
     1    YY(8)*4.*SP+YY(9)*4.*(PP-TP)-YY(10)*4.*(RP+SP)
      A23=YY(3)*(4.*RP-1.)-YY(4)*(4.*PP-1.)+YY(6)*4.*TP+YY(7)*4.*SP-
     1    YY(8)*4.*SP-YY(9)*4.*TP+YY(10)*(4.*PP-RP)
      A31=ZZ(1)*(4.*SP-1.)-ZZ(4)*(4.*PP-1.)+ZZ(5)*4.*TP+ZZ(7)*4.*RP+
     1    ZZ(8)*4.*(PP-SP)-ZZ(9)*4.*TP-ZZ(10)*4.*(RP+TP)
      A32=ZZ(2)*(4.*TP-1.)-ZZ(4)*(4.*PP-1.)+ZZ(5)*4.*SP+ZZ(6)*4.*RP-
     1    ZZ(8)*4.*SP+ZZ(9)*4.*(PP-TP)-ZZ(10)*4.*(RP+SP)
      A33=ZZ(3)*(4.*RP-1.)-ZZ(4)*(4.*PP-1.)+ZZ(6)*4.*TP+ZZ(7)*4.*SP-
     1    ZZ(8)*4.*SP-ZZ(9)*4.*TP+ZZ(10)*(4.*PP-RP)
      RETURN
C
C     6-NODE PRISM
C
  170 CONTINUE
      PP=1.-SP-TP
      F1=(XX(1)*SP+XX(2)*TP+XX(3)*PP)*.5*(1.-RP)+
     1   (XX(4)*SP+XX(5)*TP+XX(6)*PP)*.5*(1.+RP)
      F2=(YY(1)*SP+YY(2)*TP+YY(3)*PP)*.5*(1.-RP)+
     1   (YY(4)*SP+YY(5)*TP+YY(6)*PP)*.5*(1.+RP)
      F3=(ZZ(1)*SP+ZZ(2)*TP+ZZ(3)*PP)*.5*(1.-RP)+
     1   (ZZ(4)*SP+ZZ(5)*TP+ZZ(6)*PP)*.5*(1.+RP)
C
      A11=(XX(1)-XX(3))*.5*(1.-RP)+
     1    (XX(4)-XX(6))*.5*(1.+RP)
      A12=(XX(2)-XX(3))*.5*(1.-RP)+
     1    (XX(5)-XX(6))*.5*(1.+RP)
      A13=(XX(4)-XX(1))*.5*SP+
     1    (XX(5)-XX(2))*.5*TP+
     2    (XX(6)-XX(3))*.5*PP
      A21=(YY(1)-YY(3))*.5*(1.-RP)+
     1    (YY(4)-YY(6))*.5*(1.+RP)
      A22=(YY(2)-YY(3))*.5*(1.-RP)+
     1    (YY(5)-YY(6))*.5*(1.+RP)
      A23=(YY(4)-YY(1))*.5*SP+
     1    (YY(5)-YY(2))*.5*TP+
     2    (YY(6)-YY(3))*.5*PP
      A31=(ZZ(1)-ZZ(3))*.5*(1.-RP)+
     1    (ZZ(4)-ZZ(6))*.5*(1.+RP)
      A32=(ZZ(2)-ZZ(3))*.5*(1.-RP)+
     1    (ZZ(5)-ZZ(6))*.5*(1.+RP)
      A33=(ZZ(4)-ZZ(1))*.5*SP+
     1    (ZZ(5)-ZZ(2))*.5*TP+
     2    (ZZ(6)-ZZ(3))*.5*PP
      RETURN
C
C     15-NODE PRISM
C
  180 CONTINUE
      PP=1.-SP-TP
      AA=1.-RP
      BB=1.+RP
      CC=1.-RP**2
      F1=.50*(XX(1)*SP*((2.*SP-1.)*AA-CC)+XX(2)*TP*((2.*TP-1.)*AA-CC)+
     1        XX(3)*PP*((2.*PP-1.)*AA-CC)+
     2        XX(4)*SP*((2.*SP-1.)*BB-CC)+XX(5)*TP*((2.*TP-1.)*BB-CC)+
     3        XX(6)*PP*((2.*PP-1.)*BB-CC))+
     4 2.*AA*(XX(7)*SP*TP+XX(8)*TP*PP+XX(9)*PP*SP)+
     5       (XX(10)*SP+XX(11)*TP+XX(12)*PP)*CC+
     6 2.*BB*(XX(13)*SP*TP+XX(14)*TP*PP+XX(15)*PP*SP)
      F2=.50*(YY(1)*SP*((2.*SP-1.)*AA-CC)+YY(2)*TP*((2.*TP-1.)*AA-CC)+
     1        YY(3)*PP*((2.*PP-1.)*AA-CC)+
     2        YY(4)*SP*((2.*SP-1.)*BB-CC)+YY(5)*TP*((2.*TP-1.)*BB-CC)+
     3        YY(6)*PP*((2.*PP-1.)*BB-CC))+
     4 2.*AA*(YY(7)*SP*TP+YY(8)*TP*PP+YY(9)*PP*SP)+
     5       (YY(10)*SP+YY(11)*TP+YY(12)*PP)*CC+
     6 2.*BB*(YY(13)*SP*TP+YY(14)*TP*PP+YY(15)*PP*SP)
      F3=.50*(ZZ(1)*SP*((2.*SP-1.)*AA-CC)+ZZ(2)*TP*((2.*TP-1.)*AA-CC)+
     1        ZZ(3)*PP*((2.*PP-1.)*AA-CC)+
     2        ZZ(4)*SP*((2.*SP-1.)*BB-CC)+ZZ(5)*TP*((2.*TP-1.)*BB-CC)+
     3        ZZ(6)*PP*((2.*PP-1.)*BB-CC))+
     4 2.*AA*(ZZ(7)*SP*TP+ZZ(8)*TP*PP+ZZ(9)*PP*SP)+
     5       (ZZ(10)*SP+ZZ(11)*TP+ZZ(12)*PP)*CC+
     6 2.*BB*(ZZ(13)*SP*TP+ZZ(14)*TP*PP+ZZ(15)*PP*SP)
C
      A11=.5*(XX(1)*((4.*SP-1.)*AA-CC)-XX(3)*((4.*PP-1.)*AA+CC)+
     1        XX(4)*((4.*SP-1.)*BB-CC)-XX(6)*((4.*PP-1.)*BB+CC))+
     2 2.*AA*(XX(7)*TP-XX(8)*TP+XX(9)*(PP-SP))+
     3       (XX(10)-XX(12))*CC+
     4 2.*BB*(XX(13)*TP-XX(14)*TP+XX(15)*(PP-SP))
      A12=.5*(XX(2)*((4.*TP-1.)*AA-CC)-XX(3)*((4.*PP-1.)*AA+CC)+
     1        XX(5)*((4.*TP-1.)*BB-CC)-XX(6)*((4.*PP-1.)*BB+CC))+
     2 2.*AA*(XX(7)*TP+XX(8)*(PP-TP)-XX(9)*SP)+
     3       (XX(11)-XX(12))*CC+
     4 2.*BB*(XX(13)*SP+XX(14)*(PP-TP)-XX(15)*SP)
      A13=   (XX(4)-XX(1))*.5*SP*(2.*SP-1.+RP)+
     1       (XX(5)-XX(2))*.5*TP*(2.*TP-1.+RP)+
     2       (XX(6)-XX(3))*.5*PP*(2.*PP-1.+RP)-
     3    2.*(XX(7)*SP*TP+XX(8)*TP*PP+XX(9)*PP*SP)-
     4    2.*(XX(10)*SP+XX(11)*TP+XX(12)*PP)*RP+
     5    2.*(XX(13)*SP*TP+XX(14)*TP*PP+XX(15)*PP*SP)
      A21=.5*(YY(1)*((4.*SP-1.)*AA-CC)-YY(3)*((4.*PP-1.)*AA+CC)+
     1        YY(4)*((4.*SP-1.)*BB-CC)-YY(6)*((4.*PP-1.)*BB+CC))+
     2 2.*AA*(YY(7)*TP-YY(8)*TP+YY(9)*(PP-SP))+
     3       (YY(10)-YY(12))*CC+
     4 2.*BB*(YY(13)*TP-YY(14)*TP+YY(15)*(PP-SP))
      A22=.5*(YY(2)*((4.*TP-1.)*AA-CC)-YY(3)*((4.*PP-1.)*AA+CC)+
     1        YY(5)*((4.*TP-1.)*BB-CC)-YY(6)*((4.*PP-1.)*BB+CC))+
     2 2.*AA*(YY(7)*TP+YY(8)*(PP-TP)-YY(9)*SP)+
     3       (YY(11)-YY(12))*CC+
     4 2.*BB*(YY(13)*SP+YY(14)*(PP-TP)-YY(15)*SP)
      A23=   (YY(4)-YY(1))*.5*SP*(2.*SP-1.+RP)+
     1       (YY(5)-YY(2))*.5*TP*(2.*TP-1.+RP)+
     2       (YY(6)-YY(3))*.5*PP*(2.*PP-1.+RP)-
     3    2.*(YY(7)*SP*TP+YY(8)*TP*PP+YY(9)*PP*SP)-
     4    2.*(YY(10)*SP+YY(11)*TP+YY(12)*PP)*RP+
     5    2.*(YY(13)*SP*TP+YY(14)*TP*PP+YY(15)*PP*SP)
      A31=.5*(ZZ(1)*((4.*SP-1.)*AA-CC)-ZZ(3)*((4.*PP-1.)*AA+CC)+
     1        ZZ(4)*((4.*SP-1.)*BB-CC)-ZZ(6)*((4.*PP-1.)*BB+CC))+
     2 2.*AA*(ZZ(7)*TP-ZZ(8)*TP+ZZ(9)*(PP-SP))+
     3       (ZZ(10)-ZZ(12))*CC+
     4 2.*BB*(ZZ(13)*TP-ZZ(14)*TP+ZZ(15)*(PP-SP))
      A32=.5*(ZZ(2)*((4.*TP-1.)*AA-CC)-ZZ(3)*((4.*PP-1.)*AA+CC)+
     1        ZZ(5)*((4.*TP-1.)*BB-CC)-ZZ(6)*((4.*PP-1.)*BB+CC))+
     2 2.*AA*(ZZ(7)*TP+ZZ(8)*(PP-TP)-ZZ(9)*SP)+
     3       (ZZ(11)-ZZ(12))*CC+
     4 2.*BB*(ZZ(13)*SP+ZZ(14)*(PP-TP)-ZZ(15)*SP)
      A33=   (ZZ(4)-ZZ(1))*.5*SP*(2.*SP-1.+RP)+
     1       (ZZ(5)-ZZ(2))*.5*TP*(2.*TP-1.+RP)+
     2       (ZZ(6)-ZZ(3))*.5*PP*(2.*PP-1.+RP)-
     3    2.*(ZZ(7)*SP*TP+ZZ(8)*TP*PP+ZZ(9)*PP*SP)-
     4    2.*(ZZ(10)*SP+ZZ(11)*TP+ZZ(12)*PP)*RP+
     5    2.*(ZZ(13)*SP*TP+ZZ(14)*TP*PP+ZZ(15)*PP*SP)
      RETURN
C
C     8-NODE HEX
C
  190 CONTINUE
      F1=.125*(XX(1)*(1.-SP)*(1.-TP)*(1.-RP)+
     1         XX(2)*(1.+SP)*(1.-TP)*(1.-RP)+
     2         XX(3)*(1.+SP)*(1.+TP)*(1.-RP)+
     3         XX(4)*(1.-SP)*(1.+TP)*(1.-RP)+
     4         XX(5)*(1.-SP)*(1.-TP)*(1.+RP)+
     5         XX(6)*(1.+SP)*(1.-TP)*(1.+RP)+
     6         XX(7)*(1.+SP)*(1.+TP)*(1.+RP)+
     7         XX(8)*(1.-SP)*(1.+TP)*(1.+RP))
      F2=.125*(YY(1)*(1.-SP)*(1.-TP)*(1.-RP)+
     1         YY(2)*(1.+SP)*(1.-TP)*(1.-RP)+
     2         YY(3)*(1.+SP)*(1.+TP)*(1.-RP)+
     3         YY(4)*(1.-SP)*(1.+TP)*(1.-RP)+
     4         YY(5)*(1.-SP)*(1.-TP)*(1.+RP)+
     5         YY(6)*(1.+SP)*(1.-TP)*(1.+RP)+
     6         YY(7)*(1.+SP)*(1.+TP)*(1.+RP)+
     7         YY(8)*(1.-SP)*(1.+TP)*(1.+RP))
      F3=.125*(ZZ(1)*(1.-SP)*(1.-TP)*(1.-RP)+
     1         ZZ(2)*(1.+SP)*(1.-TP)*(1.-RP)+
     2         ZZ(3)*(1.+SP)*(1.+TP)*(1.-RP)+
     3         ZZ(4)*(1.-SP)*(1.+TP)*(1.-RP)+
     4         ZZ(5)*(1.-SP)*(1.-TP)*(1.+RP)+
     5         ZZ(6)*(1.+SP)*(1.-TP)*(1.+RP)+
     6         ZZ(7)*(1.+SP)*(1.+TP)*(1.+RP)+
     7         ZZ(8)*(1.-SP)*(1.+TP)*(1.+RP))
C
      A11=.125*((XX(2)-XX(1))*(1.-TP)*(1.-RP)+
     1          (XX(3)-XX(4))*(1.+TP)*(1.-RP)+
     2          (XX(6)-XX(5))*(1.-TP)*(1.+RP)+
     3          (XX(7)-XX(8))*(1.+TP)*(1.+RP))
      A12=.125*((XX(4)-XX(1))*(1.-SP)*(1.-RP)+
     1          (XX(3)-XX(2))*(1.+SP)*(1.-RP)+
     2          (XX(8)-XX(5))*(1.-SP)*(1.+RP)+
     3          (XX(7)-XX(6))*(1.+SP)*(1.+RP))
      A13=.125*((XX(5)-XX(1))*(1.-SP)*(1.-TP)+
     1          (XX(6)-XX(2))*(1.+SP)*(1.-TP)+
     2          (XX(7)-XX(3))*(1.+SP)*(1.+TP)+
     3          (XX(8)-XX(4))*(1.-SP)*(1.+TP))
      A21=.125*((YY(2)-YY(1))*(1.-TP)*(1.-RP)+
     1          (YY(3)-YY(4))*(1.+TP)*(1.-RP)+
     2          (YY(6)-YY(5))*(1.-TP)*(1.+RP)+
     3          (YY(7)-YY(8))*(1.+TP)*(1.+RP))
      A22=.125*((YY(4)-YY(1))*(1.-SP)*(1.-RP)+
     1          (YY(3)-YY(2))*(1.+SP)*(1.-RP)+
     2          (YY(8)-YY(5))*(1.-SP)*(1.+RP)+
     3          (YY(7)-YY(6))*(1.+SP)*(1.+RP))
      A23=.125*((YY(5)-YY(1))*(1.-SP)*(1.-TP)+
     1          (YY(6)-YY(2))*(1.+SP)*(1.-TP)+
     2          (YY(7)-YY(3))*(1.+SP)*(1.+TP)+
     3          (YY(8)-YY(4))*(1.-SP)*(1.+TP))
      A31=.125*((ZZ(2)-ZZ(1))*(1.-TP)*(1.-RP)+
     1          (ZZ(3)-ZZ(4))*(1.+TP)*(1.-RP)+
     2          (ZZ(6)-ZZ(5))*(1.-TP)*(1.+RP)+
     3          (ZZ(7)-ZZ(8))*(1.+TP)*(1.+RP))
      A32=.125*((ZZ(4)-ZZ(1))*(1.-SP)*(1.-RP)+
     1          (ZZ(3)-ZZ(2))*(1.+SP)*(1.-RP)+
     2          (ZZ(8)-ZZ(5))*(1.-SP)*(1.+RP)+
     3          (ZZ(7)-ZZ(6))*(1.+SP)*(1.+RP))
      A33=.125*((ZZ(5)-ZZ(1))*(1.-SP)*(1.-TP)+
     1          (ZZ(6)-ZZ(2))*(1.+SP)*(1.-TP)+
     2          (ZZ(7)-ZZ(3))*(1.+SP)*(1.+TP)+
     3          (ZZ(8)-ZZ(4))*(1.-SP)*(1.+TP))
      RETURN
C
C     20-NODE HEX
C
  200 CONTINUE
      AA=.125*(XX(1)*(1.-SP)*(1.-TP)*(1.-RP)*(-SP-TP-RP-2.)+
     1         XX(2)*(1.+SP)*(1.-TP)*(1.-RP)*(SP-TP-RP-2.)+
     2         XX(3)*(1.+SP)*(1.+TP)*(1.-RP)*(SP+TP-RP-2.)+
     3         XX(4)*(1.-SP)*(1.+TP)*(1.-RP)*(-SP+TP-RP-2.)+
     4         XX(5)*(1.-SP)*(1.-TP)*(1.+RP)*(-SP-TP+RP-2.)+
     5         XX(6)*(1.+SP)*(1.-TP)*(1.+RP)*(SP-TP+RP-2.)+
     6         XX(7)*(1.+SP)*(1.+TP)*(1.+RP)*(SP+TP+RP-2.)+
     7         XX(8)*(1.-SP)*(1.+TP)*(1.+RP)*(-SP+TP+RP-2.))
      BB=.250*(XX(9)*(1.-SP**2)*(1.-TP)*(1.-RP)+
     1         XX(10)*(1.+SP)*(1.-TP**2)*(1.-RP)+
     2         XX(11)*(1.-SP**2)*(1.+TP)*(1.-RP)+
     3         XX(12)*(1.-SP)*(1.-TP**2)*(1.-RP)+
     4         XX(13)*(1.-SP)*(1.-TP)*(1.-RP**2)+
     5         XX(14)*(1.+SP)*(1.-TP)*(1.-RP**2)+
     6         XX(15)*(1.+SP)*(1.+TP)*(1.-RP**2)+
     7         XX(16)*(1.-SP)*(1.+TP)*(1.-RP**2)+
     8         XX(17)*(1.-SP**2)*(1.-TP)*(1.+RP)+
     9         XX(18)*(1.+SP)*(1.-TP**2)*(1.+RP)+
     #         XX(19)*(1.-SP**2)*(1.+TP)*(1.+RP)+
     1         XX(20)*(1.-SP)*(1.-TP**2)*(1.+RP))
      F1=AA+BB
      AA=.125*(YY(1)*(1.-SP)*(1.-TP)*(1.-RP)*(-SP-TP-RP-2.)+
     1         YY(2)*(1.+SP)*(1.-TP)*(1.-RP)*(SP-TP-RP-2.)+
     2         YY(3)*(1.+SP)*(1.+TP)*(1.-RP)*(SP+TP-RP-2.)+
     3         YY(4)*(1.-SP)*(1.+TP)*(1.-RP)*(-SP+TP-RP-2.)+
     4         YY(5)*(1.-SP)*(1.-TP)*(1.+RP)*(-SP-TP+RP-2.)+
     5         YY(6)*(1.+SP)*(1.-TP)*(1.+RP)*(SP-TP+RP-2.)+
     6         YY(7)*(1.+SP)*(1.+TP)*(1.+RP)*(SP+TP+RP-2.)+
     7         YY(8)*(1.-SP)*(1.+TP)*(1.+RP)*(-SP+TP+RP-2.))
      BB=.250*(YY(9)*(1.-SP**2)*(1.-TP)*(1.-RP)+
     1         YY(10)*(1.+SP)*(1.-TP**2)*(1.-RP)+
     2         YY(11)*(1.-SP**2)*(1.+TP)*(1.-RP)+
     3         YY(12)*(1.-SP)*(1.-TP**2)*(1.-RP)+
     4         YY(13)*(1.-SP)*(1.-TP)*(1.-RP**2)+
     5         YY(14)*(1.+SP)*(1.-TP)*(1.-RP**2)+
     6         YY(15)*(1.+SP)*(1.+TP)*(1.-RP**2)+
     7         YY(16)*(1.-SP)*(1.+TP)*(1.-RP**2)+
     8         YY(17)*(1.-SP**2)*(1.-TP)*(1.+RP)+
     9         YY(18)*(1.+SP)*(1.-TP**2)*(1.+RP)+
     #         YY(19)*(1.-SP**2)*(1.+TP)*(1.+RP)+
     1         YY(20)*(1.-SP)*(1.-TP**2)*(1.+RP))
      F2=AA+BB
      AA=.125*(ZZ(1)*(1.-SP)*(1.-TP)*(1.-RP)*(-SP-TP-RP-2.)+
     1         ZZ(2)*(1.+SP)*(1.-TP)*(1.-RP)*(SP-TP-RP-2.)+
     2         ZZ(3)*(1.+SP)*(1.+TP)*(1.-RP)*(SP+TP-RP-2.)+
     3         ZZ(4)*(1.-SP)*(1.+TP)*(1.-RP)*(-SP+TP-RP-2.)+
     4         ZZ(5)*(1.-SP)*(1.-TP)*(1.+RP)*(-SP-TP+RP-2.)+
     5         ZZ(6)*(1.+SP)*(1.-TP)*(1.+RP)*(SP-TP+RP-2.)+
     6         ZZ(7)*(1.+SP)*(1.+TP)*(1.+RP)*(SP+TP+RP-2.)+
     7         ZZ(8)*(1.-SP)*(1.+TP)*(1.+RP)*(-SP+TP+RP-2.))
      BB=.250*(ZZ(9)*(1.-SP**2)*(1.-TP)*(1.-RP)+
     1         ZZ(10)*(1.+SP)*(1.-TP**2)*(1.-RP)+
     2         ZZ(11)*(1.-SP**2)*(1.+TP)*(1.-RP)+
     3         ZZ(12)*(1.-SP)*(1.-TP**2)*(1.-RP)+
     4         ZZ(13)*(1.-SP)*(1.-TP)*(1.-RP**2)+
     5         ZZ(14)*(1.+SP)*(1.-TP)*(1.-RP**2)+
     6         ZZ(15)*(1.+SP)*(1.+TP)*(1.-RP**2)+
     7         ZZ(16)*(1.-SP)*(1.+TP)*(1.-RP**2)+
     8         ZZ(17)*(1.-SP**2)*(1.-TP)*(1.+RP)+
     9         ZZ(18)*(1.+SP)*(1.-TP**2)*(1.+RP)+
     #         ZZ(19)*(1.-SP**2)*(1.+TP)*(1.+RP)+
     1         ZZ(20)*(1.-SP)*(1.-TP**2)*(1.+RP))
      F3=AA+BB
C
      AA=.125*(-XX(1)*(1.-TP)*(1.-RP)*(-2.*SP-TP-RP-1.)
     1         +XX(2)*(1.-TP)*(1.-RP)*(2.*SP-TP-RP-1.)
     2         +XX(3)*(1.+TP)*(1.-RP)*(-2.*SP+TP-RP-1.)
     3         -XX(4)*(1.+TP)*(1.-RP)*(-2.*SP+TP-RP-1.)
     4         -XX(5)*(1.-TP)*(1.+RP)*(-2.*SP-TP+RP-1.)
     5         +XX(6)*(1.-TP)*(1.+RP)*(2.*SP-TP+RP-1.)
     6         +XX(7)*(1.+TP)*(1.+RP)*(2.*SP+TP+RP-1.)
     7         -XX(8)*(1.+TP)*(1.+RP)*(-2.*SP+TP+RP-1.))
      BB=.250*(+XX(9)*(-2.*SP)*(1.-TP)*(1.-RP)
     1         +XX(10)*(1.-TP**2)*(1.-RP)
     2         +XX(11)*(-2.*SP)*(1.+TP)*(1.-RP)
     3         -XX(12)*(1.-TP**2)*(1.-RP)
     4         -XX(13)*(1.-TP)*(1.-RP**2)
     5         +XX(14)*(1.-TP)*(1.-RP**2)
     6         +XX(15)*(1.+TP)*(1.-RP**2)
     7         -XX(16)*(1.+TP)*(1.-RP**2)
     8         +XX(17)*(-2.*SP)*(1.-TP)*(1.+RP)
     9         +XX(18)*(1.-TP**2)*(1.+RP)
     #         +XX(19)*(-2.*SP)*(1.+TP)*(1.+RP)
     1         -XX(20)*(1.-TP**2)*(1.+RP))
      A11=AA+BB
      AA=.125*(-XX(1)*(1.-SP)*(1.-RP)*(-SP-2.*TP-RP-1.)
     1         -XX(2)*(1.+SP)*(1.-RP)*(SP-2.*TP-RP-1.)
     2         +XX(3)*(1.+SP)*(1.-RP)*(SP+2.*TP-RP-1.)
     3         +XX(4)*(1.-SP)*(1.-RP)*(-SP+2.*TP-RP-1.)
     4         -XX(5)*(1.-SP)*(1.+RP)*(-SP-2.*TP+RP-1.)
     5         -XX(6)*(1.+SP)*(1.+RP)*(SP-2.*TP+RP-1.)
     6         +XX(7)*(1.+SP)*(1.+RP)*(SP+2.*TP+RP-1.)
     7         +XX(8)*(1.-SP)*(1.+RP)*(-SP+2.*TP+RP-1.))
      BB=.250*(-XX(9)*(1.-SP**2)*(1.-RP)
     1         +XX(10)*(1.+SP)*(-2.*TP)*(1.-RP)
     2         +XX(11)*(1.-SP**2)*(1.-RP)
     3         +XX(12)*(1.-SP)*(-2.*TP)*(1.-RP)
     4         -XX(13)*(1.-SP)*(1.-RP**2)
     5         -XX(14)*(1.+SP)*(1.-RP**2)
     6         +XX(15)*(1.+SP)*(1.-RP**2)
     7         +XX(16)*(1.-SP)*(1.-RP**2)
     8         -XX(17)*(1.-SP**2)*(1.+RP)
     9         +XX(18)*(1.+SP)*(-2.*TP)*(1.+RP)
     #         +XX(19)*(1.-SP**2)*(1.+RP)
     1         +XX(20)*(1.-SP)*(-2.*TP)*(1.+RP))
      A12=AA+BB
      AA=.125*(-XX(1)*(1.-SP)*(1.-TP)*(-SP-TP-2.*RP-1.)
     1         -XX(2)*(1.+SP)*(1.-TP)*(SP-TP-2.*RP-1.)
     2         -XX(3)*(1.+SP)*(1.+TP)*(SP+TP-2.*RP-1.)
     3         -XX(4)*(1.-SP)*(1.+TP)*(-SP+TP-2.*RP-1.)
     4         +XX(5)*(1.-SP)*(1.-TP)*(-SP-TP+2.*RP-1.)
     5         +XX(6)*(1.+SP)*(1.-TP)*(SP-TP+2.*RP-1.)
     6         +XX(7)*(1.+SP)*(1.+TP)*(SP+TP+2.*RP-1.)
     7         +XX(8)*(1.-SP)*(1.+TP)*(-SP+TP+2.*RP-1.))
      BB=.250*(-XX(9)*(1.-SP**2)*(1.-TP)
     1         -XX(10)*(1.+SP)*(1.-TP**2)
     2         -XX(11)*(1.-SP**2)*(1.+TP)
     3         -XX(12)*(1.-SP)*(1.-TP**2)
     4         +XX(13)*(1.-SP)*(1.-TP)*(-2.*RP)
     5         +XX(14)*(1.+SP)*(1.-TP)*(-2.*RP)
     6         +XX(15)*(1.+SP)*(1.+TP)*(-2.*RP)
     7         +XX(16)*(1.-SP)*(1.+TP)*(-2.*RP)
     8         +XX(17)*(1.-SP**2)*(1.-TP)
     9         +XX(18)*(1.+SP)*(1.-TP**2)
     #         +XX(19)*(1.-SP**2)*(1.+TP)
     1         +XX(20)*(1.-SP)*(1.-TP**2))
      A13=AA+BB
C
      AA=.125*(-YY(1)*(1.-TP)*(1.-RP)*(-2.*SP-TP-RP-1.)
     1         +YY(2)*(1.-TP)*(1.-RP)*(2.*SP-TP-RP-1.)
     2         +YY(3)*(1.+TP)*(1.-RP)*(-2.*SP+TP-RP-1.)
     3         -YY(4)*(1.+TP)*(1.-RP)*(-2.*SP+TP-RP-1.)
     4         -YY(5)*(1.-TP)*(1.+RP)*(-2.*SP-TP+RP-1.)
     5         +YY(6)*(1.-TP)*(1.+RP)*(2.*SP-TP+RP-1.)
     6         +YY(7)*(1.+TP)*(1.+RP)*(2.*SP+TP+RP-1.)
     7         -YY(8)*(1.+TP)*(1.+RP)*(-2.*SP+TP+RP-1.))
      BB=.250*(+YY(9)*(-2.*SP)*(1.-TP)*(1.-RP)
     1         +YY(10)*(1.-TP**2)*(1.-RP)
     2         +YY(11)*(-2.*SP)*(1.+TP)*(1.-RP)
     3         -YY(12)*(1.-TP**2)*(1.-RP)
     4         -YY(13)*(1.-TP)*(1.-RP**2)
     5         +YY(14)*(1.-TP)*(1.-RP**2)
     6         +YY(15)*(1.+TP)*(1.-RP**2)
     7         -YY(16)*(1.+TP)*(1.-RP**2)
     8         +YY(17)*(-2.*SP)*(1.-TP)*(1.+RP)
     9         +YY(18)*(1.-TP**2)*(1.+RP)
     #         +YY(19)*(-2.*SP)*(1.+TP)*(1.+RP)
     1         -YY(20)*(1.-TP**2)*(1.+RP))
      A21=AA+BB
      AA=.125*(-YY(1)*(1.-SP)*(1.-RP)*(-SP-2.*TP-RP-1.)
     1         -YY(2)*(1.+SP)*(1.-RP)*(SP-2.*TP-RP-1.)
     2         +YY(3)*(1.+SP)*(1.-RP)*(SP+2.*TP-RP-1.)
     3         +YY(4)*(1.-SP)*(1.-RP)*(-SP+2.*TP-RP-1.)
     4         -YY(5)*(1.-SP)*(1.+RP)*(-SP-2.*TP+RP-1.)
     5         -YY(6)*(1.+SP)*(1.+RP)*(SP-2.*TP+RP-1.)
     6         +YY(7)*(1.+SP)*(1.+RP)*(SP+2.*TP+RP-1.)
     7         +YY(8)*(1.-SP)*(1.+RP)*(-SP+2.*TP+RP-1.))
      BB=.250*(-YY(9)*(1.-SP**2)*(1.-RP)
     1         +YY(10)*(1.+SP)*(-2.*TP)*(1.-RP)
     2         +YY(11)*(1.-SP**2)*(1.-RP)
     3         +YY(12)*(1.-SP)*(-2.*TP)*(1.-RP)
     4         -YY(13)*(1.-SP)*(1.-RP**2)
     5         -YY(14)*(1.+SP)*(1.-RP**2)
     6         +YY(15)*(1.+SP)*(1.-RP**2)
     7         +YY(16)*(1.-SP)*(1.-RP**2)
     8         -YY(17)*(1.-SP**2)*(1.+RP)
     9         +YY(18)*(1.+SP)*(-2.*TP)*(1.+RP)
     #         +YY(19)*(1.-SP**2)*(1.+RP)
     1         +YY(20)*(1.-SP)*(-2.*TP)*(1.+RP))
      A22=AA+BB
      AA=.125*(-YY(1)*(1.-SP)*(1.-TP)*(-SP-TP-2.*RP-1.)
     1         -YY(2)*(1.+SP)*(1.-TP)*(SP-TP-2.*RP-1.)
     2         -YY(3)*(1.+SP)*(1.+TP)*(SP+TP-2.*RP-1.)
     3         -YY(4)*(1.-SP)*(1.+TP)*(-SP+TP-2.*RP-1.)
     4         +YY(5)*(1.-SP)*(1.-TP)*(-SP-TP+2.*RP-1.)
     5         +YY(6)*(1.+SP)*(1.-TP)*(SP-TP+2.*RP-1.)
     6         +YY(7)*(1.+SP)*(1.+TP)*(SP+TP+2.*RP-1.)
     7         +YY(8)*(1.-SP)*(1.+TP)*(-SP+TP+2.*RP-1.))
      BB=.250*(-YY(9)*(1.-SP**2)*(1.-TP)
     1         -YY(10)*(1.+SP)*(1.-TP**2)
     2         -YY(11)*(1.-SP**2)*(1.+TP)
     3         -YY(12)*(1.-SP)*(1.-TP**2)
     4         +YY(13)*(1.-SP)*(1.-TP)*(-2.*RP)
     5         +YY(14)*(1.+SP)*(1.-TP)*(-2.*RP)
     6         +YY(15)*(1.+SP)*(1.+TP)*(-2.*RP)
     7         +YY(16)*(1.-SP)*(1.+TP)*(-2.*RP)
     8         +YY(17)*(1.-SP**2)*(1.-TP)
     9         +YY(18)*(1.+SP)*(1.-TP**2)
     #         +YY(19)*(1.-SP**2)*(1.+TP)
     1         +YY(20)*(1.-SP)*(1.-TP**2))
      A23=AA+BB
C
      AA=.125*(-ZZ(1)*(1.-TP)*(1.-RP)*(-2.*SP-TP-RP-1.)
     1         +ZZ(2)*(1.-TP)*(1.-RP)*(2.*SP-TP-RP-1.)
     2         +ZZ(3)*(1.+TP)*(1.-RP)*(-2.*SP+TP-RP-1.)
     3         -ZZ(4)*(1.+TP)*(1.-RP)*(-2.*SP+TP-RP-1.)
     4         -ZZ(5)*(1.-TP)*(1.+RP)*(-2.*SP-TP+RP-1.)
     5         +ZZ(6)*(1.-TP)*(1.+RP)*(2.*SP-TP+RP-1.)
     6         +ZZ(7)*(1.+TP)*(1.+RP)*(2.*SP+TP+RP-1.)
     7         -ZZ(8)*(1.+TP)*(1.+RP)*(-2.*SP+TP+RP-1.))
      BB=.250*(+ZZ(9)*(-2.*SP)*(1.-TP)*(1.-RP)
     1         +ZZ(10)*(1.-TP**2)*(1.-RP)
     2         +ZZ(11)*(-2.*SP)*(1.+TP)*(1.-RP)
     3         -ZZ(12)*(1.-TP**2)*(1.-RP)
     4         -ZZ(13)*(1.-TP)*(1.-RP**2)
     5         +ZZ(14)*(1.-TP)*(1.-RP**2)
     6         +ZZ(15)*(1.+TP)*(1.-RP**2)
     7         -ZZ(16)*(1.+TP)*(1.-RP**2)
     8         +ZZ(17)*(-2.*SP)*(1.-TP)*(1.+RP)
     9         +ZZ(18)*(1.-TP**2)*(1.+RP)
     #         +ZZ(19)*(-2.*SP)*(1.+TP)*(1.+RP)
     1         -ZZ(20)*(1.-TP**2)*(1.+RP))
      A31=AA+BB
      AA=.125*(-ZZ(1)*(1.-SP)*(1.-RP)*(-SP-2.*TP-RP-1.)
     1         -ZZ(2)*(1.+SP)*(1.-RP)*(SP-2.*TP-RP-1.)
     2         +ZZ(3)*(1.+SP)*(1.-RP)*(SP+2.*TP-RP-1.)
     3         +ZZ(4)*(1.-SP)*(1.-RP)*(-SP+2.*TP-RP-1.)
     4         -ZZ(5)*(1.-SP)*(1.+RP)*(-SP-2.*TP+RP-1.)
     5         -ZZ(6)*(1.+SP)*(1.+RP)*(SP-2.*TP+RP-1.)
     6         +ZZ(7)*(1.+SP)*(1.+RP)*(SP+2.*TP+RP-1.)
     7         +ZZ(8)*(1.-SP)*(1.+RP)*(-SP+2.*TP+RP-1.))
      BB=.250*(-ZZ(9)*(1.-SP**2)*(1.-RP)
     1         +ZZ(10)*(1.+SP)*(-2.*TP)*(1.-RP)
     2         +ZZ(11)*(1.-SP**2)*(1.-RP)
     3         +ZZ(12)*(1.-SP)*(-2.*TP)*(1.-RP)
     4         -ZZ(13)*(1.-SP)*(1.-RP**2)
     5         -ZZ(14)*(1.+SP)*(1.-RP**2)
     6         +ZZ(15)*(1.+SP)*(1.-RP**2)
     7         +ZZ(16)*(1.-SP)*(1.-RP**2)
     8         -ZZ(17)*(1.-SP**2)*(1.+RP)
     9         +ZZ(18)*(1.+SP)*(-2.*TP)*(1.+RP)
     #         +ZZ(19)*(1.-SP**2)*(1.+RP)
     1         +ZZ(20)*(1.-SP)*(-2.*TP)*(1.+RP))
      A32=AA+BB
      AA=.125*(-ZZ(1)*(1.-SP)*(1.-TP)*(-SP-TP-2.*RP-1.)
     1         -ZZ(2)*(1.+SP)*(1.-TP)*(SP-TP-2.*RP-1.)
     2         -ZZ(3)*(1.+SP)*(1.+TP)*(SP+TP-2.*RP-1.)
     3         -ZZ(4)*(1.-SP)*(1.+TP)*(-SP+TP-2.*RP-1.)
     4         +ZZ(5)*(1.-SP)*(1.-TP)*(-SP-TP+2.*RP-1.)
     5         +ZZ(6)*(1.+SP)*(1.-TP)*(SP-TP+2.*RP-1.)
     6         +ZZ(7)*(1.+SP)*(1.+TP)*(SP+TP+2.*RP-1.)
     7         +ZZ(8)*(1.-SP)*(1.+TP)*(-SP+TP+2.*RP-1.))
      BB=.250*(-ZZ(9)*(1.-SP**2)*(1.-TP)
     1         -ZZ(10)*(1.+SP)*(1.-TP**2)
     2         -ZZ(11)*(1.-SP**2)*(1.+TP)
     3         -ZZ(12)*(1.-SP)*(1.-TP**2)
     4         +ZZ(13)*(1.-SP)*(1.-TP)*(-2.*RP)
     5         +ZZ(14)*(1.+SP)*(1.-TP)*(-2.*RP)
     6         +ZZ(15)*(1.+SP)*(1.+TP)*(-2.*RP)
     7         +ZZ(16)*(1.-SP)*(1.+TP)*(-2.*RP)
     8         +ZZ(17)*(1.-SP**2)*(1.-TP)
     9         +ZZ(18)*(1.+SP)*(1.-TP**2)
     #         +ZZ(19)*(1.-SP**2)*(1.+TP)
     1         +ZZ(20)*(1.-SP)*(1.-TP**2))
      A33=AA+BB
      RETURN
C
C     27-NODE HEX
C
  210 CONTINUE
      AA=.125*(-XX(1)*SP*TP*RP*(1.-SP)*(1.-TP)*(1.-RP)
     1         +XX(2)*SP*TP*RP*(1.+SP)*(1.-TP)*(1.-RP)
     2         -XX(3)*SP*TP*RP*(1.+SP)*(1.+TP)*(1.-RP)
     3         +XX(4)*SP*TP*RP*(1.-SP)*(1.+TP)*(1.-RP)
     4         +XX(5)*SP*TP*RP*(1.-SP)*(1.-TP)*(1.+RP)
     5         -XX(6)*SP*TP*RP*(1.+SP)*(1.-TP)*(1.+RP)
     6         +XX(7)*SP*TP*RP*(1.+SP)*(1.+TP)*(1.+RP)
     7         -XX(8)*SP*TP*RP*(1.-SP)*(1.+TP)*(1.+RP))
      BB=.250*(+XX (9)*TP*RP*(1.-SP**2)*(1.-TP)*(1.-RP)
     1         -XX(10)*SP*RP*(1.+SP)*(1.-TP**2)*(1.-RP)
     2         -XX(11)*TP*RP*(1.-SP**2)*(1.+TP)*(1.-RP)
     3         +XX(12)*SP*RP*(1.-SP)*(1.-TP**2)*(1.-RP)
     4         +XX(13)*SP*TP*(1.-SP)*(1.-TP)*(1.-RP**2)
     5         -XX(14)*SP*TP*(1.+SP)*(1.-TP)*(1.-RP**2)
     6         +XX(15)*SP*TP*(1.+SP)*(1.+TP)*(1.-RP**2)
     7         -XX(16)*SP*TP*(1.-SP)*(1.+TP)*(1.-RP**2))
      CC=.250*(-XX(17)*TP*RP*(1.-SP**2)*(1.-TP)*(1.+RP)
     1         +XX(18)*SP*RP*(1.+SP)*(1.-TP**2)*(1.+RP)
     2         +XX(19)*TP*RP*(1.-SP**2)*(1.+TP)*(1.+RP)
     3         -XX(20)*SP*RP*(1.-SP)*(1.-TP**2)*(1.+RP))+
     &   .500*(-XX(21)*TP*(1.-SP**2)*(1.-TP)*(1.-RP**2) +
     &          XX(22)*SP*(1.+SP)*(1.-TP**2)*(1.-RP**2) +
     &          XX(23)*TP*(1.-SP**2)*(1.+TP)*(1.-RP**2) -
     &          XX(24)*SP*(1.-SP)*(1.-TP**2)*(1.-RP**2) -
     &          XX(25)*RP*(1.-SP**2)*(1.-TP**2)*(1.-RP) +
     &          XX(26)*RP*(1.-SP**2)*(1.-TP**2)*(1.+RP))+
     &          XX(27)*(1.-SP**2)*(1.-TP**2)*(1.-RP**2)
      F1=AA+BB+CC
C
      AA=.125*(-YY(1)*SP*TP*RP*(1.-SP)*(1.-TP)*(1.-RP)
     1         +YY(2)*SP*TP*RP*(1.+SP)*(1.-TP)*(1.-RP)
     2         -YY(3)*SP*TP*RP*(1.+SP)*(1.+TP)*(1.-RP)
     3         +YY(4)*SP*TP*RP*(1.-SP)*(1.+TP)*(1.-RP)
     4         +YY(5)*SP*TP*RP*(1.-SP)*(1.-TP)*(1.+RP)
     5         -YY(6)*SP*TP*RP*(1.+SP)*(1.-TP)*(1.+RP)
     6         +YY(7)*SP*TP*RP*(1.+SP)*(1.+TP)*(1.+RP)
     7         -YY(8)*SP*TP*RP*(1.-SP)*(1.+TP)*(1.+RP))
      BB=.250*(+YY (9)*TP*RP*(1.-SP**2)*(1.-TP)*(1.-RP)
     1         -YY(10)*SP*RP*(1.+SP)*(1.-TP**2)*(1.-RP)
     2         -YY(11)*TP*RP*(1.-SP**2)*(1.+TP)*(1.-RP)
     3         +YY(12)*SP*RP*(1.-SP)*(1.-TP**2)*(1.-RP)
     4         +YY(13)*SP*TP*(1.-SP)*(1.-TP)*(1.-RP**2)
     5         -YY(14)*SP*TP*(1.+SP)*(1.-TP)*(1.-RP**2)
     6         +YY(15)*SP*TP*(1.+SP)*(1.+TP)*(1.-RP**2)
     7         -YY(16)*SP*TP*(1.-SP)*(1.+TP)*(1.-RP**2))
      CC=.250*(-YY(17)*TP*RP*(1.-SP**2)*(1.-TP)*(1.+RP)
     1         +YY(18)*SP*RP*(1.+SP)*(1.-TP**2)*(1.+RP)
     2         +YY(19)*TP*RP*(1.-SP**2)*(1.+TP)*(1.+RP)
     3         -YY(20)*SP*RP*(1.-SP)*(1.-TP**2)*(1.+RP))+
     &   .500*(-YY(21)*TP*(1.-SP**2)*(1.-TP)*(1.-RP**2) +
     &          YY(22)*SP*(1.+SP)*(1.-TP**2)*(1.-RP**2) +
     &          YY(23)*TP*(1.-SP**2)*(1.+TP)*(1.-RP**2) -
     &          YY(24)*SP*(1.-SP)*(1.-TP**2)*(1.-RP**2) -
     &          YY(25)*RP*(1.-SP**2)*(1.-TP**2)*(1.-RP) +
     &          YY(26)*RP*(1.-SP**2)*(1.-TP**2)*(1.+RP))+
     &          YY(27)*(1.-SP**2)*(1.-TP**2)*(1.-RP**2)
      F2=AA+BB+CC
      AA=.125*(-ZZ(1)*SP*TP*RP*(1.-SP)*(1.-TP)*(1.-RP)
     1         +ZZ(2)*SP*TP*RP*(1.+SP)*(1.-TP)*(1.-RP)
     2         -ZZ(3)*SP*TP*RP*(1.+SP)*(1.+TP)*(1.-RP)
     3         +ZZ(4)*SP*TP*RP*(1.-SP)*(1.+TP)*(1.-RP)
     4         +ZZ(5)*SP*TP*RP*(1.-SP)*(1.-TP)*(1.+RP)
     5         -ZZ(6)*SP*TP*RP*(1.+SP)*(1.-TP)*(1.+RP)
     6         +ZZ(7)*SP*TP*RP*(1.+SP)*(1.+TP)*(1.+RP)
     7         -ZZ(8)*SP*TP*RP*(1.-SP)*(1.+TP)*(1.+RP))
      BB=.250*(+ZZ (9)*TP*RP*(1.-SP**2)*(1.-TP)*(1.-RP)
     1         -ZZ(10)*SP*RP*(1.+SP)*(1.-TP**2)*(1.-RP)
     2         -ZZ(11)*TP*RP*(1.-SP**2)*(1.+TP)*(1.-RP)
     3         +ZZ(12)*SP*RP*(1.-SP)*(1.-TP**2)*(1.-RP)
     4         +ZZ(13)*SP*TP*(1.-SP)*(1.-TP)*(1.-RP**2)
     5         -ZZ(14)*SP*TP*(1.+SP)*(1.-TP)*(1.-RP**2)
     6         +ZZ(15)*SP*TP*(1.+SP)*(1.+TP)*(1.-RP**2)
     7         -ZZ(16)*SP*TP*(1.-SP)*(1.+TP)*(1.-RP**2))
      CC=.250*(-ZZ(17)*TP*RP*(1.-SP**2)*(1.-TP)*(1.+RP)
     1         +ZZ(18)*SP*RP*(1.+SP)*(1.-TP**2)*(1.+RP)
     2         +ZZ(19)*TP*RP*(1.-SP**2)*(1.+TP)*(1.+RP)
     3         -ZZ(20)*SP*RP*(1.-SP)*(1.-TP**2)*(1.+RP))+
     &   .500*(-ZZ(21)*TP*(1.-SP**2)*(1.-TP)*(1.-RP**2) +
     &          ZZ(22)*SP*(1.+SP)*(1.-TP**2)*(1.-RP**2) +
     &          ZZ(23)*TP*(1.-SP**2)*(1.+TP)*(1.-RP**2) -
     &          ZZ(24)*SP*(1.-SP)*(1.-TP**2)*(1.-RP**2) -
     &          ZZ(25)*RP*(1.-SP**2)*(1.-TP**2)*(1.-RP) +
     &          ZZ(26)*RP*(1.-SP**2)*(1.-TP**2)*(1.+RP))+
     &          ZZ(27)*(1.-SP**2)*(1.-TP**2)*(1.-RP**2)
      F3=AA+BB+CC
C
C     **** SHAPE FUNCTION DERIVATIVES FOR 27 NODE HEX GO HERE
C
      RETURN
C
      END
