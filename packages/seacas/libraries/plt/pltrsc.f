C Copyright (C) 2009 Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software
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

C $Id: pltrsc.f,v 1.1 1993/07/16 16:49:19 gdsjaar Exp $ 
C $Log: pltrsc.f,v $
C Revision 1.1  1993/07/16 16:49:19  gdsjaar
C Changed plt to library rather than single source file.
C 
C=======================================================================
      SUBROUTINE PLTRSC
      REAL DEVCAP(23)
      REAL DEFOUT(7)
      COMMON /STATUS/DEVCAP,DEFOUT
      REAL DEVP(5)
      COMMON /DEVICE/DEVP
      REAL COLP(3)
      REAL PALETT(3,16)
      COMMON /COLOR/COLP,PALETT
      REAL TEXTP(40)
      COMMON /TEXT/TEXTP
      REAL VECTP(5)
      REAL XCUR
      REAL YCUR
      COMMON /VECTRC/VECTP,XCUR,YCUR
      INTEGER IDEX(200,2)
      INTEGER NVECT(200,2)
      REAL XSIZE(200,2)
      REAL YSIZE(200,2)
      REAL X0(2300,2)
      REAL Y0(2300,2)
      REAL X1(2300,2)
      REAL Y1(2300,2)
      COMMON /FONT/IDEX,NVECT,XSIZE,YSIZE,X0,Y0,X1,Y1
      REAL GRAPHP(100)
      COMMON /GRAPH/GRAPHP
      COMMON /MAPPAR/MAPP(11)
      REAL MAPP
      COMMON /STORAG/MEMORY(1000)

      COLP(2) = MIN(16.,DEVCAP(4))
      COLP(3) = DEVCAP(4) - COLP(2)
      PALETT(1,1) = 0.
      PALETT(2,1) = 0.
      PALETT(3,1) = 0.
      CALL VDSTCO(1,0,PALETT(1,1),0)
      PALETT(1,2) = 1.
      PALETT(2,2) = 0.
      PALETT(2,3) = 0.
      CALL VDSTCO(1,1,PALETT(1,2),0)
      PALETT(1,3) = 0.
      PALETT(2,3) = 1.
      PALETT(3,3) = 0.
      CALL VDSTCO(1,2,PALETT(1,3),0)
      PALETT(1,4) = 1.
      PALETT(2,4) = 1.
      PALETT(3,4) = 0.
      CALL VDSTCO(1,3,PALETT(1,4),0)
      PALETT(1,5) = 0.
      PALETT(2,5) = 0.
      PALETT(3,5) = 1.
      CALL VDSTCO(1,4,PALETT(1,5),0)
      PALETT(1,6) = 1.
      PALETT(2,6) = 0.
      PALETT(3,6) = 1.
      CALL VDSTCO(1,5,PALETT(1,6),0)
      PALETT(1,7) = 0.
      PALETT(2,7) = 1.
      PALETT(3,7) = 1.
      CALL VDSTCO(1,6,PALETT(1,7),0)
      PALETT(1,8) = 1.
      PALETT(2,8) = 1.
      PALETT(3,8) = 1.
      CALL VDSTCO(1,7,PALETT(1,8),0)
      PALETT(1,9) = .4
      PALETT(2,9) = .4
      PALETT(3,9) = .4
      CALL VDSTCO(1,8,PALETT(1,9),0)
      PALETT(1,10) = .7
      PALETT(2,10) = .7
      PALETT(3,10) = .7
      CALL VDSTCO(1,9,PALETT(1,10),0)
      PALETT(1,11) = .225
      PALETT(2,11) = .225
      PALETT(3,11) = .225
      CALL VDSTCO(1,10,PALETT(1,11),0)
      PALETT(1,12) = 1.
      PALETT(2,12) = .35
      PALETT(3,12) = .45
      CALL VDSTCO(1,11,PALETT(1,12),0)
      PALETT(1,13) = .5
      PALETT(2,13) = 1.
      PALETT(3,13) = .8
      CALL VDSTCO(1,12,PALETT(1,13),0)
      PALETT(1,14) = .4
      PALETT(2,14) = .7
      PALETT(3,14) = 1.
      CALL VDSTCO(1,13,PALETT(1,14),0)
      PALETT(1,15) = .706
      PALETT(2,15) = 0.
      PALETT(3,15) = .706
      CALL VDSTCO(1,14,PALETT(1,15),0)
      PALETT(1,16) = 1.
      PALETT(2,16) = .659
      PALETT(3,16) = 0.
      CALL VDSTCO(1,15,PALETT(1,16),0)
      CALL PLTSTC(1,0.)
      RETURN

      END
