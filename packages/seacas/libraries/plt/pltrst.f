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

C $Id: pltrst.f,v 1.1 1993/07/16 16:49:23 gdsjaar Exp $ 
C $Log: pltrst.f,v $
C Revision 1.1  1993/07/16 16:49:23  gdsjaar
C Changed plt to library rather than single source file.
C 
C=======================================================================
      SUBROUTINE PLTRST
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

      IF ((DEVCAP(15)+DEVCAP(16))/2..LE.800.) THEN
         CALL PLTSTT(12,2.)

      ELSE
         CALL PLTSTT(12,3.)
      END IF

      BUFF = 1./58.
      CALL VDSTCS(BUFF)
      TEXTP(35) = DEFOUT(6)
      TEXTP(36) = DEFOUT(7)
      TEXTP(1) = .015
      TEXTP(2) = 0.
      TEXTP(3) = 0.
      DO 2340 I = 4,11
         TEXTP(I) = 0.
 2340 CONTINUE
      TEXTP(20) = 0.
      TEXTP(21) = 0.
      TEXTP(22) = 1.
      TEXTP(23) = 0.
      TEXTP(24) = TEXTP(22)
      TEXTP(25) = .75
      TEXTP(26) = 0.
      TEXTP(27) = TEXTP(25)
      TEXTP(30) = 1.5
      TEXTP(31) = 1.
      TEXTP(32) = .7
      TEXTP(33) = .8
      TEXTP(34) = .5
      TEXTP(37) = 160.
      CALL PLTITM
      RETURN

      END
