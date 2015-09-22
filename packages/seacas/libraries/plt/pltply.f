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

C $Id: pltply.f,v 1.1 1993/07/16 16:49:07 gdsjaar Exp $ 
C $Log: pltply.f,v $
C Revision 1.1  1993/07/16 16:49:07  gdsjaar
C Changed plt to library rather than single source file.
C 
C=======================================================================
      SUBROUTINE PLTPLY(N,XA,YA)
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
      DIMENSION XA(*),YA(*)
      LOGICAL STATUS,MEMFRE
      REAL MEMORY

      IF (N.LE.0) THEN
         RETURN

      END IF

      IVX = MEMALL(N,MEMORY)
      IVY = MEMALL(N,MEMORY)
      IF (MAPP(10).EQ.1.) THEN
         DO 2080 I = 1,N
            IMX = IVX + I - 1
            IMY = IVY + I - 1
            CALL PLTP2D(XA(I),YA(I),MEMORY(IMX),MEMORY(IMY))
 2080    CONTINUE

      ELSE
         DO 2100 I = 1,N
            IMX = IVX + I - 1
            IMY = IVY + I - 1
            MEMORY(IMX) = XA(I)
            MEMORY(IMY) = YA(I)
 2100    CONTINUE
      END IF

      ICX = MEMALL(N+10,MEMORY)
      ICY = MEMALL(N+10,MEMORY)
      IZ1 = MEMALL(N,MEMORY)
      IZ2 = MEMALL(N+10,MEMORY)
      NO = N + 10
      CALL PLTVWG(MAPP(6),MAPP(8),N,MEMORY(IVX),MEMORY(IVY),MEMORY(IZ1),
     *            NO,MEMORY(ICX),MEMORY(ICY),MEMORY(IZ2))
      IF (NO.GT.0) THEN
         CALL VDPOLY(MEMORY(ICX),MEMORY(ICY),NO)
      END IF

      STATUS = MEMFRE(ICX,MEMORY)
      STATUS = MEMFRE(ICY,MEMORY)
      STATUS = MEMFRE(IZ1,MEMORY)
      STATUS = MEMFRE(IZ2,MEMORY)
      STATUS = MEMFRE(IVX,MEMORY)
      STATUS = MEMFRE(IVY,MEMORY)
      RETURN

      END
