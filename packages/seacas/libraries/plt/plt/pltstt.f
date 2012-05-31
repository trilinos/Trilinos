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

C $Id: pltstt.f,v 1.1 1993/07/16 16:49:36 gdsjaar Exp $ 
C $Log: pltstt.f,v $
C Revision 1.1  1993/07/16 16:49:36  gdsjaar
C Changed plt to library rather than single source file.
C 
C=======================================================================
      LOGICAL FUNCTION PLTSTT(INDX,BUFF)
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
      DIMENSION BUFF(*),DUMMY(7)
      CHARACTER*16 IERROR

      PLTSTT = .TRUE.
      IF (INDX.EQ.0) THEN
         CALL PLTRST

      ELSE IF (INDX.EQ.1) THEN
         CALL VDSTCS(BUFF(1))
         CALL VDIQOS(DUMMY)
         TEXTP(35) = DUMMY(6)
         TEXTP(36) = DUMMY(7)

      ELSE IF (INDX.EQ.2) THEN
         TEXTP(1) = BUFF(1)
         CALL PLTITM

      ELSE IF (INDX.EQ.3) THEN
         TEXTP(2) = BUFF(1)
         CALL PLTITM

      ELSE IF (INDX.EQ.4) THEN
         TEXTP(3) = BUFF(1)
         CALL PLTITM

      ELSE IF (INDX.EQ.5) THEN
         DO 2300 I = 20,27
            TEXTP(I) = BUFF(I-19)
 2300    CONTINUE

      ELSE IF (INDX.EQ.6) THEN
         TEXTP(30) = BUFF(1)

      ELSE IF (INDX.EQ.7) THEN
         TEXTP(31) = BUFF(1)

      ELSE IF (INDX.EQ.8) THEN
         TEXTP(32) = BUFF(1)

      ELSE IF (INDX.EQ.9) THEN
         TEXTP(33) = BUFF(1)

      ELSE IF (INDX.EQ.10) THEN
         TEXTP(34) = BUFF(1)

      ELSE IF (INDX.EQ.11) THEN
         TEXTP(37) = BUFF(1)

      ELSE IF (INDX.EQ.12) THEN
         IF (BUFF(1).EQ.1.) THEN
            CALL PLTFNT('ROMFNT')
            TEXTP(38) = 1.

         ELSE IF (BUFF(1).EQ.2.) THEN
            CALL PLTFNT('STKFNT')
            TEXTP(38) = 2.

         ELSE IF (BUFF(1).EQ.3.) THEN
            CALL PLTFNT('SSRFNT')
            TEXTP(38) = 3.

         ELSE
            CALL CHRRVC(BUFF(1),IERROR,L)
            CALL PLTFLU
            CALL SIORPT('PLTSTT','Illegal buffer '//IERROR(1:L)//
     *                  ' passed to PLTSTT.',2)
            PLTSTT = .FALSE.
         END IF

      ELSE
         CALL CHRIC(INDX,IERROR,L)
         CALL PLTFLU
         CALL SIORPT('PLTSTT','Illegal index '//IERROR(1:L)//'.',2)
         PLTSTT = .FALSE.
         RETURN

      END IF

      RETURN

      END
