C Copyright (C) 2009-2017 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
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
C     * Neither the name of NTESS nor the names of its
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

C $Id: mpview.f,v 1.1 1993/07/16 16:47:31 gdsjaar Exp $ 
C $Log: mpview.f,v $
C Revision 1.1  1993/07/16 16:47:31  gdsjaar
C Changed plt to library rather than single source file.
C 
C=======================================================================
      LOGICAL FUNCTION MPVIEW(LEFT,RIGHT,BOTTOM,TOP)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
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
      REAL LEFT
      CHARACTER*6 SUBNAM
      PARAMETER (SUBNAM='MPVIEW')

      MPVIEW = .FALSE.
      IF (LEFT.GE.RIGHT) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,
     *  'You cannot specify the left viewport edge >= to the right edge'
     *               ,2)
         RETURN

      END IF

      IF (TOP.LE.BOTTOM) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,
     *  'You cannot specify the top viewport edge <= to the bottom edge'
     *               ,2)
         RETURN

      END IF

      IF (TOP.LT.0 .OR. TOP.GT.DEVP(5)) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,'Top viewport specification out of range',2)
         RETURN

      END IF

      IF (BOTTOM.LT.0 .OR. BOTTOM.GT.DEVP(5)) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,'Bottom viewport specification out of range'
     *               ,2)
         RETURN

      END IF

      IF (LEFT.LT.0 .OR. LEFT.GT.DEVP(4)) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,'Left viewport specification out of range',
     *               2)
         RETURN

      END IF

      IF (RIGHT.LT.0 .OR. RIGHT.GT.DEVP(4)) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,'Right viewport specification out of range',
     *               2)
         RETURN

      END IF

      MPVIEW = .TRUE.
      VWPORT(1) = LEFT
      VWPORT(2) = RIGHT
      VWPORT(3) = BOTTOM
      VWPORT(4) = TOP
      RETURN

      END
