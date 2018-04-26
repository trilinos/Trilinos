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

C $Id: mp3pt.f,v 1.1 1993/07/16 16:47:05 gdsjaar Exp $ 
C $Log: mp3pt.f,v $
C Revision 1.1  1993/07/16 16:47:05  gdsjaar
C Changed plt to library rather than single source file.
C 
C=======================================================================
      SUBROUTINE MP3PT(N,X0,Y0,Z0,PX,PY,PZ,MASK)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      LOGICAL NOCLIP
      DIMENSION PX(*),PY(*),PZ(*),MASK(*),C1(3),C2(3),X0(*),Y0(*),Z0(*)

      VCX = (VWPORT(1)+VWPORT(2))/2.
      VSX = (VWPORT(2)-VWPORT(1))/2.
      VCY = (VWPORT(3)+VWPORT(4))/2.
      VSY = (VWPORT(4)-VWPORT(3))/2.
      NOCLIP = (NCPLAN.EQ.0)
      KM = 0
      J = 0
 2260 IF (.NOT. (J.LT.N)) GO TO 2270
      JN = MIN(N-J,32)
      J1 = J
      KM = KM + 1
      J = J + JN
      MASK(KM) = -1
      IF (NOCLIP) THEN
         CALL MPMUL3(JN,X0(1+J1),Y0(1+J1),Z0(1+J1),MVP,TARR1,TARR2,
     *               TARR3,TARR4)

      ELSE
         CALL MPMUL3(JN,X0(1+J1),Y0(1+J1),Z0(1+J1),MODEL,TARR1,TARR2,
     *               TARR3,TARR4)
         DO 2280 I = 1,NCPLAN
            C1(1) = CPPLAN(1,1,I)
            C1(2) = CPPLAN(1,2,I)
            C1(3) = CPPLAN(1,3,I)
            C2(1) = CPPLAN(2,1,I)
            C2(2) = CPPLAN(2,2,I)
            C2(3) = CPPLAN(2,3,I)
            CALL PLTCP3(JN,MASK(KM),TARR1,TARR2,TARR3,C1,C2)
 2280    CONTINUE
         CALL MPMUL4(JN,MASK(KM),TARR1,TARR2,TARR3,TARR4,VP,TARR1,TARR2,
     *               TARR3,TARR4)
      END IF

      CALL PLTZCP(CPNEAR,CPFAR,JN,MASK(KM),TARR4)
      DO 2300 I = 1,JN
         TARR1(I) = (TARR1(I)/TARR4(I))*VSX + VCX
         TARR2(I) = (TARR2(I)/TARR4(I))*VSY + VCY
         TARR3(I) = TARR3(I)/TARR4(I)
 2300 CONTINUE
      C1(1) = VWPORT(1)
      C1(2) = VWPORT(3)
      C2(1) = VWPORT(2)
      C2(2) = VWPORT(4)
      CALL PLTVWP(C1,C2,JN,MASK(KM),TARR1,TARR2)
      DO 2320 I = 1,JN
         PX(J1+I) = TARR1(I)
         PY(J1+I) = TARR2(I)
         PZ(J1+I) = TARR3(I)
 2320 CONTINUE
      GO TO 2260

 2270 CONTINUE
      RETURN

      END
