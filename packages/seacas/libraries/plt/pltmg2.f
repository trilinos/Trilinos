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

C $Id: pltmg2.f,v 1.1 1993/07/16 16:48:47 gdsjaar Exp $ 
C $Log: pltmg2.f,v $
C Revision 1.1  1993/07/16 16:48:47  gdsjaar
C Changed plt to library rather than single source file.
C 
C=======================================================================
      SUBROUTINE PLTMG2(MAP,N,XV,YV,NO,XVO,YVO)
      REAL MAP(*)
      INTEGER N
      REAL XV(*),YV(*)
      INTEGER NO
      REAL XVO(*),YVO(*)
      REAL XWORK(50),YWORK(50)
      INTEGER NWORK

      NOSAVE = NO
      AXX = MAP(1)
      AYY = MAP(4)
      AXY = MAP(3)
      AYX = MAP(2)
      BX = MAP(5)
      BY = MAP(6)
      DO 2220 I = 1,N
         XVO(I) = AXX*XV(I) + AXY*YV(I) + BX
         YVO(I) = AYX*XV(I) + AYY*YV(I) + BY
 2220 CONTINUE
      NWORK = 50
      CALL PLTCG2(N,XVO,YVO,NWORK,XWORK,YWORK,MAP(7),MAP(9))
      NO = NOSAVE
      CALL PLTCG2(NWORK,XWORK,YWORK,NO,XVO,YVO,MAP(9),MAP(11))
      NWORK = 50
      CALL PLTCG2(NO,XVO,YVO,NWORK,XWORK,YWORK,MAP(11),MAP(13))
      NO = NOSAVE
      CALL PLTCG2(NWORK,XWORK,YWORK,NO,XVO,YVO,MAP(13),MAP(7))
      RETURN

      END
