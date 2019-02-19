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

C $Id: pltmv2.f,v 1.1 1993/07/16 16:48:57 gdsjaar Exp $
C $Log: pltmv2.f,v $
C Revision 1.1  1993/07/16 16:48:57  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      SUBROUTINE PLTMV2(UMAP,N,MASK,PX,PY,QX,QY,PPX,PPY,QQX,QQY)
      DIMENSION UMAP(*),MASK(*),PX(*),PY(*),QX(*),QY(*),PPX(*),PPY(*),
     *          QQX(*),QQY(*)

      AXX = UMAP(1)
      AYY = UMAP(4)
      AXY = UMAP(3)
      AYX = UMAP(2)
      BX = UMAP(5)
      BY = UMAP(6)
      DO 2120 I = 1,N
         PPX(I) = AXX*PX(I) + AXY*PY(I) + BX
         QQX(I) = AXX*QX(I) + AXY*QY(I) + BX
         PPY(I) = AYX*PX(I) + AYY*PY(I) + BY
         QQY(I) = AYX*QX(I) + AYY*QY(I) + BY
 2120 CONTINUE
      J = 0
 2140 IF (.NOT. (J.LT.N)) GO TO 2150
      JN = MIN(N-J,32)
      J1 = J + 1
      KM = 1 + J/32
      J = J + JN
      CALL PLTVWV(UMAP(7),UMAP(11),JN,MASK(KM),PPX(J1),PPY(J1),QQX(J1),
     *            QQY(J1))
      GO TO 2140

 2150 CONTINUE
      RETURN

      END
