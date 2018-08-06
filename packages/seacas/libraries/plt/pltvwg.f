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

C $Id: pltvwg.f,v 1.1 1993/07/16 16:49:49 gdsjaar Exp $ 
C $Log: pltvwg.f,v $
C Revision 1.1  1993/07/16 16:49:49  gdsjaar
C Changed plt to library rather than single source file.
C 
C=======================================================================
      SUBROUTINE PLTVWG(PLL,PUR,N,XV,YV,ZV,NO,XVO,YVO,ZVO)
      DIMENSION XV(*),YV(*),ZV(*),XVO(*),YVO(*),ZVO(*),PLL(*),PUR(*)
      DIMENSION XVT(50),YVT(50),QUR(2),QLL(2)

      QUR(1) = PUR(1) + .0001
      QUR(2) = PUR(2) + .0001
      QLL(1) = PLL(1) - .0001
      QLL(2) = PLL(2) - .0001
      NOSAVE = NO
      NT = 50
      CALL PLTLI1(QLL,QUR,N,XV,YV,NT,XVT,YVT)
      NO = NOSAVE
      CALL PLTLI2(QLL,QUR,NT,XVT,YVT,NO,XVO,YVO)
      NT = 50
      CALL PLTLI3(QLL,QUR,NO,XVO,YVO,NT,XVT,YVT)
      NO = NOSAVE
      CALL PLTLI4(QLL,QUR,NT,XVT,YVT,NO,XVO,YVO)
      DO 2180 I = 1,NO
         ZVO(I) = PLTPGZ(NO,XV,YV,ZV,XVO(I),YVO(I))
 2180 CONTINUE
      RETURN

      END
