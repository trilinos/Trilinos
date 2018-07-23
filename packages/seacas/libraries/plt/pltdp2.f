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

C $Id: pltdp2.f,v 1.2 2000/10/25 18:55:02 gdsjaar Exp $ 
C $Log: pltdp2.f,v $
C Revision 1.2  2000/10/25 18:55:02  gdsjaar
C In the pltli? functions, check for N==0 before doing any array
C accesses.
C
C Also changed all references to 'mask' to be arrays where they were
C scalars since downstream code seems to treat them as arrays.
C
C Revision 1.1  1993/07/16 16:48:01  gdsjaar
C Changed plt to library rather than single source file.
C 
C=======================================================================
      SUBROUTINE PLTDP2(MAP,N,PX,PY)
      REAL MAP(*),PX(*),PY(*)
      REAL XWORK(32),YWORK(32)
      INTEGER MASK(1)

      J = 0
 2360 IF (.NOT. (J.LT.N)) GO TO 2370
      JN = MIN(N-J,32)
      J1 = J
      J = J + JN

      do 2400 i=1, jn
        XWORK(I) = MAP(1)*PX(J1+I) + MAP(3)*PY(J1+I) + MAP(5)
        YWORK(I) = MAP(2)*PX(J1+I) + MAP(4)*PY(J1+I) + MAP(6)
 2400 CONTINUE

      MASK(1) = -1
      CALL PLTCP2(JN,MASK,XWORK,YWORK,MAP(7),MAP(9))
      CALL PLTCP2(JN,MASK,XWORK,YWORK,MAP(9),MAP(11))
      CALL PLTCP2(JN,MASK,XWORK,YWORK,MAP(11),MAP(13))
      CALL PLTCP2(JN,MASK,XWORK,YWORK,MAP(13),MAP(7))
      CALL PLTPTM(JN,MASK,XWORK,YWORK)
      GO TO 2360

 2370 CONTINUE
      RETURN

      END
