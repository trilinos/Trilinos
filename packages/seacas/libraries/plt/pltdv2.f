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

C $Id: pltdv2.f,v 1.2 2000/10/25 18:55:02 gdsjaar Exp $
C $Log: pltdv2.f,v $
C Revision 1.2  2000/10/25 18:55:02  gdsjaar
C In the pltli? functions, check for N==0 before doing any array
C accesses.
C
C Also changed all references to 'mask' to be arrays where they were
C scalars since downstream code seems to treat them as arrays.
C
C Revision 1.1  1993/07/16 16:48:04  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      SUBROUTINE PLTDV2(MAP,N,PX,PY,QX,QY)
      REAL MAP(*),PX(*),PY(*),QX(*),QY(*)
      DIMENSION PPX(32),PPY(32),QQX(32),QQY(32)
      INTEGER MASK(1)

      J = 0
 2470 IF (J.LT.N) THEN
         JN = MIN(N-J,32)
         J1 = J
         J = J + JN
         DO 2490 I = 1,JN
            PPX(I) = MAP(1)*PX(I+J1) + MAP(3)*PY(I+J1) + MAP(5)
            QQX(I) = MAP(1)*QX(I+J1) + MAP(3)*QY(I+J1) + MAP(5)
            PPY(I) = MAP(2)*PX(I+J1) + MAP(4)*PY(I+J1) + MAP(6)
            QQY(I) = MAP(2)*QX(I+J1) + MAP(4)*QY(I+J1) + MAP(6)
 2490    CONTINUE

         MASK(1) = -1
         CALL PLTVWV(MAP(7),MAP(11),JN,MASK,PPX,PPY,QQX,QQY)
         CALL PLTVCM(JN,MASK,PPX,PPY,QQX,QQY)
         GO TO 2470

      END IF
      RETURN

      END
