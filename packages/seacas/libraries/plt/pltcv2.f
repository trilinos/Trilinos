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

C $Id: pltcv2.f,v 1.3 1993/07/19 17:06:36 gdsjaar Exp $ 
C $Log: pltcv2.f,v $
C Revision 1.3  1993/07/19 17:06:36  gdsjaar
C Changed hex constants back to preceding X, --needed on cray. Works
C either way on other systems.
C
c Revision 1.2  1993/07/16  17:33:13  gdsjaar
c Integer constant too big on sun, replaced it with hexadecimal notation
c
c Revision 1.1  1993/07/16  16:47:55  gdsjaar
c Changed plt to library rather than single source file.
c 
C=======================================================================
      SUBROUTINE PLTCV2(N,MASK,PX,PY,QX,QY,PPX,PPY,QQX,QQY,C1,C2)
      DIMENSION MASK(*),PX(*),PY(*),QX(*),QY(*),PPX(*),PPY(*),QQX(*),
     *          QQY(*),C1(*),C2(*)
      INTEGER IZBIT(32)
      DATA IZBIT/1,          2,          4,          8,
     *          16,         32,         64,        128,
     *         256,        512,       1024,       2048,
     *        4096,       8192,      16384,      32768,
     *       65536,     131072,     262144,     524288,
     *     1048576,    2097152,    4194304,    8388608,
     *    16777216,   33554432,   67108864,  134217728,
     *   268435456,  536870912, 1073741824, X'80000000'/

      CX = C1(1)
      CY = C1(2)
      DX = C2(1) - CX
      DY = C2(2) - CY
      J = 0
      KM = 0
 2100 IF (.NOT. (J.LT.N)) GO TO 2110
      JN = MIN(N-J,32)
      J1 = J
      KM = KM + 1
      J = J + JN
      M = MASK(KM)
      IF (M.EQ.0) THEN
         GO TO 2100

      END IF

      DO 2120 K = 1,JN
         JB = IZBIT(K)
         IF (IAND(JB,M).EQ.0) THEN
            GO TO 2120

         END IF

         X1 = PX(J1+K)
         Y1 = PY(J1+K)
         X2 = QX(J1+K)
         Y2 = QY(J1+K)
         FP = (Y1-CY)*DX - (X1-CX)*DY
         FQ = (Y2-CY)*DX - (X2-CX)*DY
         IF (FP.LT.0. .AND. FQ.LT.0.) THEN
            M = IAND(M,NOT(JB))

         ELSE IF (FP.LT.0.) THEN
            XL = FQ/ (FQ-FP)
            X1 = X2 + XL* (X1-X2)
            Y1 = Y2 + XL* (Y1-Y2)

         ELSE IF (FQ.LT.0.) THEN
            XL = FP/ (FP-FQ)
            X2 = X1 + XL* (X2-X1)
            Y2 = Y1 + XL* (Y2-Y1)
         END IF

         PPX(K+J1) = X1
         PPY(K+J1) = Y1
         QQX(K+J1) = X2
         QQY(K+J1) = Y2
 2120 CONTINUE
      MASK(KM) = M
      GO TO 2100

 2110 CONTINUE
      RETURN

      END
