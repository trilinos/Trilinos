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

C $Id: pltdv3.f,v 1.4 2000/10/25 18:55:02 gdsjaar Exp $ 
C $Log: pltdv3.f,v $
C Revision 1.4  2000/10/25 18:55:02  gdsjaar
C In the pltli? functions, check for N==0 before doing any array
C accesses.
C
C Also changed all references to 'mask' to be arrays where they were
C scalars since downstream code seems to treat them as arrays.
C
C Revision 1.3  1993/07/19 17:06:39  gdsjaar
C Changed hex constants back to preceding X, --needed on cray. Works
C either way on other systems.
C
c Revision 1.2  1993/07/16  17:33:16  gdsjaar
c Integer constant too big on sun, replaced it with hexadecimal notation
c
c Revision 1.1  1993/07/16  16:48:05  gdsjaar
c Changed plt to library rather than single source file.
c 
C=======================================================================
      SUBROUTINE PLTDV3(MAP,N,UX,UY,UZ,VX,VY,VZ)
      REAL MAP(*),UX(*),UY(*),UZ(*),VX(*),VY(*),VZ(*)
      DIMENSION TUX(32),TUY(32),TUZ(32),TVX(32),TVY(32),TVZ(32),
     *          TTUX(32),TTUY(32),TTUZ(32),TTVX(32),TTVY(32),TTVZ(32),
     *          V1(3),Q1(3),V2(3),Q2(3)
      INTEGER MASK(1)
      INTEGER IZBIT(32)
      DATA IZBIT/1,          2,          4,          8,
     *          16,         32,         64,        128,
     *         256,        512,       1024,       2048,
     *        4096,       8192,      16384,      32768,
     *       65536,     131072,     262144,     524288,
     *     1048576,    2097152,    4194304,    8388608,
     *    16777216,   33554432,   67108864,  134217728,
     *   268435456,  536870912, 1073741824, X'80000000'/

      DO 2510 L = 1,3
         V1(L) = MAP(18+L-1) + MAP(15)*MAP(27+L-1)
         Q1(L) = MAP(27+L-1)
         V2(L) = MAP(18+L-1) + MAP(16)*MAP(27+L-1)
         Q2(L) = -MAP(27+L-1)
 2510 CONTINUE
      J = 0
 2530 IF (.NOT. (J.LT.N)) GO TO 2540
      JN = MIN(N-J,32)
      J1 = J + 1
      J = J + JN
      MASK(1) = -1
      CALL PLTCV3(JN,MASK,UX(J1),UY(J1),UZ(J1),VX(J1),VY(J1),VZ(J1),TUX,
     *            TUY,TUZ,TVX,TVY,TVZ,V1,Q1)
      CALL PLTCV3(JN,MASK,TUX,TUY,TUZ,TVX,TVY,TVZ,TTUX,TTUY,TTUZ,TTVX,
     *            TTVY,TTVZ,V2,Q2)
      IF (MAP(17).EQ.1.) THEN
         DO 2550 K = 1,JN
            JB = IZBIT(K)
            IF (IAND(JB,MASK(1)).NE.0) THEN
               PMS = (TTUX(K)-MAP(18))*MAP(27) +
     *               (TTUY(K)-MAP(19))*MAP(28) +
     *               (TTUZ(K)-MAP(20))*MAP(29)
               R = MAP(30)/PMS
               TUX(K) = R* ((TTUX(K)-MAP(18))*MAP(21)+
     *                  (TTUY(K)-MAP(19))*MAP(22)+
     *                  (TTUZ(K)-MAP(20))*MAP(23))
               TUY(K) = R* ((TTUX(K)-MAP(18))*MAP(24)+
     *                  (TTUY(K)-MAP(19))*MAP(25)+
     *                  (TTUZ(K)-MAP(20))*MAP(26))
               PMS = (TTVX(K)-MAP(18))*MAP(27) +
     *               (TTVY(K)-MAP(19))*MAP(28) +
     *               (TTVZ(K)-MAP(20))*MAP(29)
               R = MAP(30)/PMS
               TVX(K) = R* ((TTVX(K)-MAP(18))*MAP(21)+
     *                  (TTVY(K)-MAP(19))*MAP(22)+
     *                  (TTVZ(K)-MAP(20))*MAP(23))
               TVY(K) = R* ((TTVX(K)-MAP(18))*MAP(24)+
     *                  (TTVY(K)-MAP(19))*MAP(25)+
     *                  (TTVZ(K)-MAP(20))*MAP(26))
            END IF

 2550    CONTINUE

      ELSE IF (MAP(17).EQ.-1.) THEN
      END IF

      CALL PLTMV2(MAP,JN,MASK,TUX,TUY,TVX,TVY,TTUX,TTUY,TTVX,TTVY)
      CALL PLTVCM(JN,MASK,TTUX,TTUY,TTVX,TTVY)
      GO TO 2530

 2540 CONTINUE
      RETURN

      END
