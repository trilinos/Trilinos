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

C $Id: pltmp3.f,v 1.3 1993/07/19 17:06:40 gdsjaar Exp $ 
C $Log: pltmp3.f,v $
C Revision 1.3  1993/07/19 17:06:40  gdsjaar
C Changed hex constants back to preceding X, --needed on cray. Works
C either way on other systems.
C
c Revision 1.2  1993/07/16  17:33:17  gdsjaar
c Integer constant too big on sun, replaced it with hexadecimal notation
c
c Revision 1.1  1993/07/16  16:48:56  gdsjaar
c Changed plt to library rather than single source file.
c 
C=======================================================================
      SUBROUTINE PLTMP3(UMAP,N,MASK,PX,PY,PZ,QX,QY)
      DIMENSION UMAP(*),MASK(*),PX(*),PY(*),PZ(*),QX(*),QY(*)
      DIMENSION Q1(3),V1(3),Q2(3),V2(3),TPX(32),TPY(32)
      INTEGER IZBIT(32)
      DATA IZBIT/1,          2,          4,          8,
     *          16,         32,         64,        128,
     *         256,        512,       1024,       2048,
     *        4096,       8192,      16384,      32768,
     *       65536,     131072,     262144,     524288,
     *     1048576,    2097152,    4194304,    8388608,
     *    16777216,   33554432,   67108864,  134217728,
     *   268435456,  536870912, 1073741824, X'80000000'/

      DO 2060 L = 1,3
         V1(L) = UMAP(18+L-1) + UMAP(15)*UMAP(27+L-1)
         Q1(L) = UMAP(L+27-1)
         V2(L) = UMAP(18+L-1) + UMAP(16)*UMAP(27+L-1)
         Q2(L) = -UMAP(L+27-1)
 2060 CONTINUE
      J = 0
      KM = 0
 2080 IF (.NOT. (J.LT.N)) GO TO 2090
      JN = MIN(N-J,32)
      J1 = J + 1
      J = J + JN
      KM = KM + 1
      CALL PLTCP3(JN,MASK(KM),PX(J1),PY(J1),PZ(J1),V1,Q1)
      CALL PLTCP3(JN,MASK(KM),PX(J1),PY(J1),PZ(J1),V2,Q2)
      IF (UMAP(17).EQ.1.) THEN
         M = MASK(KM)
         DO 2100 K = 1,JN
            JB = IZBIT(K)
            IF (IAND(JB,M).NE.0) THEN
               PMS = (PX(K+J1-1)-UMAP(18))*UMAP(27) +
     *               (PY(K+J1-1)-UMAP(19))*UMAP(28) +
     *               (PZ(K+J1-1)-UMAP(20))*UMAP(29)
               R = UMAP(30)/PMS
               TPX(K) = R* ((PX(K+J1-1)-UMAP(18))*UMAP(21)+
     *                  (PY(K+J1-1)-UMAP(19))*UMAP(22)+
     *                  (PZ(K+J1-1)-UMAP(20))*UMAP(23))
               TPY(K) = R* ((PX(K+J1-1)-UMAP(18))*UMAP(24)+
     *                  (PY(K+J1-1)-UMAP(19))*UMAP(25)+
     *                  (PZ(K+J1-1)-UMAP(20))*UMAP(26))
            END IF

 2100    CONTINUE
      END IF

      IF (UMAP(17).EQ.-1.) THEN
      END IF

      CALL PLTMP2(UMAP(1),JN,MASK(KM),TPX,TPY,QX(J1),QY(J1))
      GO TO 2080

 2090 CONTINUE
      RETURN

      END
