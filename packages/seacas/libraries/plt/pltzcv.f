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

C $Id: pltzcv.f,v 1.3 1993/07/19 17:06:48 gdsjaar Exp $
C $Log: pltzcv.f,v $
C Revision 1.3  1993/07/19 17:06:48  gdsjaar
C Changed hex constants back to preceding X, --needed on cray. Works
C either way on other systems.
C
c Revision 1.2  1993/07/16  17:33:26  gdsjaar
c Integer constant too big on sun, replaced it with hexadecimal notation
c
c Revision 1.1  1993/07/16  16:50:04  gdsjaar
c Changed plt to library rather than single source file.
c
C=======================================================================
      SUBROUTINE PLTZCV(ZNEAR,ZFAR,N,MASK,PX,PY,PZ,QX,QY,QZ)
      INTEGER N
      INTEGER MASK(*)
      REAL PX(*)
      REAL PY(*)
      REAL PZ(*)
      REAL QX(*)
      REAL QY(*)
      REAL QZ(*)
      INTEGER IZBIT(32)
      DATA IZBIT/1,          2,          4,          8,
     *          16,         32,         64,        128,
     *         256,        512,       1024,       2048,
     *        4096,       8192,      16384,      32768,
     *       65536,     131072,     262144,     524288,
     *     1048576,    2097152,    4194304,    8388608,
     *    16777216,   33554432,   67108864,  134217728,
     *   268435456,  536870912, 1073741824, X'80000000'/

      J = 0
      KM = 0
      DZ = ZFAR - ZNEAR
 2420 IF (.NOT. (J.LT.N)) GO TO 2430
      JN = MIN(N-J,32)
      J1 = J
      J = J + JN
      KM = KM + 1
      M = MASK(KM)
      IF (M.EQ.0) THEN
         GO TO 2420

      END IF

      DO 2440 K = 1,JN
         JB = IZBIT(K)
         IF (IAND(M,JB).EQ.0) THEN
            GO TO 2440

         END IF

         X1 = PX(K+J1)
         Y1 = PY(K+J1)
         Z1 = PZ(K+J1)
         X2 = QX(K+J1)
         Y2 = QY(K+J1)
         Z2 = QZ(K+J1)
         FP = Z1 - ZNEAR
         FQ = Z2 - ZNEAR
         IF (FP.LT.0. .AND. FQ.LT.0.) THEN
            M = IAND(M,NOT(JB))
            GO TO 2440

         END IF

         IF (FP.GT.DZ .AND. FQ.GT.DZ) THEN
            M = IAND(M,NOT(JB))
            GO TO 2440

         END IF

         DF = FQ - FP
         IF (DF.GT.0.) THEN
            TN = (X2-X1)/DF
            SN = (Y2-Y1)/DF
            IF (FP.LT.0.) THEN
               Z1 = ZNEAR
               X1 = X1 - FP*TN
               Y1 = Y1 - FP*SN
            END IF

            IF (FQ.GT.DZ) THEN
               Z2 = ZFAR
               X2 = X2 + (DZ-FQ)*TN
               Y2 = Y2 + (DZ-FQ)*SN
            END IF

         ELSE IF (DF.LT.0.) THEN
            TN = (X2-X1)/DF
            SN = (Y2-Y1)/DF
            IF (FQ.LT.0.) THEN
               Z2 = ZNEAR
               X2 = X2 - FQ*TN
               Y2 = Y2 - FQ*SN
            END IF

            IF (FP.GT.DZ) THEN
               Z1 = ZFAR
               X1 = X1 + (DZ-FP)*TN
               Y1 = Y1 + (DZ-FP)*SN
            END IF

         END IF

         PX(K+J1) = X1
         PY(K+J1) = Y1
         PZ(K+J1) = Z1
         QX(K+J1) = X2
         QY(K+J1) = Y2
         QZ(K+J1) = Z2
 2440 CONTINUE
      MASK(KM) = M
      GO TO 2420

 2430 CONTINUE
      RETURN

      END
