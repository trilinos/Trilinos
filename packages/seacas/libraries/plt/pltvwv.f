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

C $Id: pltvwv.f,v 1.4 1993/07/19 17:06:46 gdsjaar Exp $
C $Log: pltvwv.f,v $
C Revision 1.4  1993/07/19 17:06:46  gdsjaar
C Changed hex constants back to preceding X, --needed on cray. Works
C either way on other systems.
C
c Revision 1.3  1993/07/19  14:38:19  gdsjaar
c Reformatted flow of control
c
c Revision 1.2  1993/07/16  17:33:24  gdsjaar
c Integer constant too big on sun, replaced it with hexadecimal notation
c
c Revision 1.1  1993/07/16  16:49:51  gdsjaar
c Changed plt to library rather than single source file.
c
C=======================================================================
      SUBROUTINE PLTVWV(PLL,PUR,N,MASK,PX,PY,QX,QY)
      REAL PLL(2)
      REAL PUR(2)
      INTEGER N
      INTEGER MASK(*)
      REAL PX(*)
      REAL PY(*)
      REAL QX(*)
      REAL QY(*)
      INTEGER IZBIT(32)
      DATA IZBIT/1,          2,          4,          8,
     *          16,         32,         64,        128,
     *         256,        512,       1024,       2048,
     *        4096,       8192,      16384,      32768,
     *       65536,     131072,     262144,     524288,
     *     1048576,    2097152,    4194304,    8388608,
     *    16777216,   33554432,   67108864,  134217728,
     *   268435456,  536870912, 1073741824, X'80000000'/

      PUR1 = PUR(1) + .0001
      PUR2 = PUR(2) + .0001
      PLL1 = PLL(1) - .0001
      PLL2 = PLL(2) - .0001
      DX = PUR1 - PLL1
      DY = PUR2 - PLL2
      J = 0
      KM = 0
   10 CONTINUE
      IF ((J.LT.N)) THEN
         JN = MIN(N-J,32)
         J1 = J
         J = J + JN
         KM = KM + 1
         M = MASK(KM)
         IF (M.NE.0) THEN

            DO 20 K = 1,JN
               JB = IZBIT(K)
               IF (IAND(M,JB).NE.0) THEN

                  X1 = PX(K+J1)
                  Y1 = PY(K+J1)
                  X2 = QX(K+J1)
                  Y2 = QY(K+J1)
                  FP = X1 - PLL1
                  FQ = X2 - PLL1
                  IF (FP.LT.0. .AND. FQ.LT.0.) THEN
                     M = IAND(M,NOT(JB))

                  ELSE IF (FP.GT.DX .AND. FQ.GT.DX) THEN
                     M = IAND(M,NOT(JB))

                  ELSE

                     DF = FQ - FP
                     IF (DF.GT.0.) THEN
                        TN = (Y2-Y1)/DF
                        IF (FP.LT.0.) THEN
                           X1 = PLL1
                           Y1 = Y1 - FP*TN
                        END IF

                        IF (FQ.GT.DX) THEN
                           X2 = PUR1
                           Y2 = Y2 + (DX-FQ)*TN
                        END IF

                     ELSE IF (DF.LT.0.) THEN
                        TN = (Y2-Y1)/DF
                        IF (FQ.LT.0.) THEN
                           X2 = PLL1
                           Y2 = Y2 - FQ*TN
                        END IF

                        IF (FP.GT.DX) THEN
                           X1 = PUR1
                           Y1 = Y1 + (DX-FP)*TN
                        END IF

                     END IF

                     FP = Y1 - PLL2
                     FQ = Y2 - PLL2
                     IF (FP.LT.0. .AND. FQ.LT.0.) THEN
                        M = IAND(M,NOT(JB))

                     ELSE IF (FP.GT.DY .AND. FQ.GT.DY) THEN
                        M = IAND(M,NOT(JB))

                     ELSE

                        DF = FQ - FP
                        IF (DF.GT.0.) THEN
                           TN = (X2-X1)/DF
                           IF (FP.LT.0.) THEN
                              Y1 = PLL2
                              X1 = X1 - FP*TN
                           END IF

                           IF (FQ.GT.DY) THEN
                              Y2 = PUR2
                              X2 = X2 + (DY-FQ)*TN
                           END IF

                        ELSE IF (DF.LT.0.) THEN
                           TN = (X2-X1)/DF
                           IF (FQ.LT.0.) THEN
                              Y2 = PLL2
                              X2 = X2 - FQ*TN
                           END IF

                           IF (FP.GT.DY) THEN
                              Y1 = PUR2
                              X1 = X1 + (DY-FP)*TN
                           END IF

                        END IF

                        PX(K+J1) = X1
                        PY(K+J1) = Y1
                        QX(K+J1) = X2
                        QY(K+J1) = Y2
                     END IF

                  END IF

               END IF

   20       CONTINUE
            MASK(KM) = M
         END IF

         GO TO 10

      END IF

      END
