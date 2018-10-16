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

C $Id: pltcur.f,v 1.5 2000/10/25 18:55:02 gdsjaar Exp $ 
C $Log: pltcur.f,v $
C Revision 1.5  2000/10/25 18:55:02  gdsjaar
C In the pltli? functions, check for N==0 before doing any array
C accesses.
C
C Also changed all references to 'mask' to be arrays where they were
C scalars since downstream code seems to treat them as arrays.
C
C Revision 1.4  2000/10/25 13:32:35  gdsjaar
C Modified intrinsic functions to use generic versions to avoid warnings on SGI 64-bit compiles
C
C Revision 1.3  1993/07/19 17:06:35  gdsjaar
C Changed hex constants back to preceding X, --needed on cray. Works
C either way on other systems.
C
c Revision 1.2  1993/07/16  17:33:12  gdsjaar
c Integer constant too big on sun, replaced it with hexadecimal notation
c
c Revision 1.1  1993/07/16  16:47:54  gdsjaar
c Changed plt to library rather than single source file.
c 
C=======================================================================
      SUBROUTINE PLTCUR(X,Y,NUM)
      REAL DEVCAP(23)
      REAL DEFOUT(7)
      COMMON /STATUS/DEVCAP,DEFOUT
      REAL DEVP(5)
      COMMON /DEVICE/DEVP
      REAL COLP(3)
      REAL PALETT(3,16)
      COMMON /COLOR/COLP,PALETT
      REAL TEXTP(40)
      COMMON /TEXT/TEXTP
      REAL VECTP(5)
      REAL XCUR
      REAL YCUR
      COMMON /VECTRC/VECTP,XCUR,YCUR
      INTEGER IDEX(200,2)
      INTEGER NVECT(200,2)
      REAL XSIZE(200,2)
      REAL YSIZE(200,2)
      REAL X0(2300,2)
      REAL Y0(2300,2)
      REAL X1(2300,2)
      REAL Y1(2300,2)
      COMMON /FONT/IDEX,NVECT,XSIZE,YSIZE,X0,Y0,X1,Y1
      REAL GRAPHP(100)
      COMMON /GRAPH/GRAPHP
      COMMON /MAPPAR/MAPP(11)
      REAL MAPP
      COMMON /STORAG/MEMORY(1000)
      LOGICAL CPUIFC
      DIMENSION X(*),Y(*),X0T(32),Y0T(32),X1T(32),Y1T(32)
      INTEGER MASKS(1)
      INTEGER IZBIT(32)
      DATA IZBIT/1,          2,          4,          8,
     *          16,         32,         64,        128,
     *         256,        512,       1024,       2048,
     *        4096,       8192,      16384,      32768,
     *       65536,     131072,     262144,     524288,
     *     1048576,    2097152,    4194304,    8388608,
     *    16777216,   33554432,   67108864,  134217728,
     *   268435456,  536870912, 1073741824, X'80000000'/
      DATA SYMRAT/120./

      IF (NUM.LE.0) THEN
         RETURN

      END IF

      CALL PLTSVV
      CALL PLTSVD
      CALL PLTSVT
      CALL PLTSTV(1,GRAPHP(5))
      CALL PLTSTD(1,GRAPHP(38))
      CALL PLTSTV(2,GRAPHP(63))
      IF (GRAPHP(6).EQ.1.) THEN
         CALL PLTSTT(3,0.)
         CALL PLTSTT(4,0.)
         SYMSIZ = (GRAPHP(3)+GRAPHP(4))*GRAPHP(46)/ (10.*SYMRAT)
         CALL PLTSTT(2,SYMSIZ)
         CALL PLTSTT(11,GRAPHP(66))
      END IF

      J = 0
      NP = 0
      IF (NUM.EQ.1) THEN
         IF (GRAPHP(6).EQ.1. .AND. SYMSIZ.GT.0.0) THEN
            IF (GRAPHP(21).EQ.1.) THEN
               X0T(1) = X(1)
               Y0T(1) = Y(1)

            ELSE IF (GRAPHP(21).EQ.2.) THEN
               IF (X(1).LE.0.) THEN
                  RETURN

               END IF

               X0T(1) = LOG10(X(1))
               Y0T(1) = Y(1)

            ELSE IF (GRAPHP(21).EQ.3.) THEN
               X0T(1) = X(1)
               IF (Y(1).LE.0.) THEN
                  RETURN

               END IF

               Y0T(1) = LOG10(Y(1))

            ELSE IF (GRAPHP(21).EQ.4.) THEN
               IF (X(1).LE.0.) THEN
                  RETURN

               END IF

               IF (Y(1).LE.0.) THEN
                  RETURN

               END IF

               X0T(1) = LOG10(X(1))
               Y0T(1) = LOG10(Y(1))
            END IF

            MASKS(1) = -1
            CALL PLTMP2(GRAPHP(7),1,MASKS,X0T,Y0T,X1T,Y1T)
            JB = IZBIT(1)
            IF (IAND(MASKS(1),JB).NE.0) THEN
               CALL PLTSTD(1,GRAPHP(75))
               CALL PLTXTS(X1T(1),Y1T(1),CHAR(INT(GRAPHP(47))))
            END IF

         END IF

      ELSE
         IF (GRAPHP(21).EQ.1.) THEN
            XSAV = X(1)
            YSAV = Y(1)

         ELSE IF (GRAPHP(21).EQ.2.) THEN
            DO 2000 I = 1,NUM
               IF (X(I).LE.0.) THEN
                  GO TO 2000

               END IF

               XSAV = X(I)
               YSAV = Y(I)
               GO TO 2010

 2000       CONTINUE
 2010       CONTINUE

         ELSE IF (GRAPHP(21).EQ.3.) THEN
            DO 2020 I = 1,NUM
               IF (Y(I).LE.0.) THEN
                  GO TO 2020

               END IF

               XSAV = X(I)
               YSAV = Y(I)
               GO TO 2030

 2020       CONTINUE
 2030       CONTINUE

         ELSE IF (GRAPHP(21).EQ.4.) THEN
            DO 2040 I = 1,NUM
               IF (X(I).LE.0. .OR. Y(I).LE.0.) THEN
                  GO TO 2040

               END IF

               XSAV = X(I)
               YSAV = Y(I)
               GO TO 2050

 2040       CONTINUE
 2050       CONTINUE
         END IF

 2060    IF (.NOT. (J.LT.NUM-1)) GO TO 2070
         K = MIN(NUM-J,32)
         IF (GRAPHP(21).EQ.1.) THEN
            NV = 0
            DO 2080 I = 1,K - 1
               NV = NV + 1
               X0T(I) = XSAV
               Y0T(I) = YSAV
               X1T(I) = X(I+J+1)
               Y1T(I) = Y(I+J+1)
               XSAV = X(I+J+1)
               YSAV = Y(I+J+1)
 2080       CONTINUE

         ELSE IF (GRAPHP(21).EQ.2.) THEN
            NV = 0
            DO 2100 I = 1,K - 1
               IF (X(J+I+1).LE.0.) THEN
                  GO TO 2100

               END IF

               NV = NV + 1
               X0T(NV) = LOG10(XSAV)
               Y0T(NV) = YSAV
               X1T(NV) = LOG10(X(J+I+1))
               Y1T(NV) = Y(J+I+1)
               XSAV = X(I+J+1)
               YSAV = Y(I+J+1)
 2100       CONTINUE

         ELSE IF (GRAPHP(21).EQ.3.) THEN
            NV = 0
            DO 2120 I = 1,K - 1
               IF (Y(J+I+1).LE.0.) THEN
                  GO TO 2120

               END IF

               NV = NV + 1
               X0T(NV) = XSAV
               Y0T(NV) = LOG10(YSAV)
               X1T(NV) = X(J+I+1)
               Y1T(NV) = LOG10(Y(J+I+1))
               XSAV = X(I+J+1)
               YSAV = Y(I+J+1)
 2120       CONTINUE

         ELSE IF (GRAPHP(21).EQ.4.) THEN
            NV = 0
            DO 2140 I = 1,K - 1
               IF (X(J+I+1).LE.0.) THEN
                  GO TO 2140

               END IF

               IF (Y(J+I+1).LE.0.) THEN
                  GO TO 2140

               END IF

               NV = NV + 1
               X0T(NV) = LOG10(XSAV)
               Y0T(NV) = LOG10(YSAV)
               X1T(NV) = LOG10(X(J+I+1))
               Y1T(NV) = LOG10(Y(J+I+1))
               XSAV = X(I+J+1)
               YSAV = Y(I+J+1)
 2140       CONTINUE
         END IF

         CALL PLTDV2(GRAPHP(7),NV,X0T,Y0T,X1T,Y1T)
         IF (GRAPHP(6).EQ.1. .AND. SYMSIZ.GT.0.0) THEN
            CALL PLTSTD(1,GRAPHP(75))
            MASKS(1) = -1
            XT = X1T(NV)
            YT = Y1T(NV)
            CALL PLTMP2(GRAPHP(7),NV,MASKS,X0T,Y0T,X1T,Y1T)
            DO 2160 L = 1,NV
               JB = IZBIT(L)
               IF (IAND(MASKS(1),JB).NE.0 .AND.
     *             MOD(L+NP+INT(GRAPHP(23))-1,INT(GRAPHP(23))).EQ.
     *             0) THEN
                  CALL PLTXTS(X1T(L),Y1T(L),CHAR(INT(GRAPHP(47))))
               END IF

 2160       CONTINUE
            CALL PLTSTD(1,GRAPHP(38))
         END IF

         NP = NP + NV
         J = J + K - 1
         IF (J+1.GE.NUM .AND. (GRAPHP(6).EQ.1..AND.SYMSIZ.GT.0.0)) THEN
            X0T(1) = XT
            Y0T(1) = YT
            CALL PLTSTD(1,GRAPHP(75))
            MASKS(1) = -1
            CALL PLTMP2(GRAPHP(7),1,MASKS,X0T,Y0T,X1T,Y1T)
            JB = IZBIT(1)
            IF (IAND(MASKS(1),JB).NE.0 .AND.
     *          MOD(NP+INT(GRAPHP(23)),INT(GRAPHP(23))).EQ.0) THEN
               CALL PLTXTS(X1T(1),Y1T(1),CHAR(INT(GRAPHP(47))))
            END IF

         END IF

         IF (CPUIFC(.FALSE.)) THEN
            GO TO 2070

         END IF

         GO TO 2060

 2070    CONTINUE
      END IF

      CALL PLTRET
      CALL PLTRED
      CALL PLTREV
      RETURN

      END
