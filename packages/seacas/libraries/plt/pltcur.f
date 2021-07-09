C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

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
      include 'izbit.inc'
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
