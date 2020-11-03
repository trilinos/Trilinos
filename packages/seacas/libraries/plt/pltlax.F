C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTLAX(X,Y,XLENG,YLENG,TYPE,MINEXP,MAXEXP,LABEL,UNITS)
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
      REAL LENMAJ,LENMIN,MAJRAT,MINRAT,NUMRAT,NUMSIZ,LABSIZ,NUMSZM,
     *     LABRAT,MAJTIC,LONNUM,LABSCA,LOGTAB(8),MINEXP,MAXEXP
      CHARACTER*(*) LABEL,UNITS,TYPE
      CHARACTER TTYPE*1,LINE1*10,LINE*10,TLABEL*132,CINPUT*132
      INTEGER PLTITL
      LOGICAL FIRST,PLTNER
      DATA MAJRAT/50./
      DATA MINRAT/2./
      DATA NUMRAT/50./
      DATA LABRAT/40./
      DATA FIRST/.TRUE./
      DATA FUDGE/.0025/

      IF (MAXEXP.LE.MINEXP) THEN
         CALL PLTFLU
         CALL SIORPT('PLTLAX',
     *         'Maximum exponent must be greater than minimum exponent.'
     *               ,2)
         RETURN

      END IF

      IF (FIRST) THEN
         DO 2530 I = 1,8
            LOGTAB(I) = LOG10(DBLE(I+1))
 2530    CONTINUE
         FIRST = .FALSE.
      END IF

      LOWEXP = PLTITL(MINEXP)
      NUM = LOWEXP
      TTYPE = TYPE
      CALL CHRUP(TTYPE,TTYPE)
      CALL PLTSVT
      CALL PLTSVD
      CALL PLTSVV
      CALL PLTSTT(3,0.)
      CALL PLTSTT(4,0.)
      CALL PLTSTV(1,1.)
#if NeedsDoubleEscape
      LINE(1:4) = '10\\^'
#else
      LINE(1:4) = '10\^'
#endif
 2550 CONTINUE
      IF (TTYPE.EQ.'X') THEN
         LABSCA = (XLENG+YLENG)/2.
         NUMSIZ = (LABSCA* (GRAPHP(44)/5.))/NUMRAT
         NUMSZM = NUMSIZ*.8
         LABSIZ = (LABSCA* (GRAPHP(45)/5.))/LABRAT
         LENMAJ = XLENG/MAJRAT
         LENMIN = LENMAJ/MINRAT
         LFLAG = 0
         MAJTIC = XLENG/ (MAXEXP-MINEXP)
         FMJTIC = X - (MINEXP-LOWEXP)*MAJTIC
         XMAJ = FMJTIC
         YOFF = Y - NUMSIZ*1.8
         YOFFM = YOFF + (NUMSIZ-NUMSZM)/2.
         CALL PLTSTV(2,GRAPHP(62))
         CALL PLTSTD(1,GRAPHP(37))
         CALL PLTVCT(1,X,Y,X+XLENG,Y)
         IF (GRAPHP(32).EQ.1.) THEN
            CALL PLTVCT(1,X,Y+YLENG,X+XLENG,Y+YLENG)
         END IF

         IF (PLTNER(FMJTIC,X) .AND. NUMSIZ.GT.0.0) THEN
            CALL PLTSTT(2,NUMSIZ)
            CALL PLTSTT(11,GRAPHP(64))
            CALL PLTSTD(1,GRAPHP(76))
            IF (GRAPHP(40).EQ.1. .OR. GRAPHP(40).EQ.3.) THEN
               CALL PLTLOA(LINE1,NUM,0)
               LINE(5:) = LINE1
               CALL PLTXSL(LINE,TLEN)
               CALL PLTXTS(X-TLEN/2.,YOFF,LINE)
               LFLAG = 1

            ELSE
               CALL PLTLOA(LINE1,NUM,1)
               CALL PLTXSL(LINE1,TLEN)
               CALL PLTXTS(X-TLEN/2.,YOFF,LINE1)
            END IF

         END IF

 2580    CONTINUE
         IF (CPUIFC(.FALSE.)) THEN
            GO TO 2600

         END IF

         DO 2610 I = 1,8
            XNEW = XMAJ + LOGTAB(I)*MAJTIC
            IF (PLTNER(XNEW,X) .AND. NUMSZM.GT.0.0 .AND.
     *          GRAPHP(48).EQ.1.) THEN
               CALL PLTSTT(11,GRAPHP(64))
               CALL PLTSTD(1,GRAPHP(76))
               IF (GRAPHP(40).EQ.1. .OR. GRAPHP(40).EQ.3.) THEN
                  CALL PLTSTT(2,NUMSZM)
                  CALL PLTLOA(LINE1,I+1,0)
                  CALL PLTXSL(LINE1,TLEN)
                  CALL PLTXTS(XNEW-TLEN/2.,YOFFM,LINE1)

               ELSE
                  CALL PLTSTT(2,NUMSIZ)
                  CALL PLTLOD(LINE1,I+1,NUM)
                  CALL PLTXSL(LINE1,TLEN)
                  CALL PLTXTS(XNEW-TLEN/2.,YOFF,LINE1)
               END IF

            END IF

            IF (XNEW-FUDGE.LE.X) THEN
               GO TO 2610

            END IF

            IF (XNEW.GT.X+XLENG-FUDGE) THEN
               GO TO 2620

            END IF

            CALL PLTSTD(1,GRAPHP(77))
            CALL PLTSTV(2,GRAPHP(67))
            CALL PLTVCT(1,XNEW,Y,XNEW,Y+LENMIN)
            IF (GRAPHP(32).EQ.1.) THEN
               CALL PLTVCT(1,XNEW,Y+YLENG,XNEW,Y+YLENG-LENMIN)
            END IF

            IF (GRAPHP(73).NE.0.) THEN
               CALL PLTSTV(1,GRAPHP(73))
               CALL PLTSTD(1,GRAPHP(74))
               CALL PLTSTV(2,GRAPHP(69))
               IF (GRAPHP(32).EQ.1.) THEN
                  CALL PLTVCT(1,XNEW,Y+LENMIN,XNEW,Y+YLENG-LENMIN)

               ELSE
                  CALL PLTVCT(1,XNEW,Y+LENMIN,XNEW,Y+YLENG)
               END IF

               CALL PLTSTV(1,1.)
            END IF

            IF (GRAPHP(48).EQ.1. .AND. NUMSIZ.GT.0.0) THEN
               CALL PLTSTT(11,GRAPHP(64))
               CALL PLTSTD(1,GRAPHP(76))
               IF (GRAPHP(40).EQ.1. .OR. GRAPHP(40).EQ.3.) THEN
                  CALL PLTSTT(2,NUMSZM)
                  CALL PLTLOA(LINE1,I+1,0)
                  CALL PLTXSL(LINE1,TLEN)
                  CALL PLTXTS(XNEW-TLEN/2.,YOFFM,LINE1)

               ELSE
                  CALL PLTSTT(2,NUMSIZ)
                  CALL PLTLOD(LINE1,I+1,NUM)
                  CALL PLTXSL(LINE1,TLEN)
                  CALL PLTXTS(XNEW-TLEN/2.,YOFF,LINE1)
               END IF

            END IF

 2610    CONTINUE
 2620    CONTINUE
         XMAJ = XMAJ + MAJTIC
         NUM = NUM + 1
         IF (XMAJ.GT.X+XLENG+FUDGE .AND. GRAPHP(32).EQ.0.) THEN
            GO TO 2600

         ELSE IF (XMAJ.GT.X+XLENG-FUDGE .AND. GRAPHP(32).EQ.1.) THEN
            GO TO 2600

         END IF

         CALL PLTSTD(1,GRAPHP(77))
         CALL PLTSTV(2,GRAPHP(67))
         CALL PLTVCT(1,XMAJ,Y,XMAJ,Y+LENMAJ)
         IF (GRAPHP(32).EQ.1.) THEN
            CALL PLTVCT(1,XMAJ,Y+YLENG,XMAJ,Y+YLENG-LENMAJ)
         END IF

         IF (GRAPHP(35).NE.0. .OR. GRAPHP(73).NE.0.) THEN
            IF (GRAPHP(35).EQ.0.) THEN
               CALL PLTSTV(1,GRAPHP(73))

            ELSE
               CALL PLTSTV(1,GRAPHP(35))
            END IF

            IF (GRAPHP(35).EQ.0.) THEN
               CALL PLTSTD(1,GRAPHP(74))

            ELSE
               CALL PLTSTD(1,GRAPHP(36))
            END IF

            IF (GRAPHP(35).EQ.0.) THEN
               CALL PLTSTV(2,GRAPHP(69))

            ELSE
               CALL PLTSTV(2,GRAPHP(68))
            END IF

            IF (GRAPHP(32).EQ.1.) THEN
               CALL PLTVCT(1,XMAJ,Y+LENMAJ,XMAJ,Y+YLENG-LENMAJ)

            ELSE
               CALL PLTVCT(1,XMAJ,Y+LENMAJ,XMAJ,Y+YLENG)
            END IF

            CALL PLTSTV(1,1.)
         END IF

         IF (NUMSIZ.GT.0.0) THEN
            CALL PLTSTT(2,NUMSIZ)
            CALL PLTSTT(11,GRAPHP(64))
            CALL PLTSTD(1,GRAPHP(76))
            IF (GRAPHP(40).EQ.1. .OR. GRAPHP(40).EQ.3.) THEN
               CALL PLTLOA(LINE1,NUM,0)
               LINE(5:) = LINE1
               CALL PLTXSL(LINE,TLEN)
               CALL PLTXTS(XMAJ-TLEN/2.,YOFF,LINE)
               LFLAG = 1

            ELSE
               CALL PLTLOA(LINE1,NUM,1)
               CALL PLTXSL(LINE1,TLEN)
               CALL PLTXTS(XMAJ-TLEN/2.,YOFF,LINE1)
            END IF

         END IF

         GO TO 2580

 2600    CONTINUE
         IF (PLTNER(XNEW,X+XLENG) .AND. NUMSZM.GT.0.0 .AND.
     *       GRAPHP(48).EQ.1.) THEN
            CALL PLTSTT(11,GRAPHP(64))
            CALL PLTSTD(1,GRAPHP(76))
            IF (GRAPHP(40).EQ.1. .OR. GRAPHP(40).EQ.3.) THEN
               CALL PLTSTT(2,NUMSZM)
               CALL PLTLOA(LINE1,I+1,0)
               CALL PLTXSL(LINE1,TLEN)
               CALL PLTXTS(XNEW-TLEN/2.,YOFFM,LINE1)

            ELSE
               CALL PLTSTT(2,NUMSIZ)
               CALL PLTLOD(LINE1,I+1,NUM-1)
               CALL PLTXSL(LINE1,TLEN)
               CALL PLTXTS(XNEW-TLEN/2.,YOFF,LINE1)
            END IF

         END IF

         IF (PLTNER(XMAJ,X+XLENG) .AND. NUMSIZ.GT.0.0) THEN
            CALL PLTSTT(2,NUMSIZ)
            CALL PLTSTT(11,GRAPHP(64))
            CALL PLTSTD(1,GRAPHP(76))
            IF (GRAPHP(40).EQ.1. .OR. GRAPHP(40).EQ.3.) THEN
               CALL PLTLOA(LINE1,NUM,0)
               LINE(5:) = LINE1
               CALL PLTXSL(LINE,TLEN)
               CALL PLTXTS(XMAJ-TLEN/2.,YOFF,LINE)
               LFLAG = 1

            ELSE
               CALL PLTLOA(LINE1,NUM,1)
               CALL PLTXSL(LINE1,TLEN)
               CALL PLTXTS(XMAJ-TLEN/2.,YOFF,LINE1)
            END IF

         END IF

         IF (CPUIFC(.FALSE.)) THEN
            GO TO 2570

         END IF

         IF (LABSIZ.GT.0.0) THEN
            TLABEL = ' '
            LL = 0
            IF (LABEL.NE.' ') THEN
               TLABEL = LABEL
               CALL CHRTRM(TLABEL,LL)
#if NeedsDoubleEscape
               TLABEL(LL+1:LL+2) = '\\-'
#else
               TLABEL(LL+1:LL+2) = '\-'
#endif
               LL = LL + 2
            END IF

            IF (UNITS.NE.' ') THEN
               TLABEL(LL+6:) = UNITS
               CALL CHRTRM(TLABEL,LL)
#if NeedsDoubleEscape
               TLABEL(LL+1:LL+2) = '\\-'
#else
               TLABEL(LL+1:LL+2) = '\-'
#endif
               LL = LL + 2
            END IF

            IF (LFLAG.EQ.0 .AND. (GRAPHP(40).EQ.1..OR.
     *          GRAPHP(40).EQ.3.) .AND. LOWEXP.NE.0 .AND.
     *          NUMSIZ.GT.0.) THEN
               CALL PLTLOA(LINE1,LOWEXP,0)
               CALL CHRTRM(LINE1,L)
#if NeedsDoubleEscape
               TLABEL(LL+1:) = ' (*10\\^'//LINE1(:L)//'\\-)'
#else
               TLABEL(LL+1:) = ' (*10\^'//LINE1(:L)//'\-)'
#endif
            END IF

            CINPUT = TLABEL
            CALL CHRSTR(CINPUT,TLABEL,LL)
            IF (LL.GT.0) THEN
               CALL PLTSTT(2,LABSIZ)
               CALL PLTSTT(11,GRAPHP(65))
               CALL PLTSTD(1,GRAPHP(39))
               CALL PLTXSL(TLABEL(1:LL),TLEN)
               XLAB = X + (XLENG-TLEN)/2.
               YLAB = YOFF - LABSIZ*2.0
               CALL PLTXTS(XLAB,YLAB,TLABEL(1:LL))
            END IF

            IF (CPUIFC(.FALSE.)) THEN
               GO TO 2570

            END IF

         END IF

      ELSE IF (TTYPE.EQ.'Y') THEN
         LABSCA = (XLENG+YLENG)/2.
         NUMSIZ = (LABSCA* (GRAPHP(88)/5.))/NUMRAT
         NUMSZM = NUMSIZ*.8
         LABSIZ = (LABSCA* (GRAPHP(89)/5.))/LABRAT
         LENMAJ = YLENG/MAJRAT
         LENMIN = LENMAJ/MINRAT
         LFLAG = 0
         MAJTIC = YLENG/ (MAXEXP-MINEXP)
         FMJTIC = Y - (MINEXP-LOWEXP)*MAJTIC
         YMAJ = FMJTIC
         XOFF = NUMSIZ*.8
         YOFF = NUMSIZ/2.
         YOFFM = NUMSZM/2.
         LONNUM = 0.
         CALL PLTSTV(2,GRAPHP(62))
         CALL PLTSTD(1,GRAPHP(37))
         CALL PLTVCT(1,X,Y,X,Y+YLENG)
         IF (GRAPHP(32).EQ.1.) THEN
            CALL PLTVCT(1,X+XLENG,Y,X+XLENG,Y+YLENG)
         END IF

         IF (PLTNER(FMJTIC,Y) .AND. NUMSIZ.GT.0.0) THEN
            CALL PLTSTT(2,NUMSIZ)
            CALL PLTSTT(11,GRAPHP(64))
            CALL PLTSTD(1,GRAPHP(76))
            IF (GRAPHP(40).EQ.2. .OR. GRAPHP(40).EQ.3.) THEN
               CALL PLTLOA(LINE1,NUM,0)
               LINE(5:) = LINE1
               CALL PLTXSL(LINE,TLEN)
               IF (GRAPHP(92).EQ.1.) THEN
                  LONNUM = NUMSIZ
                  CALL PLTSTT(3,90.)
                  CALL PLTXTS(X-XOFF,Y-TLEN/2.,LINE)
                  CALL PLTSTT(3,0.)

               ELSE
                  CALL PLTXTS(X- (TLEN+XOFF),Y-YOFF,LINE)
                  IF (TLEN.GT.LONNUM) THEN
                     LONNUM = TLEN
                  END IF

               END IF

               LFLAG = 1

            ELSE
               CALL PLTLOA(LINE1,NUM,1)
               CALL PLTXSL(LINE1,TLEN)
               IF (GRAPHP(92).EQ.1.) THEN
                  LONNUM = NUMSIZ
                  CALL PLTSTT(3,90.)
                  CALL PLTXTS(X-XOFF,Y-TLEN/2.,LINE1)
                  CALL PLTSTT(3,0.)

               ELSE
                  CALL PLTXTS(X- (TLEN+XOFF),Y-YOFF,LINE1)
                  IF (TLEN.GT.LONNUM) THEN
                     LONNUM = TLEN
                  END IF

               END IF

            END IF

         END IF

 2630    CONTINUE
         IF (CPUIFC(.FALSE.)) THEN
            GO TO 2650

         END IF

         DO 2660 I = 1,8
            YNEW = YMAJ + LOGTAB(I)*MAJTIC
            IF (PLTNER(YNEW,Y) .AND. NUMSZM.GT.0.0 .AND.
     *          GRAPHP(49).EQ.1.) THEN
               CALL PLTSTT(11,GRAPHP(64))
               CALL PLTSTD(1,GRAPHP(76))
               IF (GRAPHP(40).EQ.2. .OR. GRAPHP(40).EQ.3.) THEN
                  CALL PLTSTT(2,NUMSZM)
                  CALL PLTLOA(LINE1,I+1,0)
                  CALL PLTXSL(LINE1,TLEN)
                  IF (GRAPHP(92).EQ.1.) THEN
                     LONNUM = NUMSIZ
                     CALL PLTSTT(3,90.)
                     CALL PLTXTS(X-XOFF,YNEW-TLEN/2.,LINE1)
                     CALL PLTSTT(3,0.)

                  ELSE
                     CALL PLTXTS(X- (TLEN+XOFF),YNEW-YOFFM,LINE1)
                     IF (TLEN.GT.LONNUM) THEN
                        LONNUM = TLEN
                     END IF

                  END IF

                  CALL PLTSTT(2,NUMSIZ)

               ELSE
                  CALL PLTSTT(2,NUMSIZ)
                  CALL PLTLOD(LINE1,I+1,NUM)
                  CALL PLTXSL(LINE1,TLEN)
                  IF (GRAPHP(92).EQ.1.) THEN
                     LONNUM = NUMSIZ
                     CALL PLTSTT(3,90.)
                     CALL PLTXTS(X-XOFF,YNEW-TLEN/2.,LINE1)
                     CALL PLTSTT(3,0.)

                  ELSE
                     CALL PLTXTS(X- (TLEN+XOFF),YNEW-YOFF,LINE1)
                     IF (TLEN.GT.LONNUM) THEN
                        LONNUM = TLEN
                     END IF

                  END IF

               END IF

            END IF

            IF (YNEW-FUDGE.LE.Y) THEN
               GO TO 2660

            END IF

            IF (YNEW.GT.Y+YLENG-FUDGE) THEN
               GO TO 2670

            END IF

            CALL PLTSTD(1,GRAPHP(77))
            CALL PLTSTV(2,GRAPHP(67))
            CALL PLTVCT(1,X,YNEW,X+LENMIN,YNEW)
            IF (GRAPHP(32).EQ.1.) THEN
               CALL PLTVCT(1,X+XLENG,YNEW,X+XLENG-LENMIN,YNEW)
            END IF

            IF (GRAPHP(73).NE.0.) THEN
               CALL PLTSTV(1,GRAPHP(73))
               CALL PLTSTD(1,GRAPHP(74))
               CALL PLTSTV(2,GRAPHP(69))
               IF (GRAPHP(32).EQ.1.) THEN
                  CALL PLTVCT(1,X+LENMIN,YNEW,X+XLENG-LENMIN,YNEW)

               ELSE
                  CALL PLTVCT(1,X+LENMIN,YNEW,X+XLENG,YNEW)
               END IF

               CALL PLTSTV(1,1.)
            END IF

            IF (GRAPHP(49).EQ.1. .AND. NUMSIZ.GT.0.0) THEN
               CALL PLTSTT(11,GRAPHP(64))
               CALL PLTSTD(1,GRAPHP(76))
               IF (GRAPHP(40).EQ.2. .OR. GRAPHP(40).EQ.3.) THEN
                  CALL PLTSTT(2,NUMSZM)
                  CALL PLTLOA(LINE1,I+1,0)
                  CALL PLTXSL(LINE1,TLEN)
                  IF (GRAPHP(92).EQ.1.) THEN
                     LONNUM = NUMSIZ
                     CALL PLTSTT(3,90.)
                     CALL PLTXTS(X-XOFF,YNEW-TLEN/2.,LINE1)
                     CALL PLTSTT(3,0.)

                  ELSE
                     CALL PLTXTS(X- (TLEN+XOFF),YNEW-YOFF,LINE1)
                     IF (TLEN.GT.LONNUM) THEN
                        LONNUM = TLEN
                     END IF

                  END IF

                  CALL PLTSTT(2,NUMSIZ)

               ELSE
                  CALL PLTSTT(2,NUMSIZ)
                  CALL PLTLOD(LINE1,I+1,NUM)
                  CALL PLTXSL(LINE1,TLEN)
                  IF (GRAPHP(92).EQ.1.) THEN
                     LONNUM = NUMSIZ
                     CALL PLTSTT(3,90.)
                     CALL PLTXTS(X-XOFF,YNEW-TLEN/2.,LINE1)
                     CALL PLTSTT(3,0.)

                  ELSE
                     CALL PLTXTS(X- (TLEN+XOFF),YNEW-YOFF,LINE1)
                     IF (TLEN.GT.LONNUM) THEN
                        LONNUM = TLEN
                     END IF

                  END IF

               END IF

            END IF

 2660    CONTINUE
 2670    CONTINUE
         YMAJ = YMAJ + MAJTIC
         NUM = NUM + 1
         IF (YMAJ.GT.Y+YLENG+FUDGE .AND. GRAPHP(32).EQ.0.) THEN
            GO TO 2650

         END IF

         IF (YMAJ.GT.Y+YLENG-FUDGE .AND. GRAPHP(32).EQ.1.) THEN
            GO TO 2650

         END IF

         CALL PLTSTD(1,GRAPHP(77))
         CALL PLTSTV(2,GRAPHP(67))
         CALL PLTVCT(1,X,YMAJ,X+LENMAJ,YMAJ)
         IF (GRAPHP(32).EQ.1.) THEN
            CALL PLTVCT(1,X+XLENG,YMAJ,X+XLENG-LENMAJ,YMAJ)
         END IF

         IF (GRAPHP(35).NE.0. .OR. GRAPHP(73).NE.0.) THEN
            IF (GRAPHP(35).EQ.0.) THEN
               CALL PLTSTV(1,GRAPHP(73))

            ELSE
               CALL PLTSTV(1,GRAPHP(35))
            END IF

            IF (GRAPHP(35).EQ.0.) THEN
               CALL PLTSTD(1,GRAPHP(74))

            ELSE
               CALL PLTSTD(1,GRAPHP(36))
            END IF

            IF (GRAPHP(35).EQ.0.) THEN
               CALL PLTSTV(2,GRAPHP(69))

            ELSE
               CALL PLTSTV(2,GRAPHP(68))
            END IF

            IF (GRAPHP(32).EQ.1.) THEN
               CALL PLTVCT(1,X+LENMAJ,YMAJ,X+XLENG-LENMAJ,YMAJ)

            ELSE
               CALL PLTVCT(1,X+LENMAJ,YMAJ,X+XLENG,YMAJ)
            END IF

            CALL PLTSTV(1,1.)
         END IF

         IF (NUMSIZ.GT.0.0) THEN
            CALL PLTSTT(2,NUMSIZ)
            CALL PLTSTT(11,GRAPHP(64))
            CALL PLTSTD(1,GRAPHP(76))
            IF (GRAPHP(40).EQ.2. .OR. GRAPHP(40).EQ.3.) THEN
               CALL PLTLOA(LINE1,NUM,0)
               LINE(5:) = LINE1
               CALL PLTXSL(LINE,TLEN)
               IF (GRAPHP(92).EQ.1.) THEN
                  LONNUM = NUMSIZ
                  CALL PLTSTT(3,90.)
                  CALL PLTXTS(X-XOFF,YMAJ-TLEN/2.,LINE)
                  CALL PLTSTT(3,0.)

               ELSE
                  CALL PLTXTS(X- (TLEN+XOFF),YMAJ-YOFF,LINE)
                  IF (TLEN.GT.LONNUM) THEN
                     LONNUM = TLEN
                  END IF

               END IF

               LFLAG = 1

            ELSE
               CALL PLTLOA(LINE1,NUM,1)
               CALL PLTXSL(LINE1,TLEN)
               IF (GRAPHP(92).EQ.1.) THEN
                  LONNUM = NUMSIZ
                  CALL PLTSTT(3,90.)
                  CALL PLTXTS(X-XOFF,YMAJ-TLEN/2.,LINE1)
                  CALL PLTSTT(3,0.)

               ELSE
                  CALL PLTXTS(X- (TLEN+XOFF),YMAJ-YOFF,LINE1)
                  IF (TLEN.GT.LONNUM) THEN
                     LONNUM = TLEN
                  END IF

               END IF

            END IF

         END IF

        GO TO 2630

 2650    CONTINUE
         IF (PLTNER(YNEW,Y+YLENG) .AND. NUMSIZ.GT.0.0 .AND.
     *       GRAPHP(49).EQ.1.) THEN
            CALL PLTSTT(11,GRAPHP(64))
            CALL PLTSTD(1,GRAPHP(76))
            IF (GRAPHP(40).EQ.2. .OR. GRAPHP(40).EQ.3.) THEN
               CALL PLTSTT(2,NUMSZM)
               CALL PLTLOA(LINE1,I+1,0)
               CALL PLTXSL(LINE1,TLEN)
               IF (GRAPHP(92).EQ.1.) THEN
                  LONNUM = NUMSIZ
                  CALL PLTSTT(3,90.)
                  CALL PLTXTS(X-XOFF,YNEW-TLEN/2.,LINE1)
                  CALL PLTSTT(3,0.)

               ELSE
                  CALL PLTXTS(X- (TLEN+XOFF),YNEW-YOFFM,LINE1)
                  IF (TLEN.GT.LONNUM) THEN
                     LONNUM = TLEN
                  END IF

               END IF

            ELSE
               CALL PLTSTT(2,NUMSIZ)
               CALL PLTLOD(LINE1,I+1,NUM-1)
               CALL PLTXSL(LINE1,TLEN)
               IF (GRAPHP(92).EQ.1.) THEN
                  LONNUM = NUMSIZ
                  CALL PLTSTT(3,90.)
                  CALL PLTXTS(X-XOFF,YNEW-TLEN/2.,LINE1)
                  CALL PLTSTT(3,0.)

               ELSE
                  CALL PLTXTS(X- (TLEN+XOFF),YNEW-YOFF,LINE1)
                  IF (TLEN.GT.LONNUM) THEN
                     LONNUM = TLEN
                  END IF

               END IF

            END IF

         END IF

         IF (PLTNER(YMAJ,Y+YLENG) .AND. NUMSIZ.GT.0.0) THEN
            CALL PLTSTT(2,NUMSIZ)
            CALL PLTSTT(11,GRAPHP(64))
            CALL PLTSTD(1,GRAPHP(76))
            IF (GRAPHP(40).EQ.2. .OR. GRAPHP(40).EQ.3.) THEN
               CALL PLTLOA(LINE1,NUM,0)
               LINE(5:) = LINE1
               CALL PLTXSL(LINE,TLEN)
               IF (GRAPHP(92).EQ.1.) THEN
                  LONNUM = NUMSIZ
                  CALL PLTSTT(3,90.)
                  CALL PLTXTS(X-XOFF,YMAJ-TLEN/2.,LINE)
                  CALL PLTSTT(3,0.)

               ELSE
                  CALL PLTXTS(X- (TLEN+XOFF),YMAJ-YOFF,LINE)
                  IF (TLEN.GT.LONNUM) THEN
                     LONNUM = TLEN
                  END IF

               END IF

               LFLAG = 1

            ELSE
               CALL PLTLOA(LINE1,NUM,1)
               CALL PLTXSL(LINE1,TLEN)
               IF (GRAPHP(92).EQ.1.) THEN
                  LONNUM = NUMSIZ
                  CALL PLTSTT(3,90.)
                  CALL PLTXTS(X-XOFF,YMAJ-TLEN/2.,LINE1)
                  CALL PLTSTT(3,0.)

               ELSE
                  CALL PLTXTS(X- (TLEN+XOFF),YMAJ-YOFF,LINE1)
                  IF (TLEN.GT.LONNUM) THEN
                     LONNUM = TLEN
                  END IF

               END IF

            END IF

         END IF

         IF (CPUIFC(.FALSE.)) THEN
            GO TO 2570

         END IF

         IF (LABSIZ.GT.0.0) THEN
            TLABEL = ' '
            LL = 0
            IF (LABEL.NE.' ') THEN
               TLABEL = LABEL
               CALL CHRTRM(TLABEL,LL)
#if NeedsDoubleEscape
               TLABEL(LL+1:LL+2) = '\\-'
#else
               TLABEL(LL+1:LL+2) = '\-'
#endif
               LL = LL + 2
            END IF

            IF (UNITS.NE.' ') THEN
               TLABEL(LL+6:) = UNITS
               CALL CHRTRM(TLABEL,LL)
#if NeedsDoubleEscape
               TLABEL(LL+1:LL+2) = '\\-'
#else
               TLABEL(LL+1:LL+2) = '\-'
#endif
               LL = LL + 2
            END IF

            IF (LFLAG.EQ.0 .AND. (GRAPHP(40).EQ.2..OR.
     *          GRAPHP(40).EQ.3.) .AND. LOWEXP.NE.0 .AND.
     *          NUMSIZ.GT.0.) THEN
               CALL PLTLOA(LINE1,LOWEXP,0)
               CALL CHRTRM(LINE1,L)
#if NeedsDoubleEscape
               TLABEL(LL+1:) = ' (*10\\^'//LINE1(:L)//'\\-)'
#else
               TLABEL(LL+1:) = ' (*10\^'//LINE1(:L)//'\-)'
#endif
            END IF

            CINPUT = TLABEL
            CALL CHRSTR(CINPUT,TLABEL,LL)
            IF (LL.GT.0) THEN
               CALL PLTSTT(2,LABSIZ)
               CALL PLTSTT(11,GRAPHP(65))
               CALL PLTSTD(1,GRAPHP(39))
               CALL PLTSTT(3,90.)
               CALL PLTXSL(TLABEL(1:LL),TLEN)
               XLAB = X - LONNUM - LABSIZ*1.4
               YLAB = Y + (YLENG-TLEN)/2.
               CALL PLTXTS(XLAB,YLAB,TLABEL(1:LL))
               CALL PLTSTT(3,0.)
            END IF

            IF (CPUIFC(.FALSE.)) THEN
               GO TO 2570

            END IF

         END IF

      ELSE
         CALL PLTFLU
         CALL SIORPT('PLTLAX','Invalid axis type - '//TTYPE,2)
      END IF

      IF (.NOT. (.TRUE.)) GO TO 2550
 2570 CONTINUE
      CALL PLTRET
      CALL PLTRED
      CALL PLTREV
      RETURN

      END
