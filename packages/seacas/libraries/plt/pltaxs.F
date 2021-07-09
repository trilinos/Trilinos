C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTAXS(X,Y,XLENG,YLENG,TYPE,XMIN,XMAX,XSTART,NDEC,
     *                  INTER,MININT,LABEL,UNITS,EXP)
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
      INTEGER EXP
      REAL LENMAJ,LENMIN,MAJRAT,MINRAT,NUMRAT,NUMSIZ,LABSIZ,MINTIC,
     *     LABRAT,MAJTIC,LONNUM,LABSCA,INTER
      CHARACTER*(*) LABEL,UNITS,TYPE
      CHARACTER TTYPE*1,FORMA*32,LINE*15,TLABEL*132,EXPL*4,CINPUT*132,
     *          CNDEC*1
      LOGICAL PLTNER
      DATA MAJRAT/50./
      DATA MINRAT/2./
      DATA NUMRAT/50./
      DATA LABRAT/40./
      DATA FUDGE/.0025/

      IF (INTER.EQ.0.) THEN
         CALL PLTFLU
         CALL SIORPT('PLTAXS',
     *               'Interval between major ticks must be nonzero.',2)
         RETURN

      END IF

      TTYPE = TYPE
      CALL CHRUP(TTYPE,TTYPE)
      CALL PLTSVT
      CALL PLTSVD
      CALL PLTSVV
      CALL PLTSTT(3,0.)
      CALL PLTSTT(4,0.)
      CALL PLTSTV(1,1.)
      CALL CHRIC(NDEC,CNDEC,L)
      FORMA = '(f15.'//CNDEC//')'
 2180 CONTINUE
      IF (TTYPE.EQ.'X') THEN
         LABSCA = (XLENG+YLENG)/2.
         NUMSIZ = (LABSCA* (GRAPHP(44)/5.))/NUMRAT
         LABSIZ = (LABSCA* (GRAPHP(45)/5.))/LABRAT
         YOFF = Y - NUMSIZ*1.8
         LENMAJ = XLENG/MAJRAT
         LENMIN = LENMAJ/MINRAT
         MAJTIC = XLENG/ ((XMAX-XMIN)/INTER)
         IF (MININT.NE.0) THEN
            MINTIC = MAJTIC/MININT

         ELSE
            MINTIC = 0.
         END IF

         XLOW = XSTART
         IF (XSTART.NE.XMIN) THEN
            XLOW = XSTART - INTER
         END IF

         FMJTIC = X - (XMIN-XLOW)*MAJTIC/INTER
         XMAJ = FMJTIC
         XNUM = XSTART
         IF (PLTNER(FMJTIC,X) .AND. XSTART.NE.XMIN) THEN
            XNUM = XSTART - INTER
         END IF

         CALL PLTSTD(1,GRAPHP(37))
         CALL PLTSTV(2,GRAPHP(62))
         CALL PLTVCT(1,X,Y,X+XLENG,Y)
         IF (GRAPHP(32).EQ.1.) THEN
            CALL PLTVCT(1,X,Y+YLENG,X+XLENG,Y+YLENG)
         END IF

         IF (PLTNER(FMJTIC,X) .AND. NUMSIZ.GT.0.0) THEN
            LINE = ' '
            WRITE (LINE,FORMA,ERR=120) XNUM
  120       DO 2210 J = 1,15
               IF (LINE(J:J).NE.' ') THEN
                  GO TO 2220

               END IF

 2210       CONTINUE
 2220       CONTINUE
            ILAST = 15
            IF (NDEC.EQ.0) THEN
               ILAST = 14
            END IF

            CALL PLTSTT(2,NUMSIZ)
            CALL PLTSTT(11,GRAPHP(64))
            CALL PLTSTD(1,GRAPHP(76))
            CALL PLTXSL(LINE(J:ILAST),TLEN)
            CALL PLTXTS(X-TLEN/2.,YOFF,LINE(J:ILAST))
            XNUM = XNUM + INTER
         END IF

 2230    CONTINUE
         IF (CPUIFC(.FALSE.)) THEN
            GO TO 2250

         END IF

         DO 2260 I = 1,MININT - 1
            IF (CPUIFC(.FALSE.)) THEN
               GO TO 2270

            END IF

            XNEW = XMAJ + DBLE(I)*MINTIC
            IF (XNEW-FUDGE.LE.X) THEN
               GO TO 2260

            END IF

            IF (XNEW.GT.X+XLENG-FUDGE) THEN
               GO TO 2270

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

 2260    CONTINUE
 2270    CONTINUE
         XMAJ = XMAJ + MAJTIC
         IF (XMAJ.GT.X+XLENG+FUDGE .AND. GRAPHP(32).EQ.0.) THEN
            GO TO 2250

         ELSE IF (XMAJ.GT.X+XLENG-FUDGE .AND. GRAPHP(32).EQ.1.) THEN
            GO TO 2250

         END IF

         CALL PLTSTD(1,GRAPHP(77))
         CALL PLTSTV(2,GRAPHP(67))
         CALL PLTVCT(1,XMAJ,Y,XMAJ,Y+LENMAJ)
         IF (GRAPHP(32).EQ.1.) THEN
            CALL PLTVCT(1,XMAJ,Y+YLENG,XMAJ,Y+YLENG-LENMAJ)
         END IF

         IF (GRAPHP(73).NE.0. .OR. GRAPHP(35).NE.0.) THEN
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
            LINE = ' '
            WRITE (LINE,FORMA,ERR=121) XNUM
  121       DO 2280 J = 1,15
               IF (LINE(J:J).NE.' ') THEN
                  GO TO 2290

               END IF

 2280       CONTINUE
 2290       CONTINUE
            ILAST = 15
            IF (NDEC.EQ.0) THEN
               ILAST = 14
            END IF

            CALL PLTSTT(2,NUMSIZ)
            CALL PLTSTT(11,GRAPHP(64))
            CALL PLTSTD(1,GRAPHP(76))
            CALL PLTXSL(LINE(J:ILAST),TLEN)
            CALL PLTXTS(XMAJ-TLEN/2.,YOFF,LINE(J:ILAST))
         END IF

         XNUM = XNUM + INTER
         GO TO 2230

 2250    CONTINUE
         IF (PLTNER(XMAJ,X+XLENG) .AND. NUMSIZ.GT.0.0) THEN
            LINE = ' '
            WRITE (LINE,FORMA,ERR=122) XNUM
  122       DO 2300 J = 1,15
               IF (LINE(J:J).NE.' ') THEN
                  GO TO 2310

               END IF

 2300       CONTINUE
 2310       CONTINUE
            ILAST = 15
            IF (NDEC.EQ.0) THEN
               ILAST = 14
            END IF

            CALL PLTSTT(2,NUMSIZ)
            CALL PLTSTT(11,GRAPHP(64))
            CALL PLTSTD(1,GRAPHP(76))
            CALL PLTXSL(LINE(J:ILAST),TLEN)
            CALL PLTXTS(XMAJ-TLEN/2.,YOFF,LINE(J:ILAST))
         END IF

         IF ((XMIN.LT.0..AND.XMAX.GT.0.) .AND.
     *       (GRAPHP(35).EQ.0..AND.GRAPHP(73).EQ.0.) .AND.
     *       GRAPHP(71).NE.0.) THEN
            CALL PLTSTV(1,GRAPHP(71))
            CALL PLTSTD(1,GRAPHP(72))
            CALL PLTSTV(2,GRAPHP(70))
            X0LINE = X - (XMIN*XLENG)/ (XMAX-XMIN)
            CALL PLTVCT(1,X0LINE,Y+LENMAJ,X0LINE,Y+YLENG-LENMAJ)
            CALL PLTSTV(1,1.)
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

            IF (EXP.NE.0. .AND. NUMSIZ.GT.0.) THEN
               IF (LL.GT.0) THEN
                  LL = LL + 1
               END IF

               EXPL = ' '
               WRITE (EXPL,'(i4)',ERR=110) EXP
  110          DO 2320 I = 1,4
                  IF (EXPL(I:I).NE.' ') THEN
                     GO TO 2330

                  END IF

 2320          CONTINUE
 2330          CONTINUE
#if NeedsDoubleEscape
               TLABEL(LL+1:) = '( *10\\^'//EXPL(I:)//'\\- )'
#else
               TLABEL(LL+1:) = '( *10\^'//EXPL(I:)//'\- )'
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
               YLAB = YOFF - LABSIZ*2.
               CALL PLTXTS(XLAB,YLAB,TLABEL(1:LL))
            END IF

            IF (CPUIFC(.FALSE.)) THEN
               GO TO 2200

            END IF

         END IF

      ELSE IF (TTYPE.EQ.'Y') THEN
         LABSCA = (XLENG+YLENG)/2.
         NUMSIZ = (LABSCA* (GRAPHP(88)/5.))/NUMRAT
         LABSIZ = (LABSCA* (GRAPHP(89)/5.))/LABRAT
         LONNUM = 0.
         XOFF = NUMSIZ*.8
         YOFF = NUMSIZ/2.
         LENMAJ = YLENG/MAJRAT
         LENMIN = LENMAJ/MINRAT
         MAJTIC = YLENG/ ((XMAX-XMIN)/INTER)
         IF (MININT.NE.0) THEN
            MINTIC = MAJTIC/MININT

         ELSE
            MINTIC = 0.
         END IF

         XLOW = XSTART
         IF (XSTART.NE.XMIN) THEN
            XLOW = XSTART - INTER
         END IF

         FMJTIC = Y - (XMIN-XLOW)*MAJTIC/INTER
         YMAJ = FMJTIC
         XNUM = XSTART
         IF (PLTNER(FMJTIC,Y) .AND. XSTART.NE.XMIN) THEN
            XNUM = XSTART - INTER
         END IF

         CALL PLTSTD(1,GRAPHP(37))
         CALL PLTSTV(2,GRAPHP(62))
         CALL PLTVCT(1,X,Y,X,Y+YLENG)
         IF (GRAPHP(32).EQ.1.) THEN
            CALL PLTVCT(1,X+XLENG,Y,X+XLENG,Y+YLENG)
         END IF

         IF (PLTNER(FMJTIC,Y) .AND. NUMSIZ.GT.0.0) THEN
            LINE = ' '
            WRITE (LINE,FORMA,ERR=130) XNUM
  130       DO 2340 J = 1,15
               IF (LINE(J:J).NE.' ') THEN
                  GO TO 2350

               END IF

 2340       CONTINUE
 2350       CONTINUE
            ILAST = 15
            IF (NDEC.EQ.0) THEN
               ILAST = 14
            END IF

            CALL PLTSTT(2,NUMSIZ)
            CALL PLTSTT(11,GRAPHP(64))
            CALL PLTSTD(1,GRAPHP(76))
            CALL PLTXSL(LINE(J:ILAST),TLEN)
            IF (GRAPHP(92).EQ.1.) THEN
               LONNUM = NUMSIZ
               CALL PLTSTT(3,90.)
               CALL PLTXTS(X-XOFF,Y-TLEN/2.,LINE(J:ILAST))
               CALL PLTSTT(3,0.)

            ELSE
               CALL PLTXTS(X- (TLEN+XOFF),Y-NUMSIZ/2.,LINE(J:ILAST))
               IF (TLEN.GT.LONNUM) THEN
                  LONNUM = TLEN
               END IF

            END IF

            XNUM = XNUM + INTER
         END IF

 2360    CONTINUE
         IF (CPUIFC(.FALSE.)) THEN
            GO TO 2380

         END IF

         DO 2390 I = 1,MININT - 1
            IF (CPUIFC(.FALSE.)) THEN
               GO TO 2400

            END IF

            YNEW = YMAJ + DBLE(I)*MINTIC
            IF (YNEW-FUDGE.LE.Y) THEN
               GO TO 2390

            END IF

            IF (YNEW.GT.Y+YLENG-FUDGE) THEN
               GO TO 2400

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

 2390    CONTINUE
 2400    CONTINUE
         YMAJ = YMAJ + MAJTIC
         IF (YMAJ.GT.Y+YLENG+FUDGE .AND. GRAPHP(32).EQ.0.) THEN
            GO TO 2380

         ELSE IF (YMAJ.GT.Y+YLENG-FUDGE .AND. GRAPHP(32).EQ.1.) THEN
            GO TO 2380

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
            LINE = ' '
            WRITE (LINE,FORMA,ERR=131) XNUM
  131       DO 2410 J = 1,15
               IF (LINE(J:J).NE.' ') THEN
                  GO TO 2420

               END IF

 2410       CONTINUE
 2420       CONTINUE
            ILAST = 15
            IF (NDEC.EQ.0) THEN
               ILAST = 14
            END IF

            CALL PLTSTT(2,NUMSIZ)
            CALL PLTSTT(11,GRAPHP(64))
            CALL PLTSTD(1,GRAPHP(76))
            CALL PLTXSL(LINE(J:ILAST),TLEN)
            IF (GRAPHP(92).EQ.1.) THEN
               LONNUM = NUMSIZ
               CALL PLTSTT(3,90.)
               CALL PLTXTS(X-XOFF,YMAJ-TLEN/2.,LINE(J:ILAST))
               CALL PLTSTT(3,0.)

            ELSE
               CALL PLTXTS(X- (TLEN+XOFF),YMAJ-NUMSIZ/2.,LINE(J:ILAST))
               IF (TLEN.GT.LONNUM) THEN
                  LONNUM = TLEN
               END IF

            END IF

         END IF

         XNUM = XNUM + INTER
         GO TO 2360

 2380    CONTINUE
         IF (PLTNER(YMAJ,Y+YLENG) .AND. NUMSIZ.GT.0.0) THEN
            LINE = ' '
            WRITE (LINE,FORMA,ERR=132) XNUM
  132       DO 2430 J = 1,15
               IF (LINE(J:J).NE.' ') THEN
                  GO TO 2440

               END IF

 2430       CONTINUE
 2440       CONTINUE
            ILAST = 15
            IF (NDEC.EQ.0) THEN
               ILAST = 14
            END IF

            CALL PLTSTT(2,NUMSIZ)
            CALL PLTSTT(11,GRAPHP(64))
            CALL PLTSTD(1,GRAPHP(76))
            CALL PLTXSL(LINE(J:ILAST),TLEN)
            IF (GRAPHP(92).EQ.1.) THEN
               LONNUM = NUMSIZ
               CALL PLTSTT(3,90.)
               CALL PLTXTS(X-XOFF,YMAJ-TLEN/2.,LINE(J:ILAST))
               CALL PLTSTT(3,0.)

            ELSE
               CALL PLTXTS(X- (TLEN+XOFF),YMAJ-NUMSIZ/2.,LINE(J:ILAST))
               IF (TLEN.GT.LONNUM) THEN
                  LONNUM = TLEN
               END IF

            END IF

         END IF

         IF ((XMIN.LT.0..AND.XMAX.GT.0.) .AND.
     *       (GRAPHP(35).EQ.0..AND.GRAPHP(73).EQ.0.) .AND.
     *       GRAPHP(71).NE.0.) THEN
            CALL PLTSTV(1,GRAPHP(71))
            CALL PLTSTD(1,GRAPHP(72))
            CALL PLTSTV(2,GRAPHP(70))
            Y0LINE = Y - (XMIN*YLENG)/ (XMAX-XMIN)
            CALL PLTVCT(1,X+LENMAJ,Y0LINE,X+XLENG-LENMAJ,Y0LINE)
            CALL PLTSTV(1,1.)
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

            IF (EXP.NE.0. .AND. NUMSIZ.GT.0.) THEN
               IF (LL.GT.0) THEN
                  LL = LL + 1
               END IF

               EXPL = ' '
               WRITE (EXPL,'(i4)',ERR=134) EXP
  134          DO 2450 I = 1,4
                  IF (EXPL(I:I).NE.' ') THEN
                     GO TO 2460

                  END IF

 2450          CONTINUE
 2460          CONTINUE
#if NeedsDoubleEscape
               TLABEL(LL+1:) = '( *10\\^'//EXPL(I:)//'\\- )'
#else
               TLABEL(LL+1:) = '( *10\^'//EXPL(I:)//'\- )'
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
               GO TO 2200

            END IF

         END IF

      ELSE
         CALL PLTFLU
         CALL SIORPT('PLTAXS','Invalid axis type - '//TTYPE,2)
      END IF

      IF (.NOT. (.TRUE.)) GO TO 2180
 2200 CONTINUE
      CALL PLTRET
      CALL PLTRED
      CALL PLTREV
      RETURN

      END
