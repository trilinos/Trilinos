C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE REFRSH (MP, ML, MS, MR, MSNAP, N, COOR, ILINE, LTYPE,
     &   LCON, NLPS, IFLINE, ILLIST, NSPR, IFSIDE, ISLIST, LINKP, LINKL,
     &   LINKS, LINKR, TITLE, X1, X2, Y1, Y2, ALPHA, SNAP, SNAPDX,
     &   NSNAP, SNAPDR, TALL, X11, X22, Y11, Y22, AXIST)
C***********************************************************************

C  SUBROUTINE REFRSH = REPAINTS THE DIGITIZING SCREEN WITH SNAP LINES

C***********************************************************************

C  SUBROUTINE CALLED BY:
C     DIGIT = A SUBROUTINE TO DIGITIZE GEOMETRY

C***********************************************************************

      DIMENSION COOR(2, MP), ILINE(ML), LTYPE(ML), LCON(3, ML)
      DIMENSION NLPS(MS), IFLINE(MS), ILLIST(MS*3)
      DIMENSION NSPR(MR), IFSIDE(MR), ISLIST(MR*4)
      DIMENSION LINKP(2, MP), LINKL(2, ML), LINKS(2, MS)
      DIMENSION LINKR(2, MR)
      DIMENSION XDUM(2), YDUM(2)
      DIMENSION N(28), SNAPDX(2, MSNAP), NSNAP(2)

      CHARACTER*72 TITLE

      LOGICAL ALPHA, SNAP, SNAPDR, CPUIFC, AXIST
      LOGICAL ADDLNK, GETMAX, NUMPLT, TEST

      ADDLNK = .FALSE.

      IF (.NOT.ALPHA)THEN
         CALL PLTBGN
         XDIMR = ABS (X2-X1)
         YDIMR = ABS (Y2-Y1)
         XDIMD = 1.
         YDIMD = .75
         CALL MPVIEW (0., XDIMD, 0., YDIMD)
         XRAT = XDIMR/XDIMD
         YRAT = YDIMR/YDIMD
         IF (XRAT.LT.YRAT)THEN
            XDIMR = XDIMD*YRAT
            X11 =  (X1+X2-XDIMR)*.5
            X22 =  (X1+X2+XDIMR)*.5
            XDIMR = X22-X11
            Y11 = Y1
            Y22 = Y2
         ELSE
            YDIMR = YDIMD*XRAT
            Y11 =  (Y1+Y2-YDIMR)*.5
            Y22 =  (Y1+Y2+YDIMR)*.5
            YDIMR = Y22-Y11
            X11 = X1
            X22 = X2
         ENDIF
         IF (AXIST) THEN
            XDUM(1) = X11 - (XDIMR*.05)
            XDUM(2) = X22 + (XDIMR*.05)
            YDUM(1) = Y11 - (YDIMR*.05)
            YDUM(2) = Y22 + (YDIMR*.05)
            SHRINK=.2
         ELSE
            SHRINK=.1
         ENDIF
         X11 = X11 - (XDIMR*SHRINK)
         X22 = X22 + (XDIMR*SHRINK)
         Y11 = Y11 - (YDIMR*SHRINK)
         Y22 = Y22 + (YDIMR*SHRINK)
         CALL MPORT2 (X11, X22, Y11, Y22)
         CALL PLTFRM (0)

C  PLOT THE TITLE

         CALL STRLNG (TITLE, LEN)
         IF ((LEN.GT.1).OR. (TITLE (1:1).NE.' '))THEN
            CALL PLTXHL (TITLE(1:LEN), XLEN)
            XBEGIN = AMAX1 (0., (XDIMD*.5-XLEN*.5))
            CALL PLTXTH (XBEGIN, YDIMD*.95, TITLE(1:LEN))
            IF (CPUIFC (.TRUE.) )RETURN
         ENDIF
         CALL PLTFLU

C  PLOT THE GRID LINES IF APPLICABLE

         IF (SNAP)THEN
            SNAPDR = .TRUE.
            CALL PLTSTD (1, 4.)

C  PLOT THE X GRID LINES

            YBOT = SNAPDX(2, 1)
            YTOP = SNAPDX(2, NSNAP(2))
            DO 100 I = 1, NSNAP(1)
               CALL D2GRID (SNAPDX(1, I), YBOT, SNAPDX(1, I), YTOP)
               IF (CPUIFC (.TRUE.) )RETURN
  100       CONTINUE

C  PLOT THE Y GRID LINES

            XLEFT = SNAPDX(1, 1)
            XRIGHT = SNAPDX(1, NSNAP(1))
            DO 110 I = 1, NSNAP(2)
               CALL D2GRID (XLEFT, SNAPDX(2, I), XRIGHT, SNAPDX(2, I))
               IF (CPUIFC(.TRUE.) )RETURN
  110       CONTINUE

            CALL PLTSTD (1, 7.)
         ELSE
            SNAPDR = .FALSE.
         ENDIF

C  PLOT THE AXIS IF APPLICABLE

         IF (AXIST) CALL SETAXS(XDUM,YDUM)
         IF (CPUIFC(.TRUE.) )RETURN

C  PLOT POINTS ALREADY IN THE DATABASE

         CALL PLTSTT (2, TALL*.5)
         CALL PLTSTD (1, 3.)
         DO 120 I = 1, N(18)
            CALL LTSORT (MP, LINKP, I, II, ADDLNK)
            IF (II.GT.0)THEN
#if NeedsDoubleEscape
               CALL MPD2SY (1, COOR(1, II), COOR(2, II), '\\CX')
#else
               CALL MPD2SY (1, COOR(1, II), COOR(2, II), '\CX')
#endif
            ENDIF
            IF (CPUIFC (.TRUE.) )RETURN
  120    CONTINUE

C  PLOT LINES ALREADY IN THE DATABASE

         KNUM = 0
         CALL PLTSTD (1, 7.)
         COLOR = 7.
         DO 130 I = 1, N(19)
            CALL LTSORT (ML, LINKL, I, II, ADDLNK)
            IF (II.GT.0)THEN
               CALL PLTSTD (1, COLOR)
               LT = LTYPE(II)
               IP1 = LCON(1, II)
               IP2 = LCON(2, II)
               IP3 = LCON(3, II)
               CALL LTSORT (MP, LINKP, IP1, IPNTR1, ADDLNK)
               CALL LTSORT (MP, LINKP, IP2, IPNTR2, ADDLNK)
               CALL LTSORT (MP, LINKP, IABS (IP3), IPNTR3, ADDLNK)
               IF ((IPNTR1.GT.0).AND. (IPNTR2.GT.0).AND.
     +            ((LT.EQ.1).OR. (IPNTR3.GT.0)))THEN
                  CALL DLINE (MP, ML, COOR, LINKP, KNUM, LT, IP1, IP2,
     &               IP3, NUMPLT, DUM1, DUM2, TEST, GETMAX, DUM1, DUM2,
     &               DUM3, DUM4)
                  IF (CPUIFC(.TRUE.) )RETURN
               ENDIF
            ENDIF
  130    CONTINUE

C  PLOT REGIONS ALREADY IN THE DATABASE

         CALL PLTSTD (1, 3.)
         DO 140 I = 1, N(22)
            CALL LTSORT (MR, LINKR, I, II, ADDLNK)
            IF (II.GT.0)THEN
               CALL REGEXT (MP, ML, MS, MR, N, II, COOR, ILINE, LTYPE,
     &            LCON, NLPS, IFLINE, ILLIST, NSPR, IFSIDE, ISLIST,
     &            LINKP, LINKL, LINKS, LINKR, XMIN, XMAX, YMIN, YMAX)
               XMID = (XMAX+XMIN)/2.
               YMID = (YMAX+YMIN)/2.
#if NeedsDoubleEscape
               CALL MPD2SY (1, XMID, YMID, '\\CDI')
#else
               CALL MPD2SY (1, XMID, YMID, '\CDI')
#endif
               IF (CPUIFC(.TRUE.) )RETURN
            ENDIF
  140    CONTINUE
      ENDIF

      RETURN

      END
