C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
      SUBROUTINE DIGIT (MP, ML, MS, MR, MSNAP, MCOM, ICOM, JCOM, CIN,
     &   RIN, IIN, KIN, IDUMP, N, IPOINT, COOR, IPBOUN, ILINE, LTYPE,
     &   NINT, FACTOR, LCON, ILBOUN, ISBOUN, ISIDE, NLPS, IFLINE,
     &   ILLIST, IREGN, IMAT, NSPR, IFSIDE, ISLIST, IRPB, IPBF, NPPF,
     &   IFPB, LISTPB, ILBF, NLPF, IFLB, LISTLB, ISBF, NSPF, IFSB,
     &   LISTSB, LINKP, LINKL, LINKS, LINKR, LINKM, LINKPB, LINKLB,
     &   LINKSB, IHOLDP, IHOLDL, IHOLDS, IHOLDR, IHOLDM, IHOLD2, IHOLD3,
     &   IWTPBF, IWTLBF, IWTSBF, IRGFLG, TITLE, NOROOM, XX1, YY1, SCALE,
     &   CT, ST, X1, X2, Y1, Y2, X11, X22, Y11, Y22, XMIN1, XMAX1,
     &   YMIN1, YMAX1, XMIN2, XMAX2, YMIN2, YMAX2, X1OLD, X2OLD, Y1OLD,
     &   Y2OLD, ALPHA, DEV1, SNAP, SNAPDX, NSNAP, DRWTAB, AXIST)
C***********************************************************************

C  SUBROUTINE DIGIT = SUBROUTINE TO DIGITIZE THE GEOMETRY

C***********************************************************************

C  SUBROUTINE CALLED BY:
C     TABLET = A SUBROUTINE TO CONTROL THE GEOMETRY DIGITIZATION

C***********************************************************************

C  SUBROUTINES CALLED:
C     DREAD  = SETS ALL PARAMETERS UP FOR READING FROM A DIGI-PAD
C     DPREAD = READS INPUT FROM A DIGI-PAD DIGITIZER
C     CLOSE  = FINDS THE CLOSEST EXISTING POINT TO THE MOUSE
C     REFRSH = REFRESHES THE SCREEN AND DRAWS EXISTING GEOMETRY

C***********************************************************************

C  VARIABLES USED:
C     IANS   = LOGICAL RESPONSE FROM YES-NO QUESTION
C     TITLE  = MESH TITLE
C     XX1    = DIGITIZATION PAD X COORDINATE OF POINT 1 (PAD INIT)
C     YY1    = DIGITIZATION PAD Y COORDINATE OF POINT 1 (PAD INIT)
C     X1     = USER X COORDINATE OF POINT 1 (PAD INIT)
C     Y1     = USER Y COORDINATE OF POINT 1 (PAD INIT)
C     X2     = USER X COORDINATE OF POINT 2 (PAD INIT)
C     Y2     = USER Y COORDINATE OF POINT 2 (PAD INIT)
C     X      = THE X LOCATION IN USER COORDINATES
C     Y      = THE Y LOCATION IN USER COORDINATES
C     BUTTON = THE MOUSE BUTTON PUSHED
C     SCALE  = THE SCALE FACTOR FROM DIGITIZED TO USER COORDINATES
C     CT     = THE COSINE OF THE ANGLE OF THE DRAWING ON THE PAD
C     ST     = THE SINE OF THE ANGLE OF THE DRAWING ON THE PAD
C     CHANGE = .TRUE. IF THE ZERO BUTTON WAS PUSHED LAST
C     SLIDE  = .TRUE. IF THE NEXT POINT IS TO HAVE THE CLOSEST POINT'S
C              COORDINATES, BUT NEW NUMBERING (SLIDE LINE USE) (C BUTTON
C              HAS BEEN PUSHED)
C     NOROOM = .TRUE. IF THE AMOUNT OF DATA EXCEEDS DIMENSIONED LIMITS

C***********************************************************************

      DIMENSION IPOINT(MP), COOR(2, MP), IPBOUN(MP)
      DIMENSION ILINE(ML), LTYPE(ML), NINT(ML), FACTOR(ML), LCON(3, ML)
      DIMENSION ILBOUN(ML), ISBOUN(ML)
      DIMENSION ISIDE(MS), NLPS(MS), IFLINE(MS), ILLIST(MS*3)
      DIMENSION IREGN(MR), IMAT(MR), NSPR(MR), IFSIDE(MR), ISLIST(MR*4)
      DIMENSION IRPB(MR), IPBF(MP), NPPF(MP), IFPB(MP), LISTPB(2, MP)
      DIMENSION IWTPBF(3, MP), ILBF(ML), NLPF(ML), IFLB(ML)
      DIMENSION LISTLB(2, ML), IWTLBF(3, ML), ISBF(ML), NSPF(ML)
      DIMENSION IFSB(ML), LISTSB(2, ML), IWTSBF(3, ML)
      DIMENSION LINKP(2, MP), LINKL(2, ML), LINKS(2, MS), LINKR(2, MR)
      DIMENSION LINKM(2, (MS + MR)), LINKPB(2, MP), LINKLB(2, ML)
      DIMENSION LINKSB(2, ML)
      DIMENSION IHOLDP(2, MP), IHOLDL(ML*2), IHOLDR(2, MR)
      DIMENSION IHOLDM(2, (MS + MR)), IHOLD2(2, ML), IHOLD3(2, ML)
      DIMENSION IHOLDS(2, MS), IRGFLG(MR)
      DIMENSION N(29), SNAPDX(2, MSNAP), NSNAP(2)
      DIMENSION KIN(MCOM), IIN(MCOM), RIN(MCOM)

      CHARACTER*72 TITLE, CIN(MCOM)
      CHARACTER*1 BUTTON, HOLD
      CHARACTER DEV1*3

      LOGICAL CHANGE, SLIDE, NOROOM, BOXED, ADDOLD, NUMPLT
      LOGICAL DRWTAB, OLDCUR, MERGE, ALPHA, ADDLNK, BIFIND, TEST
      LOGICAL GETMAX, SNAP, SNAPDR, AXIST, ADDCEN

C  INITIALIZE VARIABLES

      XMIN2 = XMIN1
      XMAX2 = XMAX1
      YMIN2 = YMIN1
      YMAX2 = YMAX1
      X1OLD = X1
      X2OLD = X2
      Y1OLD = Y1
      Y2OLD = Y2
      ADDLNK = .FALSE.
      MERGE = .FALSE.
      TEST = .FALSE.
      GETMAX = .FALSE.
      NUMPLT = .FALSE.
      NOLD7 = N(7)
      CALL PLTGTT (2, TALL)
      knum = 0

      IF (NOROOM) THEN
         CALL PLTBEL
         CALL PLTFLU
         NOROOM = .FALSE.
      END IF

C  INITIALIZE THE PLOTTING SURFACE

      IF ((DEV1 .EQ. 'LS5') .OR. (DEV1 .EQ. 'ls5')) THEN
         WRITE (*, 10000) CHAR(27)//'RA0'
         WRITE (*, 10000) CHAR(27)//'RR1'
      END IF

C  REFRESH THE SCREEN AND DRAW EXISTING DATA WITH SNAP LINES

      IF (DRWTAB) THEN
         CALL REFRSH(MP, ML, MS, MR, MSNAP, N, COOR, ILINE, LTYPE,
     &      LCON, NLPS, IFLINE, ILLIST, NSPR, IFSIDE, ISLIST, LINKP,
     &      LINKL, LINKS,LINKR, TITLE, X1OLD, X2OLD, Y1OLD, Y2OLD,
     &      ALPHA, SNAP, SNAPDX, NSNAP, SNAPDR, TALL, X11, X22, Y11,
     &      Y22, AXIST)
      ELSE
         CALL REFRSH(MP, ML, MS, MR, MSNAP, N, COOR, ILINE, LTYPE,
     &      LCON, NLPS, IFLINE, ILLIST, NSPR, IFSIDE, ISLIST, LINKP,
     &      LINKL, LINKS, LINKR, TITLE, X1, X2, Y1, Y2, ALPHA, SNAP,
     &      SNAPDX, NSNAP, SNAPDR, TALL, X11, X22, Y11, Y22, AXIST)
      ENDIF

C  INITIALIZE VARIABLES FOR DIGITIZATION INPUT

      LASTP = 0
      OLDCUR = .FALSE.

  100 CONTINUE
      CHANGE = .FALSE.
      SLIDE = .FALSE.
  110 CONTINUE
      IF (.NOT.ALPHA) THEN
         CALL PLTSTD (1, 7.)
         CALL PLTFLU
      END IF

C  INPUT THE DATA

C  GET THE CURSOR LOCATION IN USER SYSTEM COORDINATES
C  AND THE BUTTON PUSHED

      CALL DREAD (X, Y, BUTTON, X1, Y1, XX1, YY1, SCALE, CT, ST)

C  SHOW CURRENT CURSOR LOCATION

      IF ((.NOT.ALPHA) .AND. (BUTTON .EQ. '3')) THEN

C  ERASE ANY BOXED IN POINTS THAT ARE NOT NEEDED FOR THE
C  CURRENT SEQUENCE

         IF (ISQR .LT. 0) THEN
            ISQR = IABS (ISQR)
            CALL PLTSTT (2, TALL)
            CALL PLTSTD (1, 0.)
#if NeedsDoubleEscape
            CALL MPD2SY (1, COOR(1, ISQR), COOR(2, ISQR), '\\CSQ')
#else
            CALL MPD2SY (1, COOR(1, ISQR), COOR(2, ISQR), '\CSQ')
#endif
            CALL PLTSTD (1, 7.)
            CALL PLTSTT (2, TALL*.5)
            CALL PLTFLU
            ISQR = 0
         ENDIF
         CALL PUTCRS (X, Y, OLDCUR)
         GOTO 110

C  ERASE ANY BOXED IN POINTS ON THE SCREEN

      ELSEIF (ISQR .NE. 0) THEN
         ISQR = IABS (ISQR)
         CALL PLTSTT (2, TALL)
         CALL PLTSTD (1, 0.)
#if NeedsDoubleEscape
         CALL MPD2SY (1, COOR(1, ISQR), COOR(2, ISQR), '\\CSQ')
#else
         CALL MPD2SY (1, COOR(1, ISQR), COOR(2, ISQR), '\CSQ')
#endif
         CALL PLTSTD (1, 7.)
         CALL PLTSTT (2, TALL*.5)
         CALL PLTFLU
         ISQR = 0
      END IF

      IF ((DEV1 .EQ. 'LS5') .OR. (DEV1 .EQ. 'ls5')) THEN
         WRITE (*, 10000) CHAR(27)//'[2J'
         WRITE (*, 10000) CHAR(27)//'RA0'
      END IF

C  INPUT A POINT

      IF (BUTTON .EQ. '1') THEN
         IF (.NOT.ALPHA) CALL ERASEC(OLDCUR)

C  FIND THE CLOSEST POINT TO THE CURRENT POSITION

         IF (CHANGE) THEN
            CALL CLOSEP (MP, N(18), X, Y, IPOINT, COOR, LINKP, I)

C  INPUT A POINT AT THE CLOSEST LOCATION, BUT WITH A NEW POINT NUMBER

            CALL LTSORT (MP, LINKP, I, II, ADDLNK)
            IF (SLIDE) THEN
               JJ = N(18) + 1
               CALL INPOIN (MP, N(1), N(18), JJ, COOR(1, II),
     &            COOR(2, II), NHOLDP, IHOLDP, IPOINT, COOR, IPBOUN,
     &            LINKP, MERGE, NOROOM)
               IF (NOROOM) GO TO 170
               I = JJ
               IF (.NOT.ALPHA) THEN
                  CALL PLTSTD (1, 3.)
#if NeedsDoubleEscape
                  CALL MPD2SY (1, COOR(1, II), COOR(2, II), '\\CX')
#else
                  CALL MPD2SY (1, COOR(1, II), COOR(2, II), '\CX')
#endif
                  CALL PLTSTD (1, 7.)
                  CALL PLTFLU
               END IF
               CALL DMESS(DEV1,' A NEW POINT HAS BEEN INPUT ON TOP '//
     &            'OF CLOSEST POINT')
            ELSE
               CALL PLTSTT (2, TALL)
               CALL PLTSTD (1, 3.)
#if NeedsDoubleEscape
               CALL MPD2SY (1, COOR(1, II), COOR(2, II), '\\CSQ')
#else
               CALL MPD2SY (1, COOR(1, II), COOR(2, II), '\CSQ')
#endif
               CALL PLTSTD (1, 7.)
               CALL PLTSTT (2, TALL*.5)
               CALL PLTFLU
               ISQR = II
               CALL DMESS(DEV1,' CLOSEST POINT LOCATED AND STORED')
            END IF
            LASTP = I

C  INPUT A COMPLETELY NEW POINT

         ELSE
            IF (SNAP) CALL SNAPPT (MSNAP, SNAPDX, NSNAP, X, Y)
            JJ = N(18) + 1
            CALL INPOIN (MP, N(1), N(18), JJ, X, Y, NHOLDP, IHOLDP,
     &         IPOINT, COOR, IPBOUN, LINKP, MERGE, NOROOM)
            IF (NOROOM) GO TO 170
            LASTP = JJ
            CALL LTSORT (MP, LINKP, JJ, IPNTR, ADDLNK)
            IF (.NOT.ALPHA) THEN
               CALL PLTSTD (1, 3.)
#if NeedsDoubleEscape
               CALL MPD2SY (1, COOR(1, IPNTR), COOR(2, IPNTR), '\\CX')
#else
               CALL MPD2SY (1, COOR(1, IPNTR), COOR(2, IPNTR), '\CX')
#endif
               CALL PLTSTD (1, 7.)
               CALL PLTFLU
            END IF
            CALL DMESS(DEV1,' A POINT HAS BEEN INPUT')
         END IF

C  INPUT A STRAIGHT LINE

      ELSE IF (BUTTON .EQ. '2') THEN
         IF (.NOT.ALPHA) CALL ERASEC (OLDCUR)

C  CONNECT THIS LINE FROM THE LAST POINT TO THE CLOSEST POINT

         IF (CHANGE) THEN
            CALL CLOSEP (MP, N(18), X, Y, IPOINT, COOR, LINKP, I)
            IF (I .EQ. LASTP) GO TO 100

C  RENUMBER THIS LINE BUT PUT IT ON TOP OF THE EXISTING LINE
C  THIS IS FOR SLIDE LINE USE

            IF (SLIDE) THEN
               JJ = N(18) + 1
               CALL LTSORT (MP, LINKP, I, II, ADDLNK)
               CALL INPOIN (MP, N(1), N(18), JJ, COOR(1, II),
     &            COOR(2, II), NHOLDP, IHOLDP, IPOINT, COOR, IPBOUN,
     &            LINKP, MERGE, NOROOM)
               IF (NOROOM) GO TO 170
               I = JJ
               IF (.NOT.ALPHA) THEN
                  CALL PLTSTD (1, 3.)
#if NeedsDoubleEscape
                  CALL MPD2SY (1, COOR(1, II), COOR(2, II), '\\CX')
#else
                  CALL MPD2SY (1, COOR(1, II), COOR(2, II), '\CX')
#endif
                  CALL PLTSTD (1, 7.)
                  CALL PLTFLU
               END IF
               CALL DMESS(DEV1,' A STRAIGHT LINE HAS BEEN INPUT '//
     &            'TO NEW POINT ON TOP OF CLOSEST POINT')
            ELSE
               CALL DMESS(DEV1,' A STRAIGHT LINE HAS BEEN INPUT '//
     &            'TO CLOSEST POINT')
            END IF

C  INPUT A COMPLETELY NEW STRAIGHT LINE

         ELSE
            JJ = N(18) + 1
            IF (SNAP) CALL SNAPPT (MSNAP, SNAPDX, NSNAP, X, Y)
            CALL INPOIN (MP, N(1), N(18), JJ, X, Y, NHOLDP, IHOLDP,
     &         IPOINT, COOR, IPBOUN, LINKP, MERGE, NOROOM)
            IF (NOROOM) GO TO 170
            I = JJ
            CALL LTSORT (MP, LINKP, I, II, ADDLNK)
            IF (.NOT.ALPHA) THEN
               CALL PLTSTD (1, 3.)
#if NeedsDoubleEscape
               CALL MPD2SY (1, COOR(1, II), COOR(2, II), '\\CX')
#else
               CALL MPD2SY (1, COOR(1, II), COOR(2, II), '\CX')
#endif
               CALL PLTSTD (1, 7.)
               CALL PLTFLU
            END IF
            CALL DMESS(DEV1,' A STRAIGHT LINE HAS BEEN INPUT '//
     &         'TO NEW POINT')
         END IF
         IF (LASTP .EQ. 0) THEN
            LASTP = I
            GO TO 100
         END IF
         LT = 1
         IP3 = 0
         NN = 0
         FACT = 1.
         JJ = N(19) + 1
         CALL INLINE (ML, N(2), N(19), JJ, LT, LASTP, I, IP3, NN, FACT,
     &      NHOLDL, IHOLDL, ILINE, LTYPE, NINT, FACTOR, LCON, ILBOUN,
     &      ISBOUN, LINKL, MERGE, NOROOM)
         IF (NOROOM) GO TO 170
         IF (.NOT.ALPHA) THEN
            CALL DLINE (MP, ML, COOR, LINKP, KNUM, LT, LASTP, I, IP3,
     &         NUMPLT, DUM1, DUM2, TEST, GETMAX, DUM3, DUM4, DUM5, DUM6)
            CALL PLTFLU
         END IF
         LASTP = I

C  INPUT AN ARC ABOUT A CENTER POINT (NEED NOT BE TRULY CIRCULAR)

      ELSE IF ((BUTTON .EQ. '5') .OR. (BUTTON .EQ. '6')) THEN
         IF (.NOT.ALPHA) CALL ERASEC (OLDCUR)
         HOLD = BUTTON

C  GO FROM THE LAST POINT TO THE CLOSEST EXISTING POINT FOR THE ARC

         IF (CHANGE) THEN
            CALL CLOSEP (MP, N(18), X, Y, IPOINT, COOR, LINKP, I)
            IF (.NOT. ALPHA) THEN
               CALL PLTSTT (2, TALL)
               CALL PLTSTD (1, 3.)
               CALL LTSORT (MP, LINKP, I, IPNTR, ADDLNK)
#if NeedsDoubleEscape
               CALL MPD2SY (1, COOR(1, IPNTR), COOR(2, IPNTR), '\\CSQ')
#else
               CALL MPD2SY (1, COOR(1, IPNTR), COOR(2, IPNTR), '\CSQ')
#endif
               CALL PLTSTD (1, 7.)
               CALL PLTSTT (2, TALL*.5)
               CALL PLTFLU
               ISQR = IPNTR
            ENDIF

C  INPUT THE LINE AS NEW NUMBER IN THE SAME PLACE AS THE OLD ARC
C  THIS IS FOR SLIDE LINE USE

            IF (SLIDE) THEN
               JJ = N(18) + 1
               CALL LTSORT (MP, LINKP, I, II, ADDLNK)
               CALL INPOIN (MP, N(1), N(18), JJ, COOR(1, II),
     &            COOR(2, II), NHOLDP, IHOLDP, IPOINT, COOR, IPBOUN,
     &            LINKP, MERGE, NOROOM)
               IF (NOROOM) GO TO 170
               I = JJ
               IF (.NOT.ALPHA) THEN
                  CALL PLTSTD (1, 3.)
#if NeedsDoubleEscape
                  CALL MPD2SY (1, COOR(1, II), COOR(2, II), '\\CX')
#else
                  CALL MPD2SY (1, COOR(1, II), COOR(2, II), '\CX')
#endif
                  CALL PLTSTD (1, 7.)
                  CALL PLTFLU
               END IF
            END IF

C  INPUT A COMPLETELY NEW ARC

         ELSE
            JJ = N(18) + 1
            IF (SNAP) CALL SNAPPT (MSNAP, SNAPDX, NSNAP, X, Y)
            CALL INPOIN (MP, N(1), N(18), JJ, X, Y, NHOLDP, IHOLDP,
     &         IPOINT, COOR, IPBOUN, LINKP, MERGE, NOROOM)
            IF (NOROOM) GO TO 170
            I = JJ
            CALL LTSORT (MP, LINKP, I, II, ADDLNK)
            IF (.NOT.ALPHA) THEN
               CALL PLTSTD (1, 3.)
#if NeedsDoubleEscape
               CALL MPD2SY (1, COOR(1, II), COOR(2, II), '\\CX')
#else
               CALL MPD2SY (1, COOR(1, II), COOR(2, II), '\CX')
#endif
               CALL PLTSTD (1, 7.)
               CALL PLTFLU
            END IF
         END IF
         IF (LASTP .EQ. 0) THEN
            LASTP = I
            GO TO 100
         END IF

C  INPUT THE CENTER POINT OF THE ARC

  120    CONTINUE
         CALL DMESS (DEV1, ' INPUT A CENTER POINT TO FINISH THE ARC'//
     &      ' DEFINITION')
         CALL MESAGE(' "PUCK-1" FOR NEW OR '//
     &      '"PUCK-0 PUCK-1" FOR EXISTING CENTER POINT')
         CALL DREAD (X, Y, BUTTON, X1, Y1, XX1, YY1, SCALE, CT, ST)

C  SHOW CURRENT CURSOR LOCATION AS AN AID IN LOCATION OF THE CENTER

         IF ((.NOT.ALPHA) .AND. (BUTTON .EQ. '3')) THEN
            CALL PUTCRS (X, Y, OLDCUR)
            GO TO 120

C  USE THE CLOSEST POINT TO CURRENT LOCATION AS THE CENTER

         ELSE IF (BUTTON .EQ. '0') THEN
            CALL DREAD (X, Y, BUTTON, X1, Y1, XX1, YY1, SCALE, CT, ST)
            WRITE(*,10000)' '//CHAR(27)//'[2J'
            IF (BUTTON .NE. '1') GO TO 100
            CALL CLOSEP (MP, N(18), X, Y, IPOINT, COOR, LINKP, II)
            IF (.NOT. ALPHA) THEN
               IF (ISQR .NE. 0) THEN
                  CALL PLTSTT (2, TALL)
                  CALL PLTSTD (1, 0.)
#if NeedsDoubleEscape
                  CALL MPD2SY (1, COOR(1, ISQR), COOR(2, ISQR), '\\CSQ')
#else
                  CALL MPD2SY (1, COOR(1, ISQR), COOR(2, ISQR), '\CSQ')
#endif
                  CALL PLTSTD (1, 7.)
                  CALL PLTSTT (2, TALL*.5)
                  CALL PLTFLU
                  ISQR = 0
               END IF
               CALL PLTSTT (2, TALL)
               CALL PLTSTD (1, 3.)
               CALL LTSORT (MP, LINKP, II, IPNTR, ADDLNK)
#if NeedsDoubleEscape
               CALL MPD2SY (1, COOR(1, IPNTR), COOR(2, IPNTR), '\\CSQ')
#else
               CALL MPD2SY (1, COOR(1, IPNTR), COOR(2, IPNTR), '\CSQ')
#endif
               CALL PLTSTD (1, 7.)
               CALL PLTSTT (2, TALL*.5)
               CALL PLTFLU
               ISQR = IPNTR
            ENDIF

C  INPUT A NEW POINT FOR THE CENTER

         ELSE IF (BUTTON .EQ. '1') THEN
            WRITE(*,10000)' '//CHAR(27)//'[2J'
            IF (.NOT.ALPHA) CALL ERASEC (OLDCUR)
            JJ = N(18) + 1
            IF (SNAP) CALL SNAPPT (MSNAP, SNAPDX, NSNAP, X, Y)
            CALL INPOIN (MP, N(1), N(18), JJ, X, Y, NHOLDP, IHOLDP,
     &         IPOINT, COOR, IPBOUN, LINKP, MERGE, NOROOM)
            IF (NOROOM) GO TO 170
            II = JJ
            CALL LTSORT (MP, LINKP, II, IPNTR, ADDLNK)
            IF (.NOT.ALPHA) THEN
               CALL PLTSTD (1, 3.)
#if NeedsDoubleEscape
               CALL MPD2SY (1, COOR(1, IPNTR), COOR(2, IPNTR), '\\CX')
#else
               CALL MPD2SY (1, COOR(1, IPNTR), COOR(2, IPNTR), '\CX')
#endif
               CALL PLTSTD (1, 7.)
               CALL PLTFLU
            END IF
         ELSE
            WRITE(*,10000)' '//CHAR(27)//'[2J'
            GO TO 100
         END IF

C  IF THE ARC IS INPUT AS CW, SWITCH END POINTS TO RECORD IT AS CCW

         IF (HOLD .EQ. '5') THEN
            IP1 = LASTP
            IP2 = I
         ELSE
            IP1 = I
            IP2 = LASTP
         END IF
         LT = 3
         IP3 = II
         IF ((IP3 .EQ. IP2) .OR. (IP3 .EQ. IP1)) THEN
            CALL PLTBEL
            CALL PLTFLU
            GO TO 100
         END IF
         NN = 0
         FACT = 1.
         JJ = N(19) + 1
         CALL INLINE (ML, N(2), N(19), JJ, LT, IP1, IP2, IP3, NN, FACT,
     &      NHOLDL, IHOLDL, ILINE, LTYPE, NINT, FACTOR, LCON, ILBOUN,
     &      ISBOUN, LINKL, MERGE, NOROOM)
         IF (NOROOM) GO TO 170
         IF (.NOT.ALPHA) CALL DLINE (MP, ML, COOR, LINKP, KNUM, LT, IP1,
     &      IP2, IP3, NUMPLT, DUM1, DUM2, TEST, GETMAX, DUM3, DUM4,
     &      DUM5, DUM6)
         LASTP = I

C  WRITE OUT THE APPROPRIATE MESSAGE

         IF (HOLD .EQ. '5') THEN
            IF (CHANGE) THEN
               IF (SLIDE) THEN
                  CALL DMESS(DEV1,' A CCW ARC LINE HAS BEEN INPUT '//
     &               'TO NEW POINT ON TOP OF CLOSEST POINT')
               ELSE
                  CALL DMESS(DEV1,' A CCW ARC LINE HAS BEEN INPUT '//
     &               'TO CLOSEST POINT')
               END IF
            ELSE
               CALL DMESS(DEV1,' A CCW ARC LINE HAS BEEN INPUT '//
     &            'TO NEW POINT')
            ENDIF
         ELSE
            IF (CHANGE) THEN
               IF (SLIDE) THEN
                  CALL DMESS(DEV1,' A CW ARC LINE HAS BEEN INPUT '//
     &               'TO NEW POINT ON TOP OF CLOSEST POINT')
               ELSE
                  CALL DMESS(DEV1,' A CW ARC LINE HAS BEEN INPUT '//
     &               'TO CLOSEST POINT')
               END IF
            ELSE
               CALL DMESS(DEV1,' A CW ARC LINE HAS BEEN INPUT '//
     &            'TO NEW POINT')
            ENDIF
         ENDIF

C  BISECT THE CLOSEST LINE

      ELSE IF (BUTTON .EQ. '4') THEN
         IF (.NOT.ALPHA) CALL ERASEC (OLDCUR)
         CALL BISECT (MP, ML, MS, MR, MSNAP, N, IPOINT, COOR, IPBOUN,
     &      ILINE, LTYPE, NINT, FACTOR, LCON, ILBOUN, ISBOUN, ISIDE,
     &      NLPS, IFLINE, ILLIST, IREGN, IMAT, NSPR, IFSIDE, ISLIST,
     &      IRPB, IPBF, NPPF, IFPB, LISTPB, ILBF, NLPF, IFLB, LISTLB,
     &      ISBF, NSPF, IFSB, LISTSB, LINKP, LINKL, LINKS, LINKR, LINKM,
     &      LINKPB, LINKLB, LINKSB, IHOLDP, IHOLDL, IHOLDS, IHOLDR,
     &      IHOLDM, IHOLD2, IHOLD3, IWTPBF, IWTLBF, IWTSBF, IRGFLG, X,
     &      Y, MERGE, LASTP, ALPHA, NOROOM, SNAP, SNAPDX, NSNAP)
         CALL  DMESS(DEV1,' AN EXISTING LINE HAS BEEN BISECTED')

C  MOVE A POINT AND ALL ENTITIES ATTACHED TO THAT POINT

      ELSE IF (BUTTON .EQ. '8') THEN
         CALL CLOSEP (MP, N(18), X, Y, IPOINT, COOR, LINKP, II)
         IF (.NOT.ALPHA) THEN
            CALL PLTSTT (2, TALL)
            CALL PLTSTD (1, 3.)
            CALL LTSORT (MP, LINKP, II, IPNTR, ADDLNK)
#if NeedsDoubleEscape
            CALL MPD2SY (1, COOR(1, IPNTR), COOR(2, IPNTR), '\\CSQ')
#else
            CALL MPD2SY (1, COOR(1, IPNTR), COOR(2, IPNTR), '\CSQ')
#endif
            CALL PLTSTT (2, TALL*.5)
            ISQR = IPNTR
            CALL PLTFLU
         END IF
  130    CONTINUE
         CALL DMESS(DEV1,' ENTER "PUCK-8" TO MOVE BOXED POINT - '//
     &      'ANY OTHER KEY TO ABANDON')

C  INPUT THE NEW POINT LOCATION

         CALL DREAD (X, Y, BUTTON, X1, Y1, XX1, YY1, SCALE, CT, ST)

C  SHOW CURRENT CURSOR LOCATION AS AN AID IN LOCATION OF THE POINT

         IF ((.NOT.ALPHA) .AND. (BUTTON .EQ. '3')) THEN
            CALL PUTCRS (X, Y, OLDCUR)
            GO TO 130

C  INPUT A NEW X, Y FOR THE POINT

         ELSE IF (BUTTON .EQ. '8') THEN
            WRITE(*,10000)' '//CHAR(27)//'[2J'
            IF (.NOT.ALPHA) CALL ERASEC (OLDCUR)
            IF (SNAP) CALL SNAPPT (MSNAP, SNAPDX, NSNAP, X, Y)
            CALL PLTSTT (2, TALL)
            CALL PLTSTD (1, 0.)
#if NeedsDoubleEscape
            CALL MPD2SY (1, COOR(1, ISQR), COOR(2, ISQR), '\\CSQ')
#else
            CALL MPD2SY (1, COOR(1, ISQR), COOR(2, ISQR), '\CSQ')
#endif
            CALL PLTSTD (1, 7.)
            CALL PLTSTT (2, TALL*.5)
            CALL PLTFLU
            ISQR = 0
            CALL MOVEP (MP, ML, MS, MR, N, COOR, LTYPE, LCON, NLPS,
     &         IFLINE, ILLIST, NSPR, IFSIDE, ISLIST, LINKP, LINKL,
     &         LINKS, LINKR, II, ALPHA, X, Y)
            LASTP = II
         ELSE
            CALL PLTSTT (2, TALL)
            CALL PLTSTD (1, 0.)
#if NeedsDoubleEscape
            CALL MPD2SY (1, COOR(1, ISQR), COOR(2, ISQR), '\\CSQ')
#else
            CALL MPD2SY (1, COOR(1, ISQR), COOR(2, ISQR), '\CSQ')
#endif
            CALL PLTSTD (1, 7.)
            CALL PLTSTT (2, TALL*.5)
            CALL PLTFLU
            WRITE(*,10000)' '//CHAR(27)//'[2J'
            ISQR = 0
         END IF

C  REFRESH THE SCREEN AND DRAW EXISTING DATA WITH SNAP LINES

      ELSE IF (BUTTON .EQ. '9') THEN
         IF (.NOT.ALPHA) CALL ERASEC (OLDCUR)
         IF (DRWTAB) THEN
            CALL REFRSH (MP, ML, MS, MR, MSNAP, N, COOR, ILINE, LTYPE,
     &         LCON, NLPS, IFLINE, ILLIST, NSPR, IFSIDE, ISLIST, LINKP,
     &         LINKL, LINKS, LINKR, TITLE, X1OLD, X2OLD, Y1OLD, Y2OLD,
     &         ALPHA, SNAP, SNAPDX, NSNAP, SNAPDR, TALL, X11, X22, Y11,
     &         Y22, AXIST)
         ELSE
            CALL REFRSH (MP, ML, MS, MR, MSNAP, N, COOR, ILINE, LTYPE,
     &         LCON, NLPS, IFLINE, ILLIST, NSPR, IFSIDE, ISLIST, LINKP,
     &         LINKL, LINKS, LINKR, TITLE, X1, X2, Y1, Y2, ALPHA, SNAP,
     &         SNAPDX, NSNAP, SNAPDR, TALL, X11, X22, Y11, Y22, AXIST)
         ENDIF

C  TOGGLE THE SNAPPING TO GRID LINES

      ELSE IF (BUTTON .EQ. 'A') THEN
         IF (.NOT.ALPHA) CALL ERASEC (OLDCUR)
         IF (SNAP) THEN
            SNAP = .FALSE.
            CALL DMESS(DEV1,' SNAP TO GRID POINTS DISABLED')
            IF (SNAPDR) THEN
               IF (DRWTAB) THEN
                  CALL REFRSH (MP, ML, MS, MR, MSNAP, N, COOR, ILINE,
     &               LTYPE, LCON, NLPS, IFLINE, ILLIST, NSPR, IFSIDE,
     &               ISLIST, LINKP, LINKL, LINKS, LINKR, TITLE, X1OLD,
     &               X2OLD, Y1OLD, Y2OLD, ALPHA, SNAP, SNAPDX, NSNAP,
     &               SNAPDR, TALL, X11, X22, Y11, Y22, AXIST)
               ELSE
                  CALL REFRSH (MP, ML, MS, MR, MSNAP, N, COOR, ILINE,
     &               LTYPE, LCON, NLPS, IFLINE, ILLIST, NSPR, IFSIDE,
     &               ISLIST, LINKP, LINKL, LINKS, LINKR, TITLE, X1, X2,
     &               Y1, Y2, ALPHA, SNAP, SNAPDX, NSNAP, SNAPDR, TALL,
     &               X11, X22, Y11, Y22, AXIST)
               END IF
            END IF
         ELSE IF ((NSNAP(1) .GE. 2) .AND. (NSNAP(2) .GE. 2)) THEN
            SNAP = .TRUE.
            CALL DMESS(DEV1,' SNAP TO GRID POINTS ENABLED')
            IF (.NOT.SNAPDR) THEN
               IF (DRWTAB) THEN
                  CALL REFRSH (MP, ML, MS, MR, MSNAP, N, COOR, ILINE,
     &               LTYPE, LCON, NLPS, IFLINE, ILLIST, NSPR, IFSIDE,
     &               ISLIST, LINKP, LINKL, LINKS, LINKR, TITLE, X1OLD,
     &               X2OLD, Y1OLD, Y2OLD, ALPHA, SNAP, SNAPDX, NSNAP,
     &               SNAPDR, TALL, X11, X22, Y11, Y22, AXIST)
               ELSE
                  CALL REFRSH (MP, ML, MS, MR, MSNAP, N, COOR, ILINE,
     &               LTYPE, LCON, NLPS, IFLINE, ILLIST, NSPR, IFSIDE,
     &               ISLIST, LINKP, LINKL, LINKS, LINKR, TITLE, X1, X2,
     &               Y1, Y2, ALPHA, SNAP, SNAPDX, NSNAP, SNAPDR, TALL,
     &               X11, X22, Y11, Y22, AXIST)
               END IF
            END IF
         END IF

C  INPUT A REGION BY TRYING TO ENCLOSE THE CURRENT LOCATION AS
C  TIGHTLY AS POSSIBLE

      ELSE IF (BUTTON .EQ. '7') THEN
         IF (.NOT.ALPHA) CALL ERASEC (OLDCUR)
         CALL BOXIT (MP, ML, MS, MR, N, IPOINT, COOR, ILINE, LTYPE,
     &      LCON, IREGN, IMAT, NSPR, IFSIDE, ISLIST, LINKP, LINKL,
     &      LINKR, LINKM, NHOLDR, IHOLDR, NHOLDM, IHOLDM, IRGFLG, X, Y,
     &      Y1, Y2, BOXED, MERGE, NOROOM)
         IF (NOROOM) GO TO 170

C  CALCULATE A GOOD LOCATION FOR THE REGION MARKER

         IF ((BOXED) .AND. (.NOT.ALPHA)) THEN
            CALL LTSORT  (MR, LINKR, N(22), II, ADDLNK)
            CALL REGEXT (MP, ML, MS, MR, N, II, COOR, ILINE, LTYPE,
     &         LCON, NLPS, IFLINE, ILLIST, NSPR, IFSIDE, ISLIST, LINKP,
     &         LINKL, LINKS, LINKR, XMIN, XMAX, YMIN, YMAX)
            XMID = (XMAX + XMIN)/2.
            YMID = (YMAX + YMIN)/2.
            CALL PLTSTD (1, 3.)
#if NeedsDoubleEscape
            CALL MPD2SY (1, XMID, YMID, '\\CDI')
#else
            CALL MPD2SY (1, XMID, YMID, '\CDI')
#endif
            CALL PLTSTD (1, 7.)
            CALL DMESS(DEV1,' A REGION CONTAINING THE CURRENT POINT'//
     &         ' HAS BEEN INPUT')
         ELSE
            CALL PLTBEL
            CALL PLTFLU
            CALL DMESS(DEV1,' NO LINES CLOSE ABOUT CURRENT POINT')
         END IF

C  CHANGE THE ZOOM LIMITS WITH THE MOUSE

      ELSE IF (BUTTON .EQ. 'B') THEN

C  SWITCH THE LATEST ZOOM LIMITS

         IF (CHANGE) THEN

C  JUST CHANGE THE SCREEN IF A DRAWING IS ON THE TABLET

            IF(DRWTAB)THEN
               X1HOLD = X1OLD
               Y1HOLD = Y1OLD
               X2HOLD = X2OLD
               Y2HOLD = Y2OLD
               X1OLD = XMIN2
               Y1OLD = YMIN2
               X2OLD = XMAX2
               Y2OLD = YMAX2
               XMIN2 = X1HOLD
               XMAX2 = X2HOLD
               YMIN2 = Y1HOLD
               YMAX2 = Y2HOLD
               CALL REFRSH (MP, ML, MS, MR, MSNAP, N, COOR, ILINE,
     &            LTYPE, LCON, NLPS, IFLINE, ILLIST, NSPR, IFSIDE,
     &            ISLIST, LINKP, LINKL, LINKS, LINKR, TITLE, X1OLD,
     &            X2OLD, Y1OLD, Y2OLD, ALPHA, SNAP, SNAPDX, NSNAP,
     &            SNAPDR, TALL, X11, X22, Y11, Y22, AXIST)

C  OTHERWISE CHANGE THE TABLET AND THE DRAWING

            ELSE
               X1HOLD = X1
               Y1HOLD = Y1
               X2HOLD = X2
               Y2HOLD = Y2
               X1 = XMIN2
               Y1 = YMIN2
               X2 = XMAX2
               Y2 = YMAX2
               XMIN2 = X1HOLD
               XMAX2 = X2HOLD
               YMIN2 = Y1HOLD
               YMAX2 = Y2HOLD
               CALL TABINT (X1, X2, Y1, Y2, CT, ST, SCALE, XX1, YY1,
     &            XX2, YY2, DRWTAB)
               CALL REFRSH (MP, ML, MS, MR, MSNAP, N, COOR, ILINE,
     &            LTYPE, LCON, NLPS, IFLINE, ILLIST, NSPR, IFSIDE,
     &            ISLIST, LINKP, LINKL, LINKS, LINKR, TITLE, X1, X2,
     &            Y1, Y2, ALPHA, SNAP, SNAPDX, NSNAP, SNAPDR, TALL,
     &            X11, X22, Y11, Y22, AXIST)
            ENDIF

            CALL DMESS(DEV1, ' ZOOM HAS BEEN RESET TO THE'//
     &         ' PREVIOUS ZOOM')

C  DRAW THE LOWER LIMITS OF THE ZOOM

         ELSE
            IF (.NOT.ALPHA) THEN
               CALL ERASEC (OLDCUR)
               CALL PLTSTD (1, 3.)
               CALL D2GRID (X, Y, X, Y22)
               CALL D2GRID (X, Y, X22, Y)
               CALL PLTSTD (1, 7.)
               CALL PLTFLU
            END IF
            XHOLD = X
            YHOLD = Y

            CALL DMESS(DEV1,' ENTER "PUCK-B" AT UPPER RIGHT OF '//
     &         'ZOOM - ANY OTHER KEY TO ABANDON')

  140       CONTINUE
            CALL DREAD (X, Y, BUTTON, X1, Y1, XX1, YY1, SCALE, CT, ST)

C  SHOW CURRENT CURSOR LOCATION AS AN AID IN LOCATION OF THE POINT

            IF ((.NOT.ALPHA) .AND. (BUTTON .EQ. '3')) THEN
               CALL PUTCRS (X, Y, OLDCUR)
               GO TO 140

C  INPUT A NEW X, Y FOR THE OTHER EXTREME

            ELSE IF (BUTTON .EQ. 'B') THEN
               IF (.NOT.ALPHA) CALL ERASEC (OLDCUR)
               IF ((X .GT. XHOLD) .AND. (Y .GT. YHOLD)) THEN

C  RESET THE ZOOM LIMITS

                  IF (DRWTAB) THEN
                     XMIN2 = X1OLD
                     XMAX2 = X2OLD
                     YMIN2 = Y1OLD
                     YMAX2 = Y2OLD
                     X1OLD = XHOLD
                     Y1OLD = YHOLD
                     X2OLD = X
                     Y2OLD = Y
                     CALL REFRSH (MP, ML, MS, MR, MSNAP, N, COOR, ILINE,
     &                  LTYPE, LCON, NLPS, IFLINE, ILLIST, NSPR, IFSIDE,
     &                  ISLIST, LINKP, LINKL, LINKS, LINKR, TITLE,
     &                  X1OLD, X2OLD, Y1OLD, Y2OLD, ALPHA, SNAP, SNAPDX,
     &                  NSNAP, SNAPDR, TALL, X11, X22, Y11, Y22, AXIST)
                  ELSE
                     XMIN2 = X1
                     XMAX2 = X2
                     YMIN2 = Y1
                     YMAX2 = Y2
                     X1 = XHOLD
                     Y1 = YHOLD
                     X2 = X
                     Y2 = Y
                     CALL TABINT (X1, X2, Y1, Y2, CT, ST, SCALE, XX1,
     &                  YY1, XX2, YY2, DRWTAB)
                     CALL REFRSH (MP, ML, MS, MR, MSNAP, N, COOR, ILINE,
     &                  LTYPE, LCON, NLPS, IFLINE, ILLIST, NSPR, IFSIDE,
     &                  ISLIST, LINKP, LINKL, LINKS, LINKR, TITLE, X1,
     &                  X2, Y1, Y2, ALPHA, SNAP, SNAPDX, NSNAP, SNAPDR,
     &                  TALL, X11, X22, Y11, Y22, AXIST)
                  ENDIF
                  CALL DMESS(DEV1,' ZOOM HAS BEEN RESET')
               ELSE
                  IF (.NOT.ALPHA) THEN
                     CALL PLTSTD (1, 0.)
                     CALL D2GRID (XHOLD, YHOLD, XHOLD, YDRAW)
                     CALL D2GRID (XHOLD, YHOLD, XDRAW, YHOLD)
                     CALL PLTSTD (1, 7.)
                     CALL PLTBEL
                     CALL PLTFLU
                  END IF
                  CALL DMESS(DEV1,' ZOOM HAS BEEN ABONDONED')
               END IF
            ELSE
               IF (.NOT.ALPHA) THEN
                  CALL PLTSTD (1, 0.)
                  CALL D2GRID (XHOLD, YHOLD, XHOLD, YDRAW)
                  CALL D2GRID (XHOLD, YHOLD, XDRAW, YHOLD)
                  CALL PLTSTD (1, 7.)
                  CALL PLTBEL
                  CALL PLTFLU
               END IF
               WRITE(*,10000)' '//CHAR(27)//'[2J'
            END IF
         END IF

C  DELETE OPTION

      ELSE IF (BUTTON .EQ. 'D') THEN
         CALL DREAD (X, Y, BUTTON, X1, Y1, XX1, YY1, SCALE, CT, ST)

C  DELETE A POINT AND ALL ENTITIES ASSOCIATED WITH THAT POINT

         IF (BUTTON .EQ. '1') THEN
            IF (.NOT.ALPHA) CALL ERASEC (OLDCUR)
            CALL CLOSEP (MP, N(18), X, Y, IPOINT, COOR, LINKP, II)
            IF (.NOT.ALPHA) THEN
               CALL PLTSTT (2, TALL)
               CALL PLTSTD (1, 3.)
               CALL LTSORT (MP, LINKP, II, IPNTR, ADDLNK)
#if NeedsDoubleEscape
               CALL MPD2SY (1, COOR(1, IPNTR), COOR(2, IPNTR), '\\CSQ')
#else
               CALL MPD2SY (1, COOR(1, IPNTR), COOR(2, IPNTR), '\CSQ')
#endif
               CALL PLTSTT (2, TALL*.5)
               CALL PLTFLU
               ISQR = IPNTR

               CALL DMESS (DEV1, ' PRESS "PUCK-1" TO CONFIRM DELETE')
               CALL DREAD (X, Y, BUTTON, X1, Y1, XX1, YY1, SCALE, CT,
     &            ST)

C  ERASE BOXED IN POINT ON THE SCREEN NO MATTER WHAT BUTTON IS ENTERED

               CALL PLTSTT (2, TALL)
               CALL PLTSTD (1, 0.)
#if NeedsDoubleEscape
               CALL MPD2SY (1, COOR(1, ISQR), COOR(2, ISQR), '\\CSQ')
#else
               CALL MPD2SY (1, COOR(1, ISQR), COOR(2, ISQR), '\CSQ')
#endif
               CALL PLTSTD (1, 7.)
               CALL PLTSTT (2, TALL*.5)
               CALL PLTFLU
               ISQR = 0

C  CHECK THE BUTTON TO SEE IF THE POINT IS REALLY TO BE DELETED

               IF (BUTTON .EQ. '1') THEN
                  WRITE(*,10000)' '//CHAR(27)//'[2J'
                  CALL PLTSTD (1, 0.)
                  CALL ERASE (MP, ML, MS, MR, N, COOR, ILINE, LTYPE,
     &               LCON, NLPS, IFLINE, ILLIST, NSPR, IFSIDE, ISLIST,
     &               LINKP, LINKL, LINKS, LINKR, II, ALPHA)
                  LASTP = 0
                  CALL DMESS (DEV1, ' POINT AND ASSOCIATED DATA '//
     &               'DELETED')
               ELSE
                  CALL DMESS (DEV1, ' POINT DELETION ABORTED')
               END IF
            END IF

C  DELETE A LINE AND ALL ENTIES ASSOCIATE WITH IT

         ELSE IF (BUTTON .EQ. '2') THEN
            IF (.NOT.ALPHA) CALL ERASEC (OLDCUR)
            CALL CLOSEL (MP, ML, N, COOR, ILINE, LTYPE, LCON, LINKP,
     &         LINKL, X, Y, BIFIND, IFIND, ADDCEN, XCHOLD, YCHOLD)
            IF (.NOT.ALPHA .AND. BIFIND) THEN
               CALL PLTSTD (1, 3.)
               CALL LTSORT (ML, LINKL, IFIND, IPNTR, ADDLNK)
               LT = LTYPE(IPNTR)
               IP1 = LCON(1, IPNTR)
               IP2 = LCON(2, IPNTR)
               IP3 = LCON(3, IPNTR)
               CALL DLINE (MP, ML, COOR, LINKP, KNUM, LT, IP1, IP2,
     &            IP3, NUMPLT, DUM1, DUM2, TEST, GETMAX, DUM3, DUM4,
     &            DUM5, DUM6)
               CALL PLTFLU

               CALL DMESS (DEV1,'PRESS "PUCK-2" TO CONFIRM DELETE')
               CALL DREAD (X, Y, BUTTON, X1, Y1, XX1, YY1, SCALE, CT,
     &            ST)
               IF (BUTTON .EQ. '2') THEN
                  WRITE(*,10000)' '//CHAR(27)//'[2J'
                  CALL PLTSTD (1, 0.)
                  CALL ERASEL (MP, ML, MS, MR, N, COOR, ILINE, LTYPE,
     &               LCON, NLPS, IFLINE, ILLIST, NSPR, IFSIDE, ISLIST,
     &               LINKP, LINKL, LINKS, LINKR, IFIND, ALPHA)
                  CALL DMESS (DEV1, ' LINE AND ASSOCIATED DATA '//
     &               'DELETED')
               ELSE
                  CALL PLTSTD (1, 7.)
                  CALL DLINE (MP, ML, COOR, LINKP, KNUM, LT, IP1, IP2,
     &               IP3, NUMPLT, DUM1, DUM2, TEST, GETMAX, DUM3, DUM4,
     &               DUM5, DUM6)
                  CALL DMESS (DEV1, ' LINE DELETION ABORTED')
               END IF
               CALL PLTFLU
            ELSE
               CALL DMESS (DEV1, ' LINE NOT FOUND, TRY AGAIN')
            END IF

C  DELETE A GRID LINE

         ELSE IF (BUTTON .EQ. 'A') THEN
            IF (.NOT.ALPHA) CALL ERASEC (OLDCUR)
            IF (SNAP .AND. (NSNAP(1) .GT. 0 .AND. NSNAP(2) .GT. 0))
     &         THEN
               CALL CLOSEG (MSNAP, SNAPDX, NSNAP, X, Y, II, INDEX,
     &            XBOT, XTOP, YBOT, YTOP)
               IF (.NOT.ALPHA) THEN
                  CALL PLTSTD (1, 3.)
                  CALL D2GRID (XBOT, YBOT, XTOP, YTOP)
                  CALL PLTFLU

                  CALL DMESS (DEV1,
     &               'PRESS "PUCK-A" TO CONFIRM DELETE')
                  CALL DREAD (X, Y, BUTTON, X1, Y1, XX1, YY1, SCALE,
     &               CT, ST)
                  IF (BUTTON .EQ. 'A') THEN
                     WRITE(*,10000)' '//CHAR(27)//'[2J'
                     DO 150 I = II, NSNAP(INDEX) - 1
                        SNAPDX(INDEX, I) = SNAPDX(INDEX, I + 1)
  150                CONTINUE
                     NSNAP(INDEX) = NSNAP(INDEX) - 1
                     CALL PLTSTD (1, 0.)
                  ELSE
                     CALL PLTSTD (1, 4.)
                  END IF
                  CALL D2GRID (XBOT, YBOT, XTOP, YTOP)
                  CALL PLTSTD (1, 7.)
                  CALL PLTFLU
               ELSE
                  CALL DMESS(DEV1, ' GRID LINES NOT ENABLED')
               END IF
            END IF

C  DELETE A ZOOM - GO BACK TO THE ORIGINAL DEFAULTS

         ELSE IF (BUTTON .EQ. 'B')THEN
            IF (DRWTAB) THEN
               XMIN2 = X1OLD
               XMAX2 = X2OLD
               YMIN2 = Y1OLD
               YMAX2 = Y2OLD
               X1OLD = XMIN1
               Y1OLD = YMIN1
               X2OLD = XMAX1
               Y2OLD = YMAX1
               CALL REFRSH (MP, ML, MS, MR, MSNAP, N, COOR, ILINE,
     &            LTYPE, LCON, NLPS, IFLINE, ILLIST, NSPR, IFSIDE,
     &            ISLIST, LINKP, LINKL, LINKS, LINKR, TITLE, X1OLD,
     &            X2OLD, Y1OLD, Y2OLD, ALPHA, SNAP, SNAPDX, NSNAP,
     &            SNAPDR, TALL, X11, X22, Y11, Y22, AXIST)
            ELSE
               XMIN2 = X1
               XMAX2 = X2
               YMIN2 = Y1
               YMAX2 = Y2
               X1 = XMIN1
               Y1 = YMIN1
               X2 = XMAX1
               Y2 = YMAX1
               CALL TABINT (X1, X2, Y1, Y2, CT, ST, SCALE, XX1, YY1,
     &            XX2, YY2, DRWTAB)
               CALL REFRSH (MP, ML, MS, MR, MSNAP, N, COOR, ILINE,
     &            LTYPE, LCON, NLPS, IFLINE, ILLIST, NSPR, IFSIDE,
     &            ISLIST, LINKP, LINKL, LINKS, LINKR, TITLE, X1, X2,
     &            Y1, Y2, ALPHA, SNAP, SNAPDX, NSNAP, SNAPDR, TALL,
     &            X11, X22, Y11, Y22, AXIST)
            ENDIF

            CALL DMESS(DEV1,' ZOOM HAS BEEN RESET TO THE DEFAULT')
         END IF

C  EXIT DIGITIZATION - INITIALIZATION OF DIGITIZING PAD IS SAVED

      ELSE IF (BUTTON .EQ. 'E') THEN
         IF (.NOT.ALPHA) THEN
            CALL PLTSTT (2, TALL)
            IF (.NOT.ALPHA) CALL ERASEC (OLDCUR)
            IF ((DEV1 .EQ. 'LS5') .OR. (DEV1 .EQ. 'ls5')) THEN
               WRITE (*, 10000) CHAR(27)//'RA1'
            END IF
            CALL PLTFLU
         END IF

C  ASSUME THE BODY IS TO CONTAIN ALL THE REGIONS ENTERED

         ADDOLD = .TRUE.
         IFOUND = 1
         DO 160 I = NOLD7 + 1, N(7)
            CALL LTSORT (MR, LINKR, IREGN(I), IPNTR, ADDLNK)
            IF (IPNTR .EQ. I) THEN
               CALL INBODY (MR, N(9), IREGN(I), IFOUND, IRPB, ADDOLD,
     &            NOROOM)
               IF (NOROOM) GO TO 170
            END IF
  160    CONTINUE
         RETURN

C  INDICATE THAT THE CLOSEST POINT TO THE CURRENT LOCATION IS NEEDED

      ELSE IF (BUTTON .EQ. '0') THEN
         CHANGE = .TRUE.
         GO TO 110

C  INDICATE THAT DATA IS TO BE INPUT AT THE CLOSEST LOCATION
C  WITH A NEW NUMBERING - FOR USE WITH SLIDE LINE INPUT

      ELSE IF (BUTTON .EQ. 'C') THEN
         SLIDE = .TRUE.
         CHANGE = .TRUE.
         GO TO 110
      END IF

C  GO GET ANOTHER BUTTON AND LOCATION

      GO TO 100

C  MORE ROOM IN DIMENSIONS NEEDED

  170 CONTINUE
      CALL PLTBEL
      CALL PLTFLU
      IF (.NOT.ALPHA) THEN
         CALL PLTSTT (2, TALL)
         IF ((DEV1 .EQ. 'LS5') .OR. (DEV1 .EQ. 'ls5')) THEN
            WRITE (*, 10000) CHAR(27)//'RA1'
         END IF
         CALL PLTFLU
      END IF
      CALL MESAGE (' ')
      CALL PLTBEL
      CALL PLTFLU
      CALL MESAGE ('DIMENSIONS MUST BE INCREASED - PLEASE WAIT')
      CALL PLTBEL
      CALL PLTFLU
      RETURN

10000 FORMAT (1X, A)
      END
