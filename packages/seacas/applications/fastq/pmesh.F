C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
      SUBROUTINE PMESH (NPNODE, NPELEM, MXNFLG, MXSFLG, NPNBC, NPSBC,
     &   MAXKXN, NNN, KKK, NNXK, NBCNOD, NBCSID, NNLIST, NVLIST, NODES,
     &   NSIDEN, NNFLG, NNPTR, NSFLG, NVPTR, NVLEN, XN, YN, NXK, MAT,
     &   MAPDXG, MAPGXD, WTNODE, WTSIDE, AXIS, AREACG, LABE, LABO, LABN,
     &   LABNB, LABSB, LABM, LABW, CENTK, ILOOK, XMIN, XMAX, YMIN, YMAX,
     &   XX1, XX2, YY1, YY2, TITLE, DEV1, EIGHT, NINE, VERSN, VAXVMS)
C***********************************************************************

C  SUBROUTINE PMESH = PLOTS ALL FLAGGED ELEMENTS IN THE MESH

C***********************************************************************

C  VARIABLES USED:
C     AREACG = .TRUE. IF AREA AND C.G. REPORT ARE TO BE CALCULATED AND
C              PLOTTED
C     LABE   = .TRUE. IF ELEMENT NUMBERS ARE TO BE PLOTTED
C     LABN   = .TRUE. IF NODE NUMBERS ARE TO BE PLOTTED
C     LABNB  = .TRUE. IF NODE BOUNDARY NUMBERS ARE TO BE PLOTTED
C     LABSB  = .TRUE. IF SIDE BOUNDARY NUMBERS ARE TO BE PLOTTED

C***********************************************************************

      DIMENSION CENTK(2, NPELEM), ILOOK(NNXK*MAXKXN)
      DIMENSION MAPGXD(NPNODE), MAPDXG(NPNODE)
      DIMENSION XN(NPNODE), YN(NPNODE), NXK(NNXK, NPELEM), MAT(NPELEM)
      DIMENSION NODES(NPNBC), NSIDEN(NPSBC)
      DIMENSION WTNODE(NPNBC), WTSIDE(NPSBC)
      DIMENSION NNFLG(MXNFLG), NNPTR(MXNFLG)
      DIMENSION NSFLG(MXSFLG), NVPTR(MXSFLG), NVLEN(MXSFLG)
      DIMENSION XDUM(2), YDUM(2)

      LOGICAL LABE, LABN, LABNB, LABSB, LABO, ADD, AXIS, LABM, EIGHT
      LOGICAL NINE, NEW, LABW, CPUIFC, AREACG

      CHARACTER*72 DUMMY, TITLE
      CHARACTER DEV1*3, CDUMMY*4, DATE*8, TIME*8, VERSN*9

C  INITIALIZE THE PLOTTING SURFACE

      CALL PLTBGN
      CALL PLTSTV (2, 160.)
      XDIMR = XMAX - XMIN
      YDIMR = YMAX - YMIN
      XDIMD = 1.
      YDIMD = .75
      CALL MPVIEW (0., XDIMD, 0., YDIMD)
      XRAT = XDIMR/XDIMD
      YRAT = YDIMR/YDIMD
      IF (XRAT .LT. YRAT) THEN
         XDIMR = XDIMD*YRAT
         XX1 = (XMIN + XMAX - XDIMR)*.5
         XX2 = (XMIN + XMAX + XDIMR)*.5
         XDIMR = XX2 - XX1
         YY1 = YMIN
         YY2 = YMAX
      ELSE
         YDIMR = YDIMD*XRAT
         YY1 = (YMIN + YMAX - YDIMR)*.5
         YY2 = (YMIN + YMAX + YDIMR)*.5
         YDIMR = YY2 - YY1
         XX1 = XMIN
         XX2 = XMAX
      END IF

C  SET UP SCALING EXTREMES FOR AXIS

      IF (AXIS) THEN
         XDUM(1) = XX1 - (XDIMR*.05)
         XDUM(2) = XX2 + (XDIMR*.05)
         YDUM(1) = YY1 - (YDIMR*.05)
         YDUM(2) = YY2 + (YDIMR*.05)
         SHRINK = .2
      ELSE
         SHRINK = .1
      END IF

C  SHRINK TO FIT A BORDER ON THE PLOT

      XX1 = XX1 - (XDIMR*SHRINK)
      XX2 = XX2 + (XDIMR*SHRINK)
      YY1 = YY1 - (YDIMR*SHRINK)
      YY2 = YY2 + (YDIMR*SHRINK)
      CALL MPORT2 (XX1, XX2, YY1, YY2)
      CALL PLTFRM (0)

C  PLOT THE TITLE AND THE TRACE

      CALL STRLNG (TITLE, LEN)
      IF ((LEN .GT. 1) .OR. (TITLE(1:1) .NE. ' ')) THEN
         CALL PLTXHL (TITLE(1:LEN), XLEN)
         XBEGIN = AMAX1(0., (XDIMD*.5 - XLEN*.5))
         CALL PLTXTH (XBEGIN, YDIMD*.95, TITLE(1:LEN))
         CALL PLTXHN (X1, YBEG)
      ELSE
         YBEG = YDIMD*.95
      END IF
      DUMMY(1:10) = ' DRAWN BY '
      DUMMY(11:19) = VERSN
      DUMMY(20:21) = '  '
      CALL EXDATE (DATE)
      DUMMY(22:29) = DATE
      DUMMY(30:31) = '  '
      CALL EXTIME (TIME)
      DUMMY(32:39) = TIME
      CALL PLTXTH (0., 0., DUMMY(1:39))

C  DRAW THE AXIS IF REQUIRED AND SET CLIPPING WITHIN THE AXIS

      IF (AXIS) CALL SETAXS (XDUM, YDUM)
      IF (CPUIFC(.TRUE.)) GO TO 340

C  FLAG NODES TO BE PLOTTED

      DO 100 I = 1, NNN
         ILOOK(I) = 0
  100 CONTINUE

      DO 120 I = 1, KKK
         IF (NXK(1, I) .GT. 0) THEN
            DO 110 J = 1, NNXK
               IF (NXK(J, I) .GT. 0) ILOOK(NXK(J, I)) = 1
  110       CONTINUE
         END IF
  120 CONTINUE

C  PLOT THE NUMBERS ASSOCIATED WITH NODES

      IF ((LABN) .OR. (LABNB) .OR. (LABSB) .OR. (LABW)) THEN
         DO 240 I = 1, NNN
            IF (CPUIFC(.TRUE.)) GO TO 340
            IF (ILOOK(I) .EQ. 1) THEN
               ADD = .FALSE.
               DUMMY = ' '
               CALL MP2PT (1, XN(I), YN(I), X1, Y1, MASK)
               IF (MOD(MASK, 2) .NE. 0) THEN

C  PLOT THE NODE NUMBERS

                  IF (LABN) THEN
                     CALL PLTSTD (1, 3.)
                     CALL GETDUM (I, DUMMY, LEN)
                     CALL PLTXTH (X1, Y1, DUMMY(1:LEN))
                     CALL PLTXHE (X1, Y1)
                     ADD = .TRUE.
                  END IF

C  PLOT THE NODAL BOUNDARY CONDITION FLAGS

                  IF (LABNB) THEN
                     CALL PLTSTD (1, 5.)
                     J1 = 1
                     K1 = 1
  130                CONTINUE
                     DO 160 J = J1, NNLIST
                        IF (NODES(J) .EQ. I) THEN
                           DO 140 K = K1, NBCNOD
                              IF (((K .EQ. NBCNOD) .AND.
     &                           (NNPTR(K) .LE. J)) .OR.
     &                           ((NNPTR(K + 1) .GT. J) .AND.
     &                           (NNPTR(K) .LE. J))) THEN
                                 IF (ADD) THEN
                                    CALL PLTXTH (X1, Y1, '/')
                                    CALL PLTXHE (X1, Y1)
                                 END IF
                                 CALL GETDUM (NNFLG(K), DUMMY, LEN)
                                 CALL PLTXTH (X1, Y1, DUMMY(1:LEN))
                                 CALL PLTXHE (X1, Y1)
                                 ADD = .TRUE.
                                 IF (K .LT. NBCNOD) THEN
                                    K1 = K + 1
                                    J1 = NNPTR(K1)
                                    GO TO 130
                                 ELSE
                                    GO TO 150
                                 END IF
                              END IF
  140                      CONTINUE
  150                      CONTINUE
                        END IF
  160                CONTINUE
                  END IF

C  PLOT THE WEIGHTING FACTORS

                  IF (LABW) THEN
                     CALL PLTSTD (1, 5.)
                     DO 170 J = 1, NNLIST
                        IF (NODES(J) .EQ. I) THEN
                           IF (ADD) THEN
                              CALL PLTXTH (X1, Y1, '/')
                              CALL PLTXHE (X1, Y1)
                           END IF
                           CALL GTXDUM (WTNODE(J), DUMMY, LEN)
                           CALL PLTXTH (X1, Y1, DUMMY(1:LEN))
                           CALL PLTXHE (X1, Y1)
                           ADD = .TRUE.
                        END IF
  170                CONTINUE
                     CALL PLTSTD (1, 7.)
                     DO 200 J = 1, NBCSID
                        DO 180 K = NVPTR(J), NVPTR(J) + NVLEN(J) - 1
                           IF (NSIDEN(K) .EQ. I) THEN
                              IF (ADD) THEN
                                 CALL PLTXTH (X1, Y1, '/')
                                 CALL PLTXHE (X1, Y1)
                              END IF
                              CALL GTXDUM (WTSIDE(K), DUMMY, LEN)
                              CALL PLTXTH (X1, Y1, DUMMY(1:LEN))
                              CALL PLTXHE (X1, Y1)
                              ADD = .TRUE.
                              GO TO 190
                           END IF
  180                   CONTINUE
  190                   CONTINUE
  200                CONTINUE
                  END IF

C  PLOT THE SIDE BOUNDARY CONDITION FLAGS

                  IF (LABSB) THEN
                     J1 = 1
                     ADD = .FALSE.
                     NEW = .TRUE.
                     CALL PLTSTD (1, 7.)
                     K1 = 1
                     JE = 0
                     J1HOLD = 0
  210                CONTINUE
                     DO 230 J = J1, NVLIST, 2
                        IF ((NSIDEN(J) .EQ. I) .AND.
     &                     (ILOOK(NSIDEN(J + 1)) .GT. 0)) THEN

C  SEE IF THIS USE OF THE NODE IS THE SAME SIDE AS BEFORE

                           IF ((ADD) .AND. (NSIDEN(J + 1) .NE. JE) .AND.
     &                        (.NOT.NEW)) THEN
                              J1HOLD = J
                           ELSE
                              JE = NSIDEN(J + 1)
                              DO 220 K = K1, NBCSID
                                 IF (((K .EQ. NBCSID) .AND.
     &                              (NVPTR(K) .LE. J)) .OR.
     &                              ((NVPTR(K + 1) .GT. J) .AND.
     &                              (NVPTR(K) .LE. J))) THEN
                                    NEW = .FALSE.
                                    K1 = K
                                    J1 = J + 2
                                    IF (ADD) THEN
                                       CALL PLTXTH (X1, Y1, '/')
                                       CALL PLTXHE (X1, Y1)
                                    ELSE
                                       X2 = .5*(XN(NSIDEN(J)) +
     &                                    XN(NSIDEN(J + 1)))
                                       Y2 = .5*(YN(NSIDEN(J)) +
     &                                    YN(NSIDEN(J + 1)))
                                       CALL MP2PT (1, X2, Y2, X1, Y1,
     &                                    MASK)
                                       ADD = .TRUE.
                                    END IF
                                    CALL GETDUM (NSFLG(K), DUMMY, LEN)
                                    CALL PLTXTH (X1, Y1, DUMMY(1:LEN))
                                    CALL PLTXHE (X1, Y1)
                                    GO TO 210
                                 END IF
  220                         CONTINUE
                           END IF
                        END IF
  230                CONTINUE

C  GO BACK AND PICK UP OTHER SIDES USING THIS NODE

                     IF (J1HOLD .GT. 0) THEN
                        J1 = J1HOLD
                        J1HOLD = 0
                        NEW = .TRUE.
                        ADD = .FALSE.
                        GO TO 210
                     END IF
                  END IF

               END IF
            END IF
  240    CONTINUE
      END IF

C  PLOT THE ELEMENT NUMBERS IF NEEDED

      IF ((LABE) .OR. (LABO) .OR. (LABM)) THEN
         DO 250 IJK = 1, KKK
            IF (CPUIFC(.TRUE.)) GO TO 340
            I = MAPDXG(IJK)
            IF (NXK(1, I) .GT. 0) THEN
               CALL MP2PT (1, CENTK(1, I), CENTK(2, I), X1, Y1, MASK)
               IF (MOD(MASK, 2) .NE. 0) THEN
                  IF (LABE) THEN
                     CALL PLTSTD (1, 2.)
                     CALL GETDUM (I, DUMMY, LEN)
                     CALL PLTXTH (X1, Y1, DUMMY(1:LEN))
                  ELSE IF (LABO) THEN
                     CALL PLTSTD (1, 2.)
                     CALL GETDUM (IJK, DUMMY, LEN)
                     CALL PLTXTH (X1, Y1, DUMMY(1:LEN))
                  END IF
                  IF (LABM) THEN
                     CALL PLTSTD (1, 1.)
                     CALL GETDUM (MAT(I), DUMMY, LEN)
                     IF ((LABE) .OR. (LABO)) THEN
                        CALL PLTXHE (X1, Y1)
                        CALL PLTXTH (X1, Y1, '/')
                        CALL PLTXHE (X1, Y1)
                     END IF
                     CALL PLTXTH (X1, Y1, DUMMY(1:LEN))
                     CALL PLTSTD (1, 2.)
                  END IF
               END IF
            END IF
  250    CONTINUE
      END IF

C  PLOT THE ELEMENTS

      DO 260 I = 1, MAXKXN*NNXK
         ILOOK(I) = 0
  260 CONTINUE
      JCOLOR = -10000
      AREA = 0.0
      XCG = 0.0
      YCG = 0.0
      DO 320 IJK = 1, KKK
         IF (CPUIFC(.TRUE.)) GO TO 340
         I = MAPDXG(IJK)
         IF (NXK(1, I) .GT. 0) THEN
            COLOR = DBLE(MOD(MAT(I) + 1, 6) + 1)
            CALL PLTSTD (1, COLOR)
            IF (NXK(3, I) .EQ. 0) THEN
               JEND = 1
            ELSE IF (NXK(4, I) .EQ. 0) THEN
               JEND = 2
            ELSE IF ((EIGHT) .OR. (NINE)) THEN
               JEND = 8
            ELSE
               JEND = 4
            END IF
            IF (AREACG) THEN
               EAREA = 0.0
               XBAR = 0.0
               YBAR = 0.0
            END IF
            DO 310 J = 1, JEND
               INUM1 = NXK(J, I)
               K = J + 1
               IF (((K .EQ. 9) .AND. (EIGHT .OR. NINE)) .OR.
     &            ((K .EQ. 5) .AND. (.NOT.EIGHT) .AND. (.NOT.NINE)))
     &            K = 1
               INUM2 = NXK(K, I)

C  CALCULATE AREA AND C.G. OF ELEMENT

               IF (AREACG) THEN
                  XBAR = XBAR + XN(INUM1)
                  YBAR = YBAR + YN(INUM1)
                  EAREA = EAREA + 0.5 * (XN(INUM1) - XN(INUM2))
     &               * (YN(INUM1) + YN(INUM2))
               END IF

               IF (INUM2 .GT. INUM1) THEN
                  IHOLD = INUM1
                  INUM1 = INUM2
                  INUM2 = IHOLD
               END IF
               HOLD = ABS((XN(INUM1)*XN(INUM2))) + (INUM1*INUM2*.013927)
               DO 270 IRED = 1, 100
                  IF (HOLD .GT. 1.E7) HOLD = HOLD*1.E-6
                  IF (HOLD .LT. 1.E7) GOTO 280
  270          CONTINUE
  280          CONTINUE
               LOOKUP = INT((HOLD - INT(HOLD))*DBLE(MAXKXN*NNXK))
               IF (LOOKUP .LE. 0) LOOKUP = 1
               IHOLD = LOOKUP
  290          CONTINUE
               IF ((JEND .EQ. 1) .OR. (ILOOK(LOOKUP) .EQ. 0)) THEN
                  CALL MPD2VC (1, XN(INUM1), YN(INUM1), XN(INUM2),
     &               YN(INUM2))
                  LOOKUP = IHOLD
  300             CONTINUE
                  IF (ILOOK(LOOKUP) .LE. 0) THEN
                     ILOOK(LOOKUP) = INUM1*10000 + INUM2
                     KLOOK = 0
                  ELSE
                     LOOKUP = LOOKUP + 1
                     GO TO 300
                  END IF
               ELSE IF (ILOOK(LOOKUP) .EQ. (INUM1*10000 + INUM2)) THEN
                  KLOOK = 0
                  ILOOK(LOOKUP) = -1
               ELSE
                  LOOKUP = LOOKUP + 1
                  IF (LOOKUP .GT. MAXKXN*NNXK) LOOKUP = 1
                  KLOOK = KLOOK + 1
                  IF (KLOOK .GT. MAXKXN*NNXK) THEN
                     CALL MESAGE ('LOOKUP TABLE TOO FULL')
                     GO TO 330
                  END IF
                  GO TO 290
               END IF
  310       CONTINUE

            IF (AREACG) THEN
               AREA = AREA + EAREA
               XCG = XCG + EAREA*XBAR/REAL(JEND)
               YCG = YCG + EAREA*YBAR/REAL(JEND)
            END IF
         END IF
  320 CONTINUE
      IF (AREACG) THEN
         XCG = XCG/AREA
         YCG = YCG/AREA
      END IF
  330 CONTINUE

C  CONCLUDE THE PLOT SESSION

  340 CONTINUE
      CALL PLTSTD (1, 7.)

C  PLOT SYMBOL AT C.G. AND REPORT AREA AND C.G.

      IF (AREACG) THEN
         WRITE (DUMMY, 10040) 'AREA', AREA
         CALL STRLNG (DUMMY, LEN)
         CALL PLTXTH (0., YBEG, DUMMY(1:LEN))
         WRITE (DUMMY, 10040) 'XBAR', XCG
         CALL STRLNG (DUMMY, LEN)
         CALL PLTXHN (X1, Y1)
         SPACE = (YBEG - Y1) * 0.2
         CALL PLTXTH (0., Y1 - SPACE, DUMMY(1:LEN))
         WRITE (DUMMY, 10040) 'YBAR', YCG
         CALL STRLNG (DUMMY, LEN)
         CALL PLTXHN (X1, Y1)
         CALL PLTXTH (0., Y1 - SPACE, DUMMY(1:LEN))
         CALL MP2PT (1, XCG, YCG, X1, Y1, MASK)
         IF (MOD(MASK, 2) .NE. 0) THEN
#if NeedsDoubleEscape
            CDUMMY = '\\CCS'
#else
            CDUMMY = '\CCS'
#endif
            CALL PLTXTS (X1, Y1, CDUMMY)
#if NeedsDoubleEscape
            CDUMMY = '\\CCI'
#else
            CDUMMY = '\CCI'
#endif
            CALL PLTXTS (X1, Y1, CDUMMY)
         END IF
      END IF
      CALL PLTBEL
      CALL PLTFLU

      RETURN

10040 FORMAT (A, ': ', G11.4)
      END
