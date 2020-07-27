C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE RPLOTL (MXND, XN, YN, ZN, NXL, XMIN, XMAX, YMIN, YMAX,
     &   ZMIN, ZMAX, LLL, DEV1, KREG)
C***********************************************************************

C  SUBROUTINE RPLOTL = REPLOTS THE CURRENT MESH FROM THE NXL ARRAY

C***********************************************************************

      DIMENSION NXL (2, 3 * MXND), XN (MXND), YN (MXND), ZN (MXND)
      DIMENSION X (2), Y (2)

      CHARACTER*72 DUMMY, HOLD, DEV1*3

      LOGICAL HARD, FIGURE

      HARD = .FALSE.
      FIGURE = .FALSE.

C  INITIALIZE THE PLOTTING SURFACE

      XDIMD = 1.
      YDIMD = .75

C  TURN ON THE HARDCOPY IF NEEDED

      IF (HARD) THEN
         CALL VDIQES (10002, KAVAL2)
         IF (KAVAL2 .NE. 1) GOTO 110
         CALL VDESCP (10002, 0, 0)
      ENDIF

C  OPEN A FIGURE FILE IF NEEDED

      IF (FIGURE) THEN
         IUNIT = 98
         OPEN (UNIT = IUNIT, FILE = 'DATA.FIG',
     &      STATUS = 'NEW', ERR = 110)
      ENDIF

      CALL PLTBGN
      XDIMR = XMAX - XMIN
      YDIMR = YMAX - YMIN
      CALL MPVIEW (0., XDIMD, 0., YDIMD)
      XRAT = XDIMR/XDIMD
      YRAT = YDIMR/YDIMD
      IF (XRAT.LT.YRAT) THEN
         XDIMR = XDIMD * YRAT
         XX1 =  (XMIN + XMAX - XDIMR) * .5
         XX2 =  (XMIN + XMAX + XDIMR) * .5
         XDIMR = XX2 - XX1
         YY1 = YMIN
         YY2 = YMAX
      ELSE
         YDIMR = YDIMD * XRAT
         YY1 =  (YMIN + YMAX - YDIMR) * .5
         YY2 =  (YMIN + YMAX + YDIMR) * .5
         YDIMR = YY2 - YY1
         XX1 = XMIN
         XX2 = XMAX
      ENDIF
      XX1 = XX1 - (XDIMR * .1)
      XX2 = XX2 + (XDIMR * .1)
      YY1 = YY1 - (YDIMR * .1)
      YY2 = YY2 + (YDIMR * .1)
      CALL MPORT2 (XX1, XX2, YY1, YY2)
      CALL PLTFRM (0)
      CALL GETDUM (KREG, HOLD, LEN)
      DUMMY = ' '
      DUMMY (8:7 + LEN) = HOLD (1:LEN)
      DUMMY (1:7) = 'REGION '
      LEN = LEN + 7
      CALL PLTXTH (XDIMD * .05, YDIMD * .95, DUMMY (1:LEN))

C  PLOT THE LINES IN NXL ARRAY,  SKIPPING DELETIONS

      IF (FIGURE) THEN
         IDUM = 0
         XDUM = 0.
         YDUM = 0.
      ENDIF
      DO 100 I = 1, LLL
         IF (NXL (1, I).GT.0) THEN
            X (2) = XN (NXL (2, I))
            Y (2) = YN (NXL (2, I))
            X (1) = XN (NXL (1, I))
            Y (1) = YN (NXL (1, I))
            CALL MPD2VC (1, X (1), Y (1), X (2), Y (2))
            CALL PLTFLU
            IF ((FIGURE) .AND.
     &         ( ((X (1) .LT. XX2) .AND. (X (1) .GT. XX1) .AND.
     &         (Y (1) .LT. YY2) .AND. (Y (1) .GT. YY1))
     &         .OR.
     &         ((X (2) .LT. XX2) .AND. (X (2) .GT. XX1) .AND.
     &         (Y (2) .LT. YY2) .AND. (Y (2) .GT. YY1)) ) ) THEN
               WRITE (IUNIT, 10000) NXL (1, I) + IDUM, X(1) + XDUM,
     &            Y(1) + YDUM
               WRITE (IUNIT, 10000) NXL (2, I) + IDUM, X(2) + XDUM,
     &            Y(2) + YDUM
               WRITE (IUNIT, 10010) I + IDUM, NXL (1, I) + IDUM,
     &            NXL (2, I) + IDUM
            ENDIF
         ENDIF
  100 CONTINUE

      CALL PLTFLU
      IF (HARD) THEN
         CALL PLTFLU
         CALL VDESCP (10001, 0, 0)
      ENDIF

  110 CONTINUE
      IF (FIGURE) CLOSE (IUNIT)
      RETURN

10000 FORMAT (' POINT ', I6, 2X, 2 (1PE14.7, 2X))
10010 FORMAT (' LINE  ', I6, 2X, 'STR ', I6, 2X, I6)
      END
