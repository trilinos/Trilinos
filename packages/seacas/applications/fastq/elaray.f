C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE ELARAY (XNOLD, YNOLD, NXKOLD, MMPOLD, LINKEG, LISTEG,
     &   MLINK, NPROLD, NPNOLD, NPEOLD, NNXK, XMIN, XMAX, YMIN, YMAX,
     &   IDIVIS)
C***********************************************************************

C  SUBROUTINE ELARAY = PUTS ELEMENTS INTO AN ARRAY BASED ON THEIR
C                      PHYSICAL LOCATION

C***********************************************************************

      DIMENSION XNOLD(NPNOLD), YNOLD(NPNOLD)
      DIMENSION NXKOLD(NNXK, NPEOLD), MMPOLD(3, NPROLD)
      DIMENSION LINKEG(2, MLINK), LISTEG(4 * NPEOLD)

      LOGICAL LCROSS, INSIDE

C  FIND THE EXTREMES FOR THE MESH DATA

      XMIN = XNOLD(1)
      XMAX = XNOLD(1)
      YMIN = YNOLD(1)
      YMAX = YNOLD(1)
      DO 100 I = 2, NPNOLD
         XMIN = AMIN1 (XMIN, XNOLD(I))
         XMAX = AMAX1 (XMAX, XNOLD(I))
         YMIN = AMIN1 (YMIN, YNOLD(I))
         YMAX = AMAX1 (YMAX, YNOLD(I))
  100 CONTINUE

C  SET UP THE SIZE OF THE ARRAY BASED ON THE MLINK DIMENSION
C        IF MLINK = 55 THEN THERE ARE 5 COLUMNS AND 5 ROWS
C                 = 66 THEN THERE ARE 6 COLUMNS AND 6 ROWS, ETC.

      IF (MLINK .EQ. 22) THEN
         IDIVIS = 2
      ELSE IF (MLINK .EQ. 33) THEN
         IDIVIS = 3
      ELSE IF (MLINK .EQ. 44) THEN
         IDIVIS = 4
      ELSE IF (MLINK .EQ. 55) THEN
         IDIVIS = 5
      ELSE IF (MLINK .EQ. 66) THEN
         IDIVIS = 6
      ELSE IF (MLINK .EQ. 77) THEN
         IDIVIS = 7
      ELSE IF (MLINK .EQ. 88) THEN
         IDIVIS = 8
      ELSE IF (MLINK .EQ. 99) THEN
         IDIVIS = 9
      ENDIF

C  NOW THE ELEMENTS MUST BE SORTED INTO ANY ARRAY SPACE THAT THE ELEMENT
C  CROSSES.  THE ARRAY IS LOGICALLY A SQUARE, BUT PHYSICALLY CAN BE
C  RECTANGULAR SINCE THE X AND Y EXTREMES MAY FORM ANY SIZE RECTANGLE.
C  ROWS FIRST IN THE ARRAY AND THEN COLUMNS.

      XDELTA = (XMAX - XMIN) / DBLE(IDIVIS)
      YDELTA = (YMAX - YMIN) / DBLE(IDIVIS)
      KOUNT = 0
      DO 160 J = IDIVIS, 1, -1
         IF (J .EQ. 1) THEN
            YL = YMIN
         ELSE
            YL = YMIN + (YDELTA * DBLE(J - 1))
         ENDIF
         IF (J .EQ. IDIVIS) THEN
            YU = YMAX
         ELSE
            YU = YMIN + (YDELTA * DBLE(J))
         ENDIF
         DO 150 I = 1, IDIVIS
            IF (I .EQ. 1) THEN
               XL = XMIN
            ELSE
               XL = XMIN + (XDELTA * DBLE(I - 1))
            ENDIF
            IF (I .EQ. IDIVIS) THEN
               XU = XMAX
            ELSE
               XU = XMIN + (XDELTA * DBLE(I))
            ENDIF
            INDEX = ((IDIVIS - J + 1) * 10) + I
            LINKEG (1, INDEX) = KOUNT + 1
            LINKEG (2, INDEX) = 0

C  ONLY CHECK ELEMENTS OF THE SAME MATERIAL ID (BLOCK ID)

            DO 140 KELEM = 1, NPEOLD
               DO 120 ICON = 1, 4
                  X1 = XNOLD (NXKOLD (ICON, KELEM))
                  Y1 = YNOLD (NXKOLD (ICON, KELEM))

C  TEST TO SEE IF THE NODE FITS IN THE GRID

                  IF ( ((X1 .LE. XU) .AND. (X1 .GE. XL)) .AND.
     &               ((Y1 .LE. YU) .AND. (Y1 .GE. YL))  ) THEN
                     KOUNT = KOUNT + 1
                     IF (KOUNT .GT. NPEOLD*4) THEN
                        CALL MESAGE ('** ERROR - NOT ENOUGH ROOM '//
     &                     'IN LISTEG, SUBROUTINE ELARAY **')
                        GOTO 170
                     ENDIF
                     LINKEG (2, INDEX) = LINKEG (2, INDEX) + 1
                     LISTEG (KOUNT) = KELEM
                     GOTO 130
                  ENDIF

C  TEST TO SEE IF THE EDGE OF THE ELEMENT CROSSES THE GRID

                  IF (ICON .EQ. 4) THEN
                     JCON = 1
                  ELSE
                     JCON = ICON + 1
                  ENDIF
                  X2 = XNOLD (NXKOLD (JCON, KELEM))
                  Y2 = YNOLD (NXKOLD (JCON, KELEM))
                  CALL INTSCT (X1, Y1, X2, Y2, XL, YL, XU, YL, U, W,
     &               LCROSS)
                  IF (.NOT. LCROSS) CALL INTSCT (X1, Y1, X2, Y2,
     &               XU, YL, XU, YU, U, W, LCROSS)
                  IF (.NOT. LCROSS) CALL INTSCT (X1, Y1, X2, Y2,
     &               XU, YU, XL, YU, U, W, LCROSS)
                  IF (.NOT. LCROSS) CALL INTSCT (X1, Y1, X2, Y2,
     &               XL, YU, XL, YL, U, W, LCROSS)
                  IF (LCROSS) THEN
                     KOUNT = KOUNT + 1
                     IF (KOUNT .GT. NPEOLD*4) THEN
                        CALL MESAGE ('** ERROR - NOT ENOUGH ROOM '//
     &                     'IN LISTEG, SUBROUTINE ELARAY **')
                        GOTO 170
                     ENDIF
                     LINKEG (2, INDEX) = LINKEG (2, INDEX) + 1
                     LISTEG (KOUNT) = KELEM
                     GOTO 130
                  ENDIF

C  OTHERWISE TEST TO SEE IF THE ELEMENT COMPLETELY ENCLOSES THE GRID

                  XEMIN = XNOLD (NXKOLD (1, KELEM))
                  XEMAX = XNOLD (NXKOLD (1, KELEM))
                  YEMIN = YNOLD (NXKOLD (1, KELEM))
                  YEMAX = YNOLD (NXKOLD (1, KELEM))
                  DO 110 IC = 2, 4
                     XEMIN = AMIN1 (XEMIN, XNOLD (NXKOLD (IC, KELEM)))
                     XEMAX = AMAX1 (XEMAX, XNOLD (NXKOLD (IC, KELEM)))
                     YEMIN = AMIN1 (YEMIN, YNOLD (NXKOLD (IC, KELEM)))
                     YEMAX = AMAX1 (YEMAX, YNOLD (NXKOLD (IC, KELEM)))
  110             CONTINUE
                  IF ((XL .GT. XEMIN) .OR. (XU .LT. XEMAX) .OR.
     &               (YL .GT. YEMIN) .OR. (YU .LT. YEMAX)) THEN
                     INSIDE = .FALSE.
                  ELSE
                     CALL INVMAP (X1, Y1, XL, YL, XU, YL, XU, YU, XL,
     &                  YU, XI, ETA, INSIDE)
                  ENDIF
                  IF (INSIDE) THEN
                     KOUNT = KOUNT + 1
                     IF (KOUNT .GT. NPEOLD*4) THEN
                        CALL MESAGE ('** ERROR - NOT ENOUGH ROOM '//
     &                     'IN LISTEG, SUBROUTINE ELARAY **')
                        GOTO 170
                     ENDIF
                     LINKEG (2, INDEX) = LINKEG (2, INDEX) + 1
                     LISTEG (KOUNT) = KELEM
                     GOTO 130
                  ENDIF

  120          CONTINUE
  130          CONTINUE

  140       CONTINUE

  150    CONTINUE

  160 CONTINUE

  170 CONTINUE
      RETURN

      END
