C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE MIDNOD (NPNODE, NNUID, NPELEM, NNXK, MP, ML, KKK, NNN,
     &   NALL, NL, NXK, NUID, XN, YN, LISTN, COOR, ILINE, LTYPE, LCON,
     &   LINKP, LINKL, THREE, EIGHT, NINE)
C***********************************************************************

C  SUBROUTINE MIDNOD = GENERATES THE MIDSIDE NODE FOR ELEMENTS

C***********************************************************************

      DIMENSION XN (NPNODE), YN (NPNODE)
      DIMENSION NUID (NNUID), NXK (NNXK, NPELEM)
      DIMENSION LISTN (NNUID)
      DIMENSION COOR (2, MP), ILINE (ML), LTYPE (ML), LCON (3, ML)
      DIMENSION LINKP (2, MP), LINKL (2, ML)

      LOGICAL ADDLNK, THREE, EIGHT, NINE, ITSOK

      PI = ATAN2(0.0, -1.0)
      TUPI = 2.0 * PI

      ADDLNK = .FALSE.

      NALL = NNN
      DO 130 J = 1, KKK
         DO 120 I = 1, 4

C  SKIP DUPLICATE DESCRIPTORS AND ELEMENTS THAT ARE NOT TO
C  HAVE MIDSIDE NODES

            IF (NXK (I, J) .LT. 0) THEN

C  THIS ELEMENT IS A BAR ELEMENT

               IF (NXK (3, J) .EQ. 0) THEN
                  IF (THREE) THEN
                     ITSOK = .TRUE.
                  ELSE
                     ITSOK = .FALSE.
                  ENDIF

C  THIS ELEMENT IS A QUAD ELEMENT

               ELSEIF ((EIGHT) .OR. (NINE)) THEN
                  ITSOK = .TRUE.
               ELSE
                  ITSOK = .FALSE.
               ENDIF
            ELSE
               ITSOK = .FALSE.
            ENDIF

C  GENERATE THE MIDSIDE NODE IF APPROPRIATE

            IF (ITSOK) THEN
               II = I + 1
               IF (I .GE. 4)II = 1
               NODEA = IABS (NXK (I, J))
               NODEB = IABS (NXK (II, J))

C  CHECK TO SEE IF THE ELEMENT IS A BARSET

               IF ((NODEA .GT. 0) .AND. (NODEB .GT. 0)) THEN
                  NUIDA = NUID (NODEA)
                  NUIDB = NUID (NODEB)

C  IF ONE NODE OR THE OTHER IS INTERIOR,  USE LINEAR INTERPOLATION

                  IF ( ((NUIDA .GT. 100000) .AND.
     &               (NUIDA .LT. 1000000000)) .OR.
     &               ((NUIDB .GT. 100000) .AND.
     &               (NUIDB .LT. 1000000000)) ) THEN
                     XINT = 0.5 * (XN (NODEA) + XN (NODEB))
                     YINT = 0.5 * (YN (NODEA) + YN (NODEB))

C  BOTH ARE POINT OR LINE NODES.
C  FIND THE LINE THAT THEY BELONG TO

                  ELSE
                     LTEST = 0

C  SEE IF ONE IS NOT THE END POINT OF A LINE
C  THEN CHECK THE OTHER TO SEE IF IT IS ON THE SAME LINE
C  IF IT IS NOT ON THE SAME LINE,  THEN USE LINEAR INTERPOLATION

                     IF (NUIDA .GT. 1000000000) THEN
                        LTEST =  (NUIDA - 1000000000) / 100000
                        CALL LTSORT (ML, LINKL, LTEST, LT, ADDLNK)
                        IF (NUIDB .GT. 1000000000) THEN
                           LNEW =  (ABS (NUIDB) - 1000000000) / 100000
                           IF (LNEW.NE.LTEST)LTEST = 0
                        ELSE
                           IF ((LCON (1, LT) .NE. NUIDB) .AND.
     &                        (LCON (2, LT) .NE. NUIDB)) LTEST = 0
                        ENDIF
                     ELSEIF (NUIDB .GT. 1000000000) THEN
                        LTEST =  (NUIDB - 1000000000) / 100000
                        CALL LTSORT (ML, LINKL, LTEST, LT, ADDLNK)
                        IF (NUIDA .GT. 1000000000) THEN
                           LNEW =  (ABS (NUIDA) - 1000000000) / 100000
                           IF (LNEW .NE. LTEST) LTEST = 0
                        ELSE
                           IF ((LCON (1, LT) .NE. NUIDA) .AND.
     &                        (LCON (2, LT) .NE. NUIDA)) LTEST = 0
                        ENDIF
                     ELSE

C  BOTH ARE END POINTS OF LINES - SEE IF THEY ARE ON THE SAME LINE

                        NSUM = ABS (NUIDA) + ABS (NUIDB)
                        DO 100 L = 1, NL
                           CALL LTSORT (ML, LINKL, ILINE (L), K, ADDLNK)
                           IF ((LCON (1, K) + LCON (2, K)) .EQ. NSUM)
     &                        THEN
                              IF ( ((LCON (1, K) .EQ. ABS (NUIDA)) .AND.
     &                           (LCON (2, K) .EQ. ABS (NUIDB))) .OR.
     &                           ((LCON (1, K) .EQ. ABS (NUIDB)) .AND.
     &                           (LCON (2, K) .EQ. ABS (NUIDA)))) THEN
                                 LTEST = ILINE (L)
                                 GOTO 110
                              ENDIF
                           ENDIF
  100                   CONTINUE
  110                   CONTINUE
                     ENDIF

C  IF THEY ARE NOT ON THE SAME LINE,  IT IS NOT IN ERROR
C  THEY SPAN A SINGLE WIDTH REGION
C  ASSUME LINEAR INTERPOLATION

                     IF (LTEST .EQ. 0) THEN
                        KT = 0
                     ELSE
                        CALL LTSORT (ML, LINKL, LTEST, LT, ADDLNK)
                        KT = LTYPE (LT)
                     ENDIF

C  CALCULATE THE MID-SIDE NODE BY LINEAR INTERPOLATION
C  IF THE LINE TYPE IS STRAIGHT OR CORNER OR THE 2 END POINTS
C  SPAN A SINGLE WIDTH REGION

                     IF (KT .LT. 3) THEN
                        XINT = 0.5 * (XN (NODEA) + XN (NODEB))
                        YINT = 0.5 * (YN (NODEA) + YN (NODEB))

C  IF THE LINE IS A CIRCLE OR PARABOLA,  GENERATE THE MID-POINT BY USING
C  THE 2 NODES AS ARC ENDPOINTS AND AN INTERVAL OF 2 FOR THE LINE.
C  THE CENTER OF THE LINE MUST BE FOUND FROM THE LINE CARD ITSELF.

                     ELSE
                        CALL LTSORT (MP, LINKP, LCON (1, LT), IP1,
     &                     ADDLNK)
                        CALL LTSORT (MP, LINKP, LCON (2, LT), IP2,
     &                     ADDLNK)
                        IF (LCON (3, LT) .GT. 0) THEN
                           CALL LTSORT (MP, LINKP, LCON (3, LT), IP3,
     &                        ADDLNK)
                        ELSEIF (LCON (3, LT) .LT. 0) THEN
                           CALL LTSORT (MP, LINKP, IABS (LCON (3, LT)),
     &                        IP3, ADDLNK)
                           IP3 = - IP3
                        ELSE
                           IP3 = 0
                        ENDIF

C  ARC WITH CENTER GIVEN
C  ARC GOES FROM 1ST POINT TO 2ND IN *COUNTER-CLOCKWISE* DIRECTION.

                        IF (KT .EQ. 3) THEN
                           XCEN = COOR (1, IABS (IP3))
                           YCEN = COOR (2, IABS (IP3))

C  CIRCLE WITH THIRD POINT ON ARC.

                        ELSEIF (KT .EQ. 4) THEN
                           THETA1 = ATAN2 (COOR (2, IP3) -
     &                        COOR (2, IP1),  COOR (1, IP3) -
     &                        COOR (1, IP1)) + PI / 2.0
                           THETA2 = ATAN2 (COOR (2, IP3) -
     &                        COOR (2, IP2),  COOR (1, IP3) -
     &                        COOR (1, IP2)) + PI / 2.0
                           DET = - COS (THETA1) * SIN (THETA2) +
     &                        COS (THETA2) * SIN (THETA1)
                           X1 = 0.5 * (COOR (1, IP1) + COOR (1, IP3))
                           Y1 = 0.5 * (COOR (2, IP1) + COOR (2, IP3))
                           X2 = 0.5 * (COOR (1, IP2) + COOR (1, IP3))
                           Y2 = 0.5 * (COOR (2, IP2) + COOR (2, IP3))
                           R =  ( - SIN (THETA2) * (X2 - X1) +
     &                        COS (THETA2) * (Y2 - Y1)) / DET
                           XCEN = X1 + R * COS (THETA1)
                           YCEN = Y1 + R * SIN (THETA1)

C  PARABOLA WITH CENTER GIVEN

                        ELSEIF (KT .EQ. 5) THEN
                           XCEN = COOR (1, IABS (IP3))
                           YCEN = COOR (2, IABS (IP3))

C  CIRCLE WITH RADIUS GIVEN

                        ELSEIF (KT .EQ. 6) THEN
                           DX = 0.5 * (COOR (1, IP2) - COOR (1, IP1))
                           DY = 0.5 * (COOR (2, IP2) - COOR (2, IP1))
                           CHORD = SQRT (DX * DX + DY * DY)
                           R = ABS (COOR (1, IABS (IP3)))
                           IF (R.LE.CHORD) THEN
                              XCEN = 0.5 * (COOR (1, IP1) +
     &                           COOR (1, IP2))
                              YCEN = 0.5 * (COOR (2, IP1) +
     &                           COOR (2, IP2))
                           ELSE
                              ARM = SQRT (R * R - CHORD * CHORD)
                           ENDIF
                           IF (IP3 .LT. 0) THEN
                              XCEN = COOR (1, IP1) + DX +
     &                           ARM * DY / CHORD
                              YCEN = COOR (2, IP1) + DY -
     &                           ARM * DX / CHORD
                           ELSE
                              XCEN = COOR (1, IP1) + DX -
     &                           ARM * DY / CHORD
                              YCEN = COOR (2, IP1) + DY +
     &                           ARM * DX / CHORD
                           ENDIF
                        ENDIF

C  CALCULATE THE MIDPOINT ON THE ARC

                        R1 = SQRT ((XN (NODEA) - XCEN) **2 +
     &                     (YN (NODEA) - YCEN) **2)
                        R2 = SQRT ((XN (NODEB) - XCEN) **2 +
     &                     (YN (NODEB) - YCEN) **2)
                        RM =  (R1 + R2) * .5
                        THETA1 = ATAN2 (YN (NODEA) - YCEN,
     &                     XN (NODEA) - XCEN)
                        IF (THETA1 .LT. 0)THETA1 = THETA1 + TUPI
                        THETA2 = ATAN2 (YN (NODEB) - YCEN,
     &                     XN (NODEB) - XCEN)
                        IF (THETA2 .LT. 0)THETA2 = THETA2 + TUPI
                        THETAM =  (THETA1 + THETA2) * .5
                        IF (ABS (THETA1 - THETA2) .GT. PI)
     &                     THETAM = THETAM + PI
                        XINT = COS (THETAM) * RM + XCEN
                        YINT = SIN (THETAM) * RM + YCEN
                     ENDIF
                  ENDIF

C  ADD THIS NEW NODE TO THE NODE LIST.

                  IF (NALL .GE. NPNODE) THEN
                     WRITE (*, 10000)NPNODE
                     GOTO 140
                  ENDIF
                  NALL = NALL + 1
                  NLO = MIN0 (NODEA, NODEB)
                  NHI = MAX0 (NODEA, NODEB)
                  LISTN (NALL) = NLO * 100000 + NHI
                  XN (NALL) = XINT
                  YN (NALL) = YINT
               ENDIF
            ENDIF
  120    CONTINUE
  130 CONTINUE

  140 CONTINUE
      RETURN

10000 FORMAT (' NODE ARRAY OVERFLOW IN MIDNOD', /
     &   ' THERE ARE MORE THAN', I5, ' NODES IN THE MESH')

      END
