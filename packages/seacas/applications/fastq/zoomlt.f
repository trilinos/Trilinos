C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE ZOOMLT (MCOM, ICOM, JCOM, CIN, RIN, IIN, KIN, IDUMP,
     &   DRAWN, ALPHA, DEV1, X1, X2, Y1, Y2, XX1, XX2, YY1, YY2, XMIN1,
     &   XMAX1, YMIN1, YMAX1, XMIN, XMAX, YMIN, YMAX)
C***********************************************************************

C  ZOOMPL = SUBROUTINE TO INPUT NEW ZOOM LIMITS

C***********************************************************************

      DIMENSION KIN(MCOM), IIN(MCOM), RIN(MCOM)

      CHARACTER*72 CIN(MCOM)
      CHARACTER*3 DEV1, ANS

      LOGICAL DRAWN, ALPHA

      IF ((ICOM .LE. JCOM) .AND. (DRAWN) .AND.
     &   ((CIN(ICOM)(1:1) .EQ. 'C') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'c')) .AND.
     &   (.NOT.ALPHA)) THEN
         CIN(ICOM) = 'PLOT'

C  USE CURSOR INPUT FROM THE SCREEN

         CALL MESAGE (' ')
         CALL MESAGE (' ')
         CALL MESAGE ('LOCATE ONE CORNER WITH CURSOR')
         CALL MESAGE ('THEN HIT ANY KEY')
C         X1 = .45
C         Y1 = .325
         CALL PLTCRS (XTEST1, YTEST1, ANS)
C         XDRAW = ABS ( (XTEST1 * (XX2 - XX1))) + XX1
C         YDRAW = ABS ( ((YTEST1 / .75) * (YY2 - YY1)) ) + YY1
         CALL PLTMOV (0., YTEST1)
         CALL PLTDRW (1., YTEST1)
         CALL PLTMOV (XTEST1, 0.)
         CALL PLTDRW (XTEST1, .75)
         CALL PLTFLU
C         X2 = MAX( X1+.05, .55)
C         Y2 = MAX( Y1+.05, .425)
         CALL PLTCRS (XTEST2, YTEST2, ANS)
         CALL PLTMOV (0., YTEST2)
         CALL PLTDRW (1., YTEST2)
         CALL PLTMOV (XTEST2, 0.)
         CALL PLTDRW (XTEST2, .75)
         CALL MESAGE ('LOCATE THE OTHER CORNER WITH CURSOR')
         CALL MESAGE ('THEN HIT ANY KEY')
         CALL PLTFLU
         X1 = MIN (XTEST1, XTEST2)
         X2 = MAX (XTEST1, XTEST2)
         Y1 = MIN (YTEST1, YTEST2)
         Y2 = MAX (YTEST1, YTEST2)
         XMIN = ABS ( (X1 * (XX2 - XX1))) + XX1
         XMAX = ABS ( (X2 * (XX2 - XX1))) + XX1
         YMIN = ABS ( ((Y1 / .75) * (YY2 - YY1)) ) + YY1
         YMAX = ABS ( ((Y2 / .75) * (YY2 - YY1)) ) + YY1

C  USE USER INPUT FROM THE KEYPAD

      ELSE
         IF ((CIN(ICOM)(1:1) .EQ. 'C') .OR.
     &      (CIN(ICOM)(1:1) .EQ. 'c')) THEN
            ICOM = ICOM+1
            CALL MESAGE (' ')
            CALL MESAGE ('NO CURRENT PLOT FOR CURSOR ZOOM')
            CALL MESAGE ('CURRENT PLOT LIMITS UNCHANGED')
            CALL MESAGE ('* IN OTHER WORDS ... PLOT FIRST (P) '//
     &         'AND THEN ZOOM (Z,C) *')

C  SEE IF ANY OF THE VALUES ARE REDEFINED

         ELSE IF ( (ICOM .LE. JCOM)  .AND.
     &      ( (KIN(ICOM) .GT. 0)    .OR.  (KIN(ICOM+1) .GT. 0)  .OR.
     &      (KIN(ICOM+2) .GT. 0)  .OR.  (KIN(ICOM+3) .GT. 0) ) ) THEN
            IF (KIN(ICOM) .GT. 0) XMIN = RIN(ICOM)
            ICOM = ICOM+1
            IF (ICOM .LE. JCOM) THEN
               IF (KIN(ICOM) .GT. 0) XMAX = RIN(ICOM)
               ICOM = ICOM+1
               IF (ICOM .LE. JCOM) THEN
                  IF (KIN(ICOM) .GT. 0) YMIN = RIN(ICOM)
                  ICOM = ICOM+1
                  IF (ICOM .LE. JCOM) THEN
                     IF (KIN(ICOM) .GT. 0) YMAX = RIN(ICOM)
                     ICOM = ICOM+1
                  END IF
               END IF
            END IF
         ELSE
            XMIN = XMIN1
            YMIN = YMIN1
            XMAX = XMAX1
            YMAX = YMAX1
            CALL MESAGE (' ')
            CALL MESAGE ('ZOOM LIMITS RESET TO PLOT EXTREMES')
         END IF
      END IF

      RETURN

      END
