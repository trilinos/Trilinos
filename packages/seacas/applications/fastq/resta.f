C    Copyright(C) 1999-2020, 2024 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE RESTA (MXND, XN, YN, NUID, LXK, KXL, NXL, LXN, KKK,
     &   KKKOLD, NAVAIL, IAVAIL, NNN, LIMIT, IREST, TILT, ERR, NOROOM)
C************************************************************************

C  SUBROUTINE RESTA  =  RESTRUCTURES THE MESH TO ELIMINATE WORST ELEMENTS

C***********************************************************************

C  NOTE:
C     A RECORD IS KEPT OF UP TO 25 OF THE CURRENT WORST CONDITION NUMBERS
C     AND THE WORST ELEMENT POSSIBLE IS RESTRUCTURED
C     UNTIL NO FURTHER RESTRUCTURING CAN BE DONE.

C***********************************************************************

      DIMENSION KCND(26), CND(26)
      DIMENSION NXL(2, 3*MXND), XN(MXND), YN(MXND), NUID(MXND)
      DIMENSION LXK(4, MXND), KXL(2, 3*MXND), LXN(4, MXND)
      DIMENSION NODES(4), ANGLES(4), SIDES(4)

      LOGICAL ERR, NOROOM, LSIDE, CCW, CAREA, DONE

      ERR = .FALSE.

C  CHECK FOR IMPENDING OVERFLOW

      IF (NAVAIL .LE. 1) THEN
         NOROOM = .TRUE.
         CALL MESSAGE('INSUFFICIENT STORAGE AVAILABLE IN RESTA')
         RETURN
      ENDIF

C  INITIALIZE

      NTAB = 0
      MAXTAB = 25
      CNDTOL = 2.0
      ASUM = 0.
      NSUM = 0
      CCW = .TRUE.
      CAREA = .FALSE.
      IREST = 0

      DO 110 K  =  KKKOLD  +  1, KKK
         IF (LXK(1, K) .GT. 0) THEN
            LSIDE = .FALSE.

C  GET THE ELEMENTS COND VALUE (BASED ON ANGLE AND SIDE LENGTH)

            CALL GNXKA (MXND, XN, YN, K, NODES, AREA, LXK, NXL, CCW)
            CALL QAAVAL (MXND, NODES, ANGLES, QRAT, DUMMY, XN, YN,
     &         CAREA)
            CALL CONDNO (MXND, NODES, QRAT, SRAT, COND, SIDES, XN, YN,
     &         LSIDE)

C  ADD UP THE NUMBER OF ANGLES < PI/2

            DO 100 I = 1, 4
               IF (ANGLES(I) .LE. 1.58) THEN
                  ASUM = ASUM + ANGLES(I)
                  NSUM = NSUM + 1
               ENDIF
  100       CONTINUE

C  ADD BAD ELEMENTS TO THE LIST

            IF (COND .GE. CNDTOL) THEN
               CND(NTAB + 1) = COND
               KCND(NTAB + 1) = K
               CALL BUBBLE (CND, KCND, NTAB, NTAB + 1)
               NTAB = MIN0(NTAB + 1, MAXTAB)
            ENDIF
         ENDIF
  110 CONTINUE

C  TILT IS THE AVERAGE VALUE IN DEGREES OF ANGLES < PI/2

      IF (NSUM .GT. 0) THEN
         TILT = (ASUM/DBLE(NSUM))*57.2957795
      ELSE
         TILT = 90.
      ENDIF
      IF ((LIMIT .LE. 0) .OR. (NTAB .LE. 0)) RETURN
      CNDTOL = CND(NTAB)

C  TRY TO RESTRUCTURE ON THE 10 WORST ELEMENTS ONLY

  120 CONTINUE
      NTRY = MIN0(NTAB, 10)
      DO 130 IK = 1, NTRY
         IK1 = IK
         CALL RESTRY (MXND, KCND(IK), K2, LXK, NXL, KXL, LXN, XN, YN,
     &      NUID, NAVAIL, IAVAIL, NNN, DONE, ERR, NOROOM)
         IF (ERR) RETURN
         IF (DONE) GO TO 140
  130 CONTINUE
      RETURN
  140 CONTINUE
      IREST = IREST + 1
      IF (IREST .GE. LIMIT) RETURN

C  UPDATE THE TABLE (AFTER 1 RESTRUCTURE)

      CALL GNXKA (MXND, XN, YN, KCND(IK1), NODES, AREA, LXK, NXL, CCW)
      CALL QAAVAL (MXND, NODES, ANGLES, QRAT, DUMMY, XN, YN, CAREA)
      CALL CONDNO (MXND, NODES, QRAT, SRAT, COND1, SIDES, XN, YN, LSIDE)
      CND(IK1) = COND1
      DO 150 IK = 1, NTAB
         IK2 = IK
         IF (KCND(IK) .EQ. K2) GO TO 160
  150 CONTINUE
      IK2 = NTAB + 1
      NTAB = NTAB + 1
  160 CONTINUE
      CALL GNXKA (MXND, XN, YN, K2, NODES, AREA, LXK, NXL, CCW)
      CALL QAAVAL (MXND, NODES, ANGLES, QRAT, DUMMY, XN, YN, CAREA)
      CALL CONDNO (MXND, NODES, QRAT, SRAT, COND2, SIDES, XN, YN, LSIDE)
      CND(IK2) = COND2
      KCND(IK2) = K2

C  RE-SORT AND PRUNE

      CALL BUBBLE (CND, KCND, 1, NTAB)
      DO 170 I = 1, 2
         IF (CND(NTAB) .LT. CNDTOL)NTAB = NTAB - 1
  170 CONTINUE
      NTAB = MIN0(NTAB, MAXTAB)
      IF (NTAB .LE. 0) RETURN
      CNDTOL = CND(NTAB)
      GO TO 120

      END
