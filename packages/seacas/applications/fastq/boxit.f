C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE BOXIT (MP, ML, MS, MR, N, IPOINT, COOR, ILINE, LTYPE,
     &   LCON, IREGN, IMAT, NSPR, IFSIDE, ISLIST, LINKP, LINKL, LINKR,
     &   LINKM, NHOLDR, IHOLDR, NHOLDM, IHOLDM, IRGFLG, X, Y, Y1, Y2,
     &   BOXED, MERGE, NOROOM)
C***********************************************************************

C  SUBROUTINE BOXIT = BOXES IN A REGION SURROUNDING A POINT

C***********************************************************************

C  SUBROUTINE CALLED BY:
C     INPUT  = INPUTS MESH DEFINITIONS FROM THE LIGHT TABLE

C***********************************************************************

C  SUBROUTINES CALLED:
C     DLPARA = DETERMINES LINE PARAMETERS FROM TWO POINTS

C***********************************************************************

      PARAMETER (MHOLD = 50)

      DIMENSION IPOINT (MP)
      DIMENSION COOR (2, MP), ILINE (ML), LTYPE (ML), LCON (3, ML)
      DIMENSION IREGN (MR), IMAT (MR), NSPR (MR), IFSIDE (MR)
      DIMENSION ISLIST (MR * 4)
      DIMENSION LINKP (2, MP), LINKL (2, ML), LINKR (2, MR)
      DIMENSION LINKM (2, (MS + MR))
      DIMENSION IHOLDR (2, MR), IHOLDM (2, (MS + MR))
      DIMENSION IHOLD (MHOLD, 2), JHOLD (MHOLD), IRGFLG (MR), N (29)

      LOGICAL BOXED, NOROOM, ADDLNK, MERGE, ERR, CIRCLE

      PI = ATAN2(0.0, -1.0)
      TWOPI = PI + PI
      BOXED = .FALSE.
      NOROOM = .FALSE.
      ADDLNK = .FALSE.
      CIRCLE = .FALSE.
      IFIND = 0
      SUMTH = 0.
      THETMX = 0.0

C  FIND THE CLOSEST LINE ABOVE THE POINT INPUT

      CALL LABOVE  (MP, ML, N, IPOINT, COOR, ILINE, LTYPE, LCON, LINKP,
     &   LINKL, X, Y, Y1, Y2, IFIND, JFIND1, ISTART, NP)
      IF (IFIND .LE. 0) RETURN
      IFHOLD = IFIND

C  SET UP REGION CHECKING CONNECTIVITY

      DO 100 I = 1, MHOLD
         JHOLD (I) = 0
         IHOLD (I, 1) = 0
         IHOLD (I, 2) = 0
  100 CONTINUE
      IFOUND = 0
      LASTP = ISTART

      CALL LTSORT (ML, LINKL, IFIND, JFIND1, ADDLNK)
      JHOLD (1) = IFIND

      DO 130 I = 1, N (2) + 2
         JKOUNT = 0

C  GET ALL LINES CONTAINING THE NEW POINT "NP"

         DO 110 J = 1, N (19)
            CALL LTSORT (ML, LINKL, J, JJ, ADDLNK)
            IF (JJ .GT. 0) THEN
               J1 = LCON (1, JJ)
               J2 = LCON (2, JJ)
               IF ( ( (J1 .EQ. NP) .OR. (J2 .EQ. NP)) .AND.
     &            (JJ .NE. JFIND1)) THEN
                  JKOUNT = JKOUNT + 1
                  IHOLD (JKOUNT, 1) = JJ
                  IF (J1 .EQ. NP) THEN
                     IHOLD (JKOUNT, 2) = 2
                  ELSE
                     IHOLD (JKOUNT, 2) = 1
                  ENDIF
               ENDIF
            ENDIF
  110    CONTINUE

C  CHECK FOR A CLOSED CIRCLE WITH NO LINES ATTACHED

         IF  ( (JKOUNT .EQ. 0) .AND.  (NP .EQ. LASTP) ) THEN
            IFOUND = 1
            GOTO 140

C  CHECK FOR NO ADDITIONAL LINES ATTACHED

         ELSEIF  (JKOUNT .EQ. 0) THEN
            RETURN

C  CHECK FOR A CLOSED CIRCLE ATTACHED TO NP

         ELSEIF  (  (JKOUNT .EQ. 1) .AND.
     &      (LCON (1, IHOLD (1, 1)) .EQ. LCON (2, IHOLD (1, 1))) )
     &      THEN
            JFIND1 = IHOLD (1, 1)
            SUMTH = SUMTH + THETMX
            LASTP = NP
            NP = LCON (1, JFIND1)
            IFIND = ILINE (JFIND1)
            JHOLD (I + 1) = ILINE (JFIND1)
            IFOUND = IFOUND + 1

C  CHECK FOR CLOSING OF THE REGION

            IF  ( IFIND .EQ. IFHOLD ) THEN
               IF  (  (LASTP .EQ. ISTART) .AND.
     &            (IFIND .EQ. IFHOLD) .AND.
     &            (I  .NE.  1) ) THEN
                  GOTO 140
               ENDIF
            ENDIF

C  SET THE FLAG THAT WE ARE RETURNING FROM THIS CLOSED CIRCLE

            CIRCLE = .TRUE.

C  USING THE NP COORDINATES AS A NEW CENTER
C  CHECK TO SEE WHICH LINE HAS THE SMALLEST INTERIOR ANGLE
C   (ASSUMES WE ARE PROGRESSING IN CLOCKWISE ORDER)

         ELSE
            CALL ENDTAN (MP, ML, N, IPOINT, COOR, LTYPE, LCON, LINKP,
     &         LINKL, IFIND, JFIND1, NP, THETA1, ERR)

C  SET THE SPECIAL CASE OF A CLOSED CIRCLE WITH ATTACHED LINES
C  AND NOT A RETURN FROM A CIRCLE AS HAVING A THETMX OF PI.
C  A CLOSED CIRCLE RETURN GETS THE ENDTANGENT FLIPPED TO GET THE
C  RIGHT INTERIOR ANGLE.

            IF  ( NP .EQ. LASTP ) THEN
               IF  (CIRCLE) THEN
                  THETA1 = THETA1 - PI
                  IF (THETA1 .LT. 0.) THEN
                     THETA1 = THETA1 + TWOPI
                  ELSEIF (THETA1 .GT. TWOPI) THEN
                     THETA1 = THETA1 - TWOPI
                  ENDIF
                  THETMX = TWOPI * 2.
               ELSE
                  THETMX = PI
               ENDIF
            ELSE
               THETMX = TWOPI * 2.
            ENDIF

            DO 120 J = 1, JKOUNT

C  JJ = THE POINTER TO THE ATTACHED POINT OF THE LINE BEING TESTED

               CALL ENDTAN  (MP, ML, N, IPOINT, COOR, LTYPE, LCON,
     &            LINKP, LINKL, ILINE (IHOLD (J, 1)), IHOLD (J, 1), NP,
     &            THETA2, ERR)

               TESTTH = THETA2 - THETA1
               IF (TESTTH .LT. 0) THEN
                  TESTTH = TESTTH + TWOPI
               ELSEIF (TESTTH .GT. TWOPI) THEN
                  TESTTH = TESTTH - TWOPI
               ENDIF
               IF (  (.NOT.ERR) .AND.
     &            (  (TESTTH .LE. THETMX) .OR.
     &            (ABS (THETMX - TESTTH) .LT. .0001) ) ) THEN
                  THETMX = TESTTH
                  JFIND1 = IHOLD (J, 1)
                  JFIND2 = IHOLD (J, 2)
               ENDIF
  120       CONTINUE

C  CHECK FOR CLOSING OF THE REGION

            IF (IFIND .EQ. IFHOLD) THEN

C  FIRST THE SINGLE LINE REGION WITH LINES ATTACHED OUTSIDE THE CIRCLE

               IF  (  (NP .EQ. LASTP) .AND.
     &            ( IFIND .EQ. ILINE ( JFIND1 ) ) ) THEN
                  IFOUND = 1
                  GOTO 140

C  SECOND TEST FOR THE NORMAL CLOSING

               ELSEIF  (  (LASTP .EQ. ISTART) .AND.
     &            (IFIND .EQ. IFHOLD) .AND.
     &            (I .NE. 1) ) THEN
                  GOTO 140
               ENDIF
            ENDIF

            SUMTH = SUMTH + THETMX
            LASTP = NP
            NP = LCON (JFIND2, JFIND1)
            IFIND = ILINE (JFIND1)
            JHOLD (I + 1) = ILINE (JFIND1)
            IFOUND = IFOUND + 1
            CIRCLE = .FALSE.
         ENDIF
  130 CONTINUE
      RETURN

C  CHECK TO MAKE SURE THE REGION CLOSED CORRECTLY

C      AVETH = SUMTH / DBLE(IFOUND)
C      IF (AVETH .GT. 180.) THEN
C         CALL VDBELL
C         CALL VDBUFL
C         RETURN
C      ENDIF

C  INPUT THE REGION

  140 CONTINUE
      JJ = N (22) + 1
      IMTRL = 0
      CALL IREVER (JHOLD, IFOUND)
      DO 150 J = 1, IFOUND
         JHOLD (J) =  - JHOLD (J)
  150 CONTINUE
      CALL INREGN (MS, MR, N (7), N (8), N (22), N (23), JJ, IMTRL,
     &   JHOLD, IFOUND, IREGN, IMAT, NSPR, IFSIDE, ISLIST, LINKR,
     &   LINKM, NHOLDR, IHOLDR, NHOLDM, IHOLDM, IRGFLG, MERGE, NOROOM)
      IF (NOROOM)RETURN
      BOXED = .TRUE.
      RETURN
      END
