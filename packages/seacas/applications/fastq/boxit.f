C $Id: boxit.f,v 1.4 2004/01/26 17:28:18 gdsjaar Exp $
C $Log: boxit.f,v $
C Revision 1.4  2004/01/26 17:28:18  gdsjaar
C Removed several unused variables from getang subroutine.
C
C Initialized a variable
C
C Revision 1.3  1998/07/14 18:18:25  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.2  1991/03/21 15:44:19  gdsjaar
C Changed all 3.14159... to atan2(0.0, -1.0)
C
c Revision 1.1.1.1  1990/11/30  11:04:09  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:04:07  gdsjaar
c Initial revision
c 
C
CC* FILE: [.MAIN]BOXIT.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE BOXIT (MP, ML, MS, MR, N, IPOINT, COOR, ILINE, LTYPE,
     &   LCON, IREGN, IMAT, NSPR, IFSIDE, ISLIST, LINKP, LINKL, LINKR,
     &   LINKM, NHOLDR, IHOLDR, NHOLDM, IHOLDM, IRGFLG, X, Y, Y1, Y2,
     &   BOXED, MERGE, NOROOM)
C***********************************************************************
C
C  SUBROUTINE BOXIT = BOXES IN A REGION SURROUNDING A POINT
C
C***********************************************************************
C
C  SUBROUTINE CALLED BY:
C     INPUT  = INPUTS MESH DEFINITIONS FROM THE LIGHT TABLE
C
C***********************************************************************
C
C  SUBROUTINES CALLED:
C     DLPARA = DETERMINES LINE PARAMETERS FROM TWO POINTS
C
C***********************************************************************
C
      PARAMETER (MHOLD = 50)
C
      DIMENSION IPOINT (MP)
      DIMENSION COOR (2, MP), ILINE (ML), LTYPE (ML), LCON (3, ML)
      DIMENSION IREGN (MR), IMAT (MR), NSPR (MR), IFSIDE (MR)
      DIMENSION ISLIST (MR * 4)
      DIMENSION LINKP (2, MP), LINKL (2, ML), LINKR (2, MR)
      DIMENSION LINKM (2, (MS + MR))
      DIMENSION IHOLDR (2, MR), IHOLDM (2, (MS + MR))
      DIMENSION IHOLD (MHOLD, 2), JHOLD (MHOLD), IRGFLG (MR), N (29)
C
      LOGICAL BOXED, NOROOM, ADDLNK, MERGE, ERR, CIRCLE
C
      PI = ATAN2(0.0, -1.0)
      TWOPI = PI + PI
      BOXED = .FALSE.
      NOROOM = .FALSE.
      ADDLNK = .FALSE.
      CIRCLE = .FALSE.
      IFIND = 0
      SUMTH = 0.
      THETMX = 0.0
C
C  FIND THE CLOSEST LINE ABOVE THE POINT INPUT
C
      CALL LABOVE  (MP, ML, N, IPOINT, COOR, ILINE, LTYPE, LCON, LINKP,
     &   LINKL, X, Y, Y1, Y2, IFIND, JFIND1, ISTART, NP)
      IF (IFIND .LE. 0) RETURN
      IFHOLD = IFIND
C
C  SET UP REGION CHECKING CONNECTIVITY
C
      DO 100 I = 1, MHOLD
         JHOLD (I) = 0
         IHOLD (I, 1) = 0
         IHOLD (I, 2) = 0
  100 CONTINUE
      IFOUND = 0
      LASTP = ISTART
C
      CALL LTSORT (ML, LINKL, IFIND, JFIND1, ADDLNK)
      JHOLD (1) = IFIND
C
      DO 130 I = 1, N (2) + 2
         JKOUNT = 0
C
C  GET ALL LINES CONTAINING THE NEW POINT "NP"
C
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
C
C  CHECK FOR A CLOSED CIRCLE WITH NO LINES ATTACHED
C
         IF  ( (JKOUNT .EQ. 0) .AND.  (NP .EQ. LASTP) ) THEN
            IFOUND = 1
            GOTO 140
C
C  CHECK FOR NO ADDITIONAL LINES ATTACHED
C
         ELSEIF  (JKOUNT .EQ. 0) THEN
            RETURN
C
C  CHECK FOR A CLOSED CIRCLE ATTACHED TO NP
C
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
C
C  CHECK FOR CLOSING OF THE REGION
C
            IF  ( IFIND .EQ. IFHOLD ) THEN
               IF  (  (LASTP .EQ. ISTART) .AND.
     &            (IFIND .EQ. IFHOLD) .AND.
     &            (I  .NE.  1) ) THEN
                  GOTO 140
               ENDIF
            ENDIF
C
C  SET THE FLAG THAT WE ARE RETURNING FROM THIS CLOSED CIRCLE
C
            CIRCLE = .TRUE.
C
C  USING THE NP COORDINATES AS A NEW CENTER
C  CHECK TO SEE WHICH LINE HAS THE SMALLEST INTERIOR ANGLE
C   (ASSUMES WE ARE PROGRESSING IN CLOCKWISE ORDER)
C
         ELSE
            CALL ENDTAN (MP, ML, N, IPOINT, COOR, LTYPE, LCON, LINKP,
     &         LINKL, IFIND, JFIND1, NP, THETA1, ERR)
C
C  SET THE SPECIAL CASE OF A CLOSED CIRCLE WITH ATTACHED LINES
C  AND NOT A RETURN FROM A CIRCLE AS HAVING A THETMX OF PI.
C  A CLOSED CIRCLE RETURN GETS THE ENDTANGENT FLIPPED TO GET THE
C  RIGHT INTERIOR ANGLE.
C
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
C
            DO 120 J = 1, JKOUNT
C
C  JJ = THE POINTER TO THE ATTACHED POINT OF THE LINE BEING TESTED
C
               CALL ENDTAN  (MP, ML, N, IPOINT, COOR, LTYPE, LCON,
     &            LINKP, LINKL, ILINE (IHOLD (J, 1)), IHOLD (J, 1), NP,
     &            THETA2, ERR)
C
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
C
C  CHECK FOR CLOSING OF THE REGION
C
            IF (IFIND .EQ. IFHOLD) THEN
C
C  FIRST THE SINGLE LINE REGION WITH LINES ATTACHED OUTSIDE THE CIRCLE
C
               IF  (  (NP .EQ. LASTP) .AND.
     &            ( IFIND .EQ. ILINE ( JFIND1 ) ) ) THEN
                  IFOUND = 1
                  GOTO 140
C
C  SECOND TEST FOR THE NORMAL CLOSING
C
               ELSEIF  (  (LASTP .EQ. ISTART) .AND.
     &            (IFIND .EQ. IFHOLD) .AND.
     &            (I .NE. 1) ) THEN
                  GOTO 140
               ENDIF
            ENDIF
C
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
C
C  CHECK TO MAKE SURE THE REGION CLOSED CORRECTLY
C
C      AVETH = SUMTH / FLOAT (IFOUND)
C      IF (AVETH .GT. 180.) THEN
C         CALL VDBELL
C         CALL VDBUFL
C         RETURN
C      ENDIF
C
C  INPUT THE REGION
C
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
