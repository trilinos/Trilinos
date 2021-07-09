C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE RGNSCH (MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN, STEP,
     &   IREGN, IPNTR, N24, MSC, SCHEME, DEFSCH, SCHSTR, LENSCH, NPER,
     &   PENTAG, TRIANG, TRNSIT, HALFC, FILL, ICODE, REMESH)
C***********************************************************************

C     RGNSCH - GET A REGION'S SCHEME

C***********************************************************************

      DIMENSION CIN(MCOM), IIN(MCOM), RIN(MCOM), KIN(MCOM)
      DIMENSION SCHEME(MSC)

      CHARACTER*72 SCHEME, DEFSCH, SCHSTR, CIN

      LOGICAL STEP, PENTAG, TRIANG, TRNSIT, HALFC, FILL, IANS, REMESH

      DATA IEXIT, IOVER, IQUIT /1, 2, 3/

      ICODE = 0

C  CHECK FOR REMESHING

      IF (REMESH) THEN
         SCHSTR = 'X'
      ELSE

C  GET THE INITIAL SCHEME

         IF ((ABS(IREGN) .LE. N24) .AND. (IPNTR .GT. 0)) THEN
            SCHSTR = SCHEME(IPNTR)
         ELSE
            SCHSTR = DEFSCH
         END IF
      END IF
      CALL STRCUT (SCHSTR)
      CALL STRLNG (SCHSTR, LENSCH)

C  STEP PROCESSING

      IF (STEP) THEN
         WRITE (*, 10000) SCHSTR(1:LENSCH)
         CALL INTRUP ('USE CURRENT SCHEME TO BEGIN PROCESSING', IANS,
     &      MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN)

C  CHANGE THE SCHEME

         IF (.NOT.IANS) THEN
  100       CONTINUE
            IF (ICOM .LE. JCOM) THEN
               SCHSTR = CIN(ICOM)
               ICOM = ICOM + 1
               IANS = .TRUE.
            ELSE
               CALL INQSTR ('ENTER PROCESSING SCHEME: ', SCHSTR)
            END IF
            CALL STRCUT (SCHSTR)
            CALL STRLNG (SCHSTR, LENSCH)

C  HELP FOR SCHEME

            IF ((SCHSTR(1:1) .EQ. 'H') .OR.
     &         (SCHSTR(1:1) .EQ. 'h')) THEN
               CALL MESSAGE(' ')
               CALL HELP_FQ (13)
               CALL MESSAGE(' ')
               GO TO 100
            END IF
         END IF

C  BLANK SCHEME

         IF ((LENSCH .LE. 0) .OR. (SCHSTR(1:1) .EQ. ' ')) THEN
            CALL MESSAGE('NO INITIAL SCHEME INPUT')
            CALL MESSAGE('FORCED RECTANGLE PROCESSING USED')
            SCHSTR = ' '
            LENSCH = 1
            GO TO 120
         END IF
      END IF

C  DETERMINE MESHING SCHEME

      PENTAG = .FALSE.
      TRIANG = .FALSE.
      TRNSIT = .FALSE.
      FILL = .FALSE.
      DO 110 J = 1, LENSCH

C  SEE IF A PENTAGON REGION HAS BEEN FLAGGED

         IF ((SCHSTR(J:J) .EQ. 'U') .OR. (SCHSTR(J:J) .EQ. 'u')) THEN
            IF (NPER .GE. 10) THEN
               PENTAG = .TRUE.
               CALL MESSAGE
     &            ('PENTAGON PRIMITIVE REGION PROCESSING USED')
            ELSE
               CALL MESSAGE('PENTAGON REGION GENERATION NOT')
               CALL MESSAGE('POSSIBLE WITH NO. IN PERIMETER < 10')
               CALL MESSAGE('REGULAR PROCESSING WILL BE ATTEMPTED')
            END IF
            GO TO 120

C  SEE IF A TRANSITION REGION HAS BEEN FLAGGED

         ELSE IF ((SCHSTR(J:J) .EQ. 'B') .OR.
     &      (SCHSTR(J:J) .EQ. 'b')) THEN
            IF (NPER .GE. 8) THEN
               TRNSIT = .TRUE.
               HALFC = .FALSE.
               CALL MESSAGE
     &            ('TRANSITION PRIMITIVE REGION PROCESSING USED')
            ELSE
               CALL MESSAGE('TRANSITION REGION GENERATION NOT')
               CALL MESSAGE('POSSIBLE WITH NO. IN PERIMETER < 8')
               CALL MESSAGE('REGULAR PROCESSING WILL BE ATTEMPTED')
            END IF
            GO TO 120

C  SEE IF A SEMI-CIRCLE REGION HAS BEEN FLAGGED

         ELSE IF ((SCHSTR(J:J) .EQ. 'C') .OR.
     &      (SCHSTR(J:J) .EQ. 'c')) THEN
            IF (NPER .GE. 8) THEN
               TRNSIT = .TRUE.
               HALFC = .TRUE.
               CALL MESSAGE
     &            ('SEMICIRCLE PRIMITIVE REGION PROCESSING USED')
            ELSE
               CALL MESSAGE
     &            ('TRANSITION/SEMICIRCLE REGION GENERATION NOT')
               CALL MESSAGE('POSSIBLE WITH NO. IN PERIMETER < 8')
               CALL MESSAGE('REGULAR PROCESSING WILL BE ATTEMPTED')
            END IF
            GO TO 120

C  SEE IF A TRIANGULAR REGION HAS BEEN FLAGGED

         ELSE IF ((SCHSTR(J:J) .EQ. 'T') .OR.
     &      (SCHSTR(J:J) .EQ. 't')) THEN
            IF (NPER .GE. 6) THEN
               TRIANG = .TRUE.
               CALL MESSAGE
     &            ('TRIANGLE PRIMITIVE REGION PROCESSING USED')
            ELSE
               CALL MESSAGE('TRIANGULAR REGION GENERATION NOT')
               CALL MESSAGE('POSSIBLE WITH NO. IN PERIMETER < 6')
               CALL MESSAGE('REGULAR PROCESSING WILL BE ATTEMPTED')
            END IF
            GO TO 120

C  SEE IF A FILL REGION HAS BEEN FLAGGED

         ELSE IF ((SCHSTR(J:J) .EQ. 'X') .OR.
     &      (SCHSTR(J:J) .EQ. 'x')) THEN
            FILL = .TRUE.
            CALL MESSAGE('PAVING TECHNIQUE INITIALLY USED')
            GO TO 120

C  SEE IF A REGULAR RECTANGULAR REGION HAS BEEN FLAGGED

         ELSE IF ((SCHSTR(J:J) .EQ. 'M') .OR.
     &      (SCHSTR(J:J) .EQ. 'm')) THEN
            GO TO 120

C  OTHER POSSIBILITIES

         ELSE IF ((SCHSTR(J:J) .EQ. 'E') .OR.
     &      (SCHSTR(J:J) .EQ. 'e')) THEN
            ICODE = IEXIT
            GO TO 120
         ELSE IF ((SCHSTR(J:J) .EQ. 'O') .OR.
     &      (SCHSTR(J:J) .EQ. 'o')) THEN
            ICODE = IOVER
            GO TO 120
         ELSE IF ((SCHSTR(J:J) .EQ. 'Q') .OR.
     &      (SCHSTR(J:J) .EQ. 'q')) THEN
            ICODE = IQUIT
            GO TO 120
         END IF
  110 CONTINUE
  120 CONTINUE

      RETURN

10000 FORMAT ('0INITIAL MESH DEFINED USING THIS SCHEME:', /, 5X, A)
      END
