C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

      SUBROUTINE GETPRO (NEREPL, NNREPL, *)

      INCLUDE 'g3_xxxxx.blk'
      PARAMETER (MAXFLD = 10)
      CHARACTER*8 WORD, VERB
      INTEGER     INTYP(MAXFLD+1)
      CHARACTER*8 CFIELD(MAXFLD)
      INTEGER     IFIELD(MAXFLD)
      REAL        RFIELD(MAXFLD)

      LOGICAL MATSTR, FFEXST, DOOLD, HELP, ISHELP

      CHARACTER*8 CMDTBL(17)
      SAVE CMDTBL
C      --CMDTBL - the valid commands table

C   --Command table follows.  Remember to change the dimensioned size when
C   --changing the table.
      DATA CMDTBL /
     1   'SCALE   ', 'OFFSET  ', 'SCLCEN  ', 'TWIST   ', 'TWICEN  ',
     2   'END     ', 'EXIT    ', 'NORMAL  ', 'PLANE   ', 'WARP    ',
     3   'RESET   ', 'TORUS   ', 'SPHERE  ', 'XCYLINDE', 'YCYLINDE',
     *   'HELP    ', '        ' /

      CALL SHOCMD ('COMMANDS', CMDTBL)

      XXSCAL = 1.0
      XYSCAL = 1.0
      XXSCL0 = 0.0
      XYSCL0 = 0.0
      XXOFFS = 0.0
      XYOFFS = 0.0
      XZOFFS = 0.0
      XXA    = 0.0
      XXB    = 0.0
      XXC    = 0.0
      XWARP  = 0.0

   10 CONTINUE

C   --Read command line

      WRITE (*, *)
      CALL FREFLD (0, 0, '   PROJECT OPTION> ', MAXFLD,
     &   IOSTAT, NUMFLD, INTYP, CFIELD, IFIELD, RFIELD)
      IF (IOSTAT .LT. 0) GOTO 20
      IF (NUMFLD .EQ. 0) GOTO 10
      INTYP(MIN(MAXFLD,NUMFLD)+1) = -999

      IFLD = 1
      CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
      CALL ABRSTR (VERB, WORD, CMDTBL)
      IF (VERB .EQ. ' ') VERB = WORD

C   --Perform command
      IF (VERB .EQ. 'NORMAL' .OR. VERB .EQ. 'PLANE') THEN
         CALL FFREAL (IFLD, INTYP, RFIELD,
     &      'X normal component', XXA, XXA, *10)
         CALL FFREAL (IFLD, INTYP, RFIELD,
     &      'Y normal component', XXB, XXB, *10)
         CALL FFREAL (IFLD, INTYP, RFIELD,
     &      'Z normal component', XXC, XXC, *10)

         IF (XXC .EQ. 0.0) THEN
            CALL PRTERR ('CMDERR',
     &         '"Z" normal component must be nonzero')
            GO TO 10
         END IF
         RMAG = SQRT(XXA**2 + XXB**2 + XXC**2)
         IF (RMAG .EQ. 0.0) THEN
            CALL PRTERR ('CMDERR',
     *         'Zero length vector entered')
            GO TO 10
         ELSE

C ... NOTE: Since mesh is translated in -Z direction, Z normal to plane
C           must be negative.  If not, then reverse total vector
C           (This was done wrong originally, therefore to not screw
C            up people who figured out a correct orientation, we allow
C            the bug to continue if they enter DOOLDWAY.

            DOOLD = .FALSE.
            IF (FFEXST (IFLD, INTYP)) THEN
               CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
               IF (MATSTR (WORD, 'DOOLDWAY', 4)) THEN
                  DOOLD = .TRUE.
               END IF
            END IF
            IF (DOOLD) THEN
               XXA =  XXA / RMAG
               XXB =  XXB / RMAG
               XXC = -XXC / RMAG
            ELSE IF (XXC .GT. 0) THEN
               XXA = -XXA / RMAG
               XXB = -XXB / RMAG
               XXC = -XXC / RMAG
            ELSE
               XXA =  XXA / RMAG
               XXB =  XXB / RMAG
               XXC =  XXC / RMAG
            END IF

         END IF
         ISXWRP = ISFLAT

      ELSE IF (VERB .EQ. 'WARP' .OR. VERB .EQ. 'SPHERE') THEN
         CONVEX = .TRUE.
         CALL FFREAL (IFLD, INTYP, RFIELD,
     *      'warping distance', 0.0, XWARP, *10)
         IF (FFEXST (IFLD, INTYP)) THEN
            CALL FFCHAR (IFLD, INTYP, CFIELD, 'CONVEX', WORD)
            IF (MATSTR (WORD, 'CONVEX', 4)) THEN
               CONVEX = .TRUE.
            ELSE IF (MATSTR (WORD, 'CONCAVE', 4)) THEN
               CONVEX = .FALSE.
            ELSE
               CALL PRTERR ('CMDWARN', 'unrecognized warp option')
               GO TO 10
            END IF
         END IF
         ISXWRP = ISSPHE

      ELSE IF (VERB .EQ. 'TORUS') THEN
         CONVEX = .TRUE.
         CALL FFREAL (IFLD, INTYP, RFIELD,
     *      'major radius', 0.0, XWARP, *10)
         CALL FFREAL (IFLD, INTYP, RFIELD,
     *      'minor radius', 0.0, YWARP, *10)
         IF (FFEXST (IFLD, INTYP)) THEN
            CALL FFCHAR (IFLD, INTYP, CFIELD, 'CONVEX', WORD)
            IF (MATSTR (WORD, 'CONVEX', 4)) THEN
               CONVEX = .TRUE.
            ELSE IF (MATSTR (WORD, 'CONCAVE', 4)) THEN
               CONVEX = .FALSE.
            ELSE
               CALL PRTERR ('CMDWARN', 'unrecognized option')
               GO TO 10
            END IF
         END IF
         ISXWRP = ISTORO

      ELSE IF (VERB .EQ. 'XCYLINDE') THEN
         CONVEX = .TRUE.
         CALL FFREAL (IFLD, INTYP, RFIELD,
     *      'major radius', 0.0, XWARP, *10)
         IF (FFEXST (IFLD, INTYP)) THEN
            CALL FFCHAR (IFLD, INTYP, CFIELD, 'CONVEX', WORD)
            IF (MATSTR (WORD, 'CONVEX', 4)) THEN
               CONVEX = .TRUE.
            ELSE IF (MATSTR (WORD, 'CONCAVE', 4)) THEN
               CONVEX = .FALSE.
            ELSE
               CALL PRTERR ('CMDWARN', 'unrecognized option')
               GO TO 10
            END IF
         END IF
         ISXWRP = ISXCYL

      ELSE IF (VERB .EQ. 'YCYLINDE') THEN
         CONVEX = .TRUE.
         CALL FFREAL (IFLD, INTYP, RFIELD,
     *      'major radius', 0.0, YWARP, *10)
         IF (FFEXST (IFLD, INTYP)) THEN
            CALL FFCHAR (IFLD, INTYP, CFIELD, 'CONVEX', WORD)
            IF (MATSTR (WORD, 'CONVEX', 4)) THEN
               CONVEX = .TRUE.
            ELSE IF (MATSTR (WORD, 'CONCAVE', 4)) THEN
               CONVEX = .FALSE.
            ELSE
               CALL PRTERR ('CMDWARN', 'unrecognized option')
               GO TO 10
            END IF
         END IF
         ISXWRP = ISYCYL

      ELSE IF (VERB .EQ. 'SCLCEN') THEN
         CALL FFREAL (IFLD, INTYP, RFIELD,
     &      'X scale center', 0.0, XXSCL0, *10)
         CALL FFREAL (IFLD, INTYP, RFIELD,
     &      'Y scale center', 0.0, XYSCL0, *10)

      ELSE IF (VERB .EQ. 'SCALE') THEN
         CALL FFREAL (IFLD, INTYP, RFIELD,
     &      'X scale factor', 0.0, XXSCAL, *10)
         XXSCAL = ABS(XXSCAL)
         CALL FFREAL (IFLD, INTYP, RFIELD,
     &      'Y scale factor', 0.0, XYSCAL, *10)
         XYSCAL = ABS(XYSCAL)

      ELSE IF (VERB .EQ. 'OFFSET') THEN
         CALL FFREAL (IFLD, INTYP, RFIELD,
     &      'X coordinate offset', XXOFFS, XXOFFS, *10)
         CALL FFREAL (IFLD, INTYP, RFIELD,
     &      'Y coordinate offset', XYOFFS, XYOFFS, *10)
         CALL FFREAL (IFLD, INTYP, RFIELD,
     &      'Z coordinate offset', XZOFFS, XZOFFS, *10)

      ELSE IF (VERB .EQ. 'TWIST') THEN
         CALL PRTERR ('CMDSPEC', 'TWIST option not yet implemented')
      ELSE IF (VERB .EQ. 'TWICEN') THEN
         CALL PRTERR ('CMDSPEC', 'TWICEN option not yet implemented')

      ELSE IF (VERB .EQ. 'RESET') THEN
         XXSCAL = 1.0
         XYSCAL = 1.0
         XXSCL0 = 0.0
         XYSCL0 = 0.0
         XXOFFS = 0.0
         XYOFFS = 0.0
         XZOFFS = 0.0
         XXA    = 0.0
         XXB    = 0.0
         XXC    = 0.0
         XWARP  = 0.0
         YWARP  = 0.0
         CALL PRTERR ('CMDSPEC', 'All values have been reset.')

      ELSE IF (VERB .EQ. 'HELP') THEN
         ISHELP = HELP ('GEN3D', 'COMMANDS', CFIELD(IFLD))
         IF (.NOT. ISHELP) CALL SHOCMD ('COMMANDS', CMDTBL)
         VERB = ' '

      ELSE IF (VERB .EQ. 'EXIT'  .OR.  VERB .EQ. 'END') THEN
         IF (XXC .EQ. 0.0 .AND. XWARP .EQ. 0.0) THEN
            CALL PRTERR ('CMDERR',
     &         'either a plane or warp must be entered')
            GO TO 10
         END IF
         RETURN
      ELSE
         CALL SHOCMD ('COMMANDS', CMDTBL)
      END IF

      GO TO 10

   20 CONTINUE
      RETURN 1
      END
