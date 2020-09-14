C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE APARSE (LINE, MAXFLD,NUMFLD, ITYPE, CFIELD, RFIELD)
C=======================================================================

C   --*** APARSE *** (ALGEBRA) Parse line into fields
C   --   Written by Amy Gilkey - revised 02/23/88
C   --
C   --APARSE converts a line to uppercase and parses the line into individual
C   --fields.  A field is one of the following:
C   --   a number (starts with 0-9 or ., rest = 0-9 or . or E or E+ or E-)
C   --   a name (starts with A-Z, rest = A-Z or $ or :)
C   --   a non-number or non-name character; ** is converted to ^
C   --
C   --Parameters:
C   --   LINE - IN - the input line to be parsed, all uppercase
C   --   MAXFLD - IN - the maximum number of parsed fields
C   --   NUMFLD - OUT - the number of parsed fields
C   --   ITYPE - OUT - the field types:
C   --      -1= no value
C   --      0 = CFIELD only valid, name or invalid number string
C   --      1 = CFIELD and RFIELD only valid
C   --      3 = CFIELD only valid, character string
C   --   CFIELD - OUT - the alphanumeric fields
C   --   RFIELD - OUT - the numeric fields

      CHARACTER*(*) LINE
      INTEGER       ITYPE(*)
      CHARACTER*(*) CFIELD(*)
      REAL          RFIELD(*)

      CHARACTER*32 CONVERT
      CHARACTER CH

C   --Initialize fields

      CALL INIINT (MAXFLD, -1, ITYPE)
      CALL INISTR (MAXFLD, ' ', CFIELD)
      CALL INIREA (MAXFLD, 0.0, RFIELD)
      NUMFLD = 0

C   --Change line to upper case

      CALL EXUPCS (LINE)

C   --Repeat until done with input characters

      IEND = LENSTR (LINE)

      NCOL = 1
  100 CONTINUE
      IF (NCOL .LE. IEND) THEN

C      --Skip blank characters

  110    CONTINUE
         IF (LINE(NCOL:NCOL) .EQ. ' ') THEN
            NCOL = NCOL + 1
            GOTO 110
         END IF

         CH = LINE(NCOL:NCOL)

         NUMFLD = NUMFLD + 1

         IF (((CH .GE. '0') .AND. (CH .LE. '9'))
     &      .OR. (CH .EQ. '.')) THEN

C         --Get number string

            ILEFT = NCOL
            IRIGHT = 0
  120       CONTINUE
            IF (IRIGHT .EQ. 0) THEN
               NCOL = NCOL + 1
               CH = LINE(NCOL:NCOL)
               IF (CH .EQ. 'E') THEN
                  CH = LINE(NCOL+1:NCOL+1)
                  IF ((CH .EQ. '+') .OR. (CH .EQ. '-')) NCOL = NCOL + 1
               ELSE IF (CH .EQ. '.') THEN
                  CONTINUE
               ELSE IF ((CH .GE. '0') .AND. (CH .LE. '9')) THEN
                  CONTINUE
               ELSE
                  IRIGHT = NCOL - 1
               END IF
               GOTO 120
            END IF

C         --Convert number string and store

            IF (NUMFLD .LE. MAXFLD) THEN
               CFIELD(NUMFLD) = LINE(ILEFT:IRIGHT)
               ITYPE(NUMFLD) = 0

               CONVERT = ' '
               IJUST = 32 - (IRIGHT - ILEFT + 1) + 1
               CONVERT(IJUST:32) = LINE(ILEFT:IRIGHT)

               READ (CONVERT, '(F32.0)', IOSTAT=ITRANS) RNUM
               IF (ITRANS .EQ. 0) THEN
                  RFIELD(NUMFLD) = RNUM
                  ITYPE(NUMFLD) = 1
               END IF
            END IF

         ELSE IF ((CH .GE. 'A') .AND. (CH .LE. 'Z')) THEN

C         --Get word string and store

            ILEFT = NCOL
            IRIGHT = 0
  130       CONTINUE
            IF (IRIGHT .EQ. 0) THEN
               NCOL = NCOL + 1
               CH = LINE(NCOL:NCOL)
               IF (((CH .GE. 'A') .AND. (CH .LE. 'Z'))
     &            .OR. ((CH .GE. '0') .AND. (CH .LE. '9'))
     &            .OR. (CH .EQ. '$') .OR. (CH .EQ. ':')
     &            .OR. (CH .EQ. '_') .OR. (CH .EQ. '.')) THEN
                  CONTINUE
               ELSE
                  IRIGHT = NCOL - 1
               END IF
               GOTO 130
            END IF

            IF (NUMFLD .LE. MAXFLD) THEN
               CFIELD(NUMFLD) = LINE(ILEFT:IRIGHT)
               ITYPE(NUMFLD) = 0
            END IF

         ELSE

C         --Get single character (or **) and store
C         --Includes { , ( ) = }

            NCOL = NCOL + 1
            IF ((CH .EQ. '*') .AND. (LINE(NCOL:NCOL) .EQ. '*')) THEN
               CH = '^'
               NCOL = NCOL + 1
            END IF

            IF (NUMFLD .LE. MAXFLD) THEN
               CFIELD(NUMFLD) = CH
               ITYPE(NUMFLD) = 3
            END IF
         END IF

         GOTO 100
      END IF

      RETURN
      END
