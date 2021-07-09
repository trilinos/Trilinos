C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE FFVRNG (IFLD, INTYP, CFIELD, IFIELD, EXPECT, MAXVAL,
     &   IRANGE, *)
C=======================================================================

C   --*** FFVRNG *** (FFLIB) Parse free-field integer range
C   --   Written by Amy Gilkey - revised 02/24/86
C   --
C   --FFVRNG parses a range of integers.  A range has one of the following
C   --forms:
C   --            n1                  assume n2 = n1, n3 = 1
C   --            n1 TO n2            assume n3 = 1
C   --            n1 TO n2 STEP n3
C   --
C   --Parameters:
C   --   IFLD - IN/OUT - the index of the current field number, incremented
C   --   INTYP - IN - the input type array from the free-field reader
C   --   CFIELD - IN - the input string array from the free-field reader
C   --   IFIELD - IN - the input integer array from the free-field reader
C   --   EXPECT - IN - the type of range being parsed, for error
C   --   MAXVAL - IN - the maximum range value (ignore if < 0)
C   --   IRANGE - OUT - the input range value array:
C   --          (1) = n1, (2) = n2, (3) = n3;
C   --      partially set on error
C   --   * - return statement if the range is invalid; message is printed

C   --Routines Called:
C   --   LENSTR - (STRLIB) Find string length

      INTEGER IFLD
      INTEGER INTYP(*)
      CHARACTER*(*) CFIELD(*)
      INTEGER IFIELD(*)
      CHARACTER*(*) EXPECT
      INTEGER MAXVAL
      INTEGER IRANGE(3)

      CHARACTER*80 STRA
      CHARACTER*80 ERRMSG

      IRANGE(1) = 0
      IRANGE(2) = 0
      IRANGE(3) = 1

      IF (INTYP(IFLD) .GE. -1) THEN

C      --Get starting number

         IF (INTYP(IFLD) .NE. 2) THEN
            WRITE (ERRMSG, 10000)
     &         EXPECT, CFIELD(IFLD)(:LENSTR(CFIELD(IFLD)))
            CALL PRTERR ('CMDERR', ERRMSG(:LENSTR(ERRMSG)))
            GOTO 100
         END IF
         IRANGE(1) = IFIELD(IFLD)
         IRANGE(2) = IRANGE(1)
         IRANGE(3) = 1
         IFLD = IFLD + 1

         IF (INTYP(IFLD) .EQ. 0) THEN

C         --Get TO and ending value

            IF (CFIELD(IFLD) .NE. 'TO') THEN
               WRITE (ERRMSG, 10000)
     &            '"TO"', CFIELD(IFLD)(:LENSTR(CFIELD(IFLD)))
               CALL PRTERR ('CMDERR', ERRMSG(:LENSTR(ERRMSG)))
               GOTO 100
            END IF
            IFLD = IFLD + 1

            IF (INTYP(IFLD) .NE. 2) THEN
               STRA = 'TO ' // EXPECT
               WRITE (ERRMSG, 10000) STRA(:LENSTR(STRA)),
     &            CFIELD(IFLD)(:LENSTR(CFIELD(IFLD)))
               CALL PRTERR ('CMDERR', ERRMSG(:LENSTR(ERRMSG)))
               GOTO 100
            END IF
            IRANGE(2) = IFIELD(IFLD)
            IFLD = IFLD + 1

            IF (INTYP(IFLD) .EQ. 0) THEN

C            --Get BY and step value

               IF (CFIELD(IFLD) .NE. 'BY') THEN
                  WRITE (ERRMSG, 10000)
     &               '"BY"', CFIELD(IFLD)(:LENSTR(CFIELD(IFLD)))
                  CALL PRTERR ('CMDERR', ERRMSG(:LENSTR(ERRMSG)))
                  GOTO 100
               END IF
               IFLD = IFLD + 1

               IF (INTYP(IFLD) .NE. 2) THEN
                  WRITE (ERRMSG, 10000)
     &               'BY value', CFIELD(IFLD)(:LENSTR(CFIELD(IFLD)))
                  CALL PRTERR ('CMDERR', ERRMSG(:LENSTR(ERRMSG)))
                  GOTO 100
               END IF
               IRANGE(3) = IFIELD(IFLD)
               IFLD = IFLD + 1
            END IF
         END IF

         IF (IRANGE(3) .EQ. 0) THEN
            WRITE (ERRMSG, 10010, IOSTAT=IDUM)
     &         'Invalid BY value', IRANGE(3)
            CALL SQZSTR (ERRMSG, L)
            CALL PRTERR ('CMDERR', ERRMSG(:L))
            GOTO 110
         END IF
         IF ((IRANGE(3) .GT. 0) .AND. (IRANGE(1) .GT. IRANGE(2))) THEN
            STRA = 'Starting ' // EXPECT
            WRITE (ERRMSG, 10010, IOSTAT=IDUM) STRA(:LENSTR(STRA)),
     &         IRANGE(1), ' > ending ', IRANGE(2)
            CALL SQZSTR (ERRMSG, L)
            CALL PRTERR ('CMDERR', ERRMSG(:L))
            GOTO 110
         END IF
         IF ((IRANGE(3) .LT. 0) .AND. (IRANGE(1) .LT. IRANGE(2))) THEN
            STRA = 'Starting ' // EXPECT
            WRITE (ERRMSG, 10010, IOSTAT=IDUM) STRA(:LENSTR(STRA)),
     &         IRANGE(1), ' < ending ', IRANGE(2), ' with negative step'
            CALL SQZSTR (ERRMSG, L)
            CALL PRTERR ('CMDERR', ERRMSG(:L))
         END IF
         IF (MIN (IRANGE(1), IRANGE(2)) .GT. MAXVAL .AND.
     *     MAXVAL .GT. 0) THEN
            STRA = 'Minimum ' // EXPECT
            WRITE (ERRMSG, 10010, IOSTAT=IDUM) STRA(:LENSTR(STRA)),
     &         MIN (IRANGE(1), IRANGE(2)), ' > maximum ', MAXVAL
            CALL SQZSTR (ERRMSG, L)
            CALL PRTERR ('CMDERR', ERRMSG(:L))
            GOTO 110
         END IF
         IF ((IRANGE(1) .LE. 0) .OR. (IRANGE(2) .LE. 0)) THEN
            STRA = 'Negative or zero ' // EXPECT
            WRITE (ERRMSG, 10010, IOSTAT=IDUM) STRA(:LENSTR(STRA)),
     &         MIN (IRANGE(1), IRANGE(2)), ' > maximum ', MAXVAL
            CALL SQZSTR (ERRMSG, L)
            CALL PRTERR ('CMDERR', ERRMSG(:L))
            GOTO 110
         END IF

         IF (IRANGE(1) .GT. MAXVAL .AND. MAXVAL .GT. 0)
     *     IRANGE(1) = MAXVAL
         IF (IRANGE(2) .GT. MAXVAL .AND. MAXVAL .GT. 0)
     *     IRANGE(2) = MAXVAL
      END IF

      RETURN

  100 CONTINUE
      IF (INTYP(IFLD) .GE. -1) IFLD = IFLD + 1
  110 CONTINUE
      RETURN 1
10000  FORMAT ('Expected ', A, ', not "', A, '"')
10010  FORMAT (A, I10, A, I10, A)
      END
