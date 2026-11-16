C=======================================================================
      PROGRAM FORMTEX
C=======================================================================

C                       Written by Amy Gilkey
C                       Last Revised 06/01/88

C   --FORMTEX formats a TEX file.

      PARAMETER (MAXCOL = 72)

      CHARACTER*120 LASTLINE
      CHARACTER*120 THISLINE
      CHARACTER*20 WORD, WORD2
      LOGICAL DONELAST
      LOGICAL DONETHIS
      LOGICAL VERBATIM
      LOGICAL ASIS

      INTEGER*2 I2

C   --Open the input and output files

      NTEX = 20

      OPEN (NTEX, FILE='INFILE', READONLY, STATUS='OLD', IOSTAT=IERR)
      IF (IERR .NE. 0) THEN
         PRINT *, 'ERROR - Cannot open input file'
         GOTO 120
      END IF

      NFMT = 30

      OPEN (NFMT, FILE='OUTFILE', CARRIAGECONTROL='LIST',
     &   STATUS='NEW')

C   --Process the file

      LASTEND = -1
      VERBATIM = .FALSE.

      DO WHILE (.TRUE.)
         READ (NTEX, '(A)', IOSTAT=IOSTAT) THISLINE
         IDUM = STR$TRIM (THISLINE, THISLINE, I2)
         IEND = JZEXT (I2)
         DONELAST = .FALSE.
         DONETHIS = .FALSE.
         ASIS = VERBATIM

         IF (IOSTAT .LT. 0) THEN
            DONELAST = .TRUE.
            IEND = -1

         ELSE IF (IEND .EQ. 0) THEN
            DONELAST = .TRUE.
            DONETHIS = .TRUE.

         ELSE IF (THISLINE(1:1) .EQ. '%') THEN
            DONELAST = .TRUE.
            DONETHIS = .TRUE.
            ASIS = .TRUE.

         ELSE IF (THISLINE(1:1) .EQ. ' ') THEN
            DONELAST = .TRUE.

         ELSE IF ((THISLINE(1:1) .EQ. '}')
     &      .OR. (THISLINE(1:1) .EQ. '{')) THEN
            DONELAST = .TRUE.
            IF (IEND .EQ. 1) DONETHIS = .TRUE.

         ELSE IF (THISLINE(1:1) .EQ. '&') THEN
            DONELAST = .TRUE.

         ELSE IF (THISLINE(1:1) .EQ. '\') THEN
            IL = 1
            CALL SCANSLASH (THISLINE, IL, IEND, WORD)

            IF ((WORD .EQ. '\caps')
     &         .OR. (WORD .EQ. '\cmd')
     &         .OR. (WORD .EQ. '\param')
     &         .OR. (WORD .EQ. '\bold')
     &         .OR. (WORD .EQ. '\mbox')
     &         .OR. (WORD .EQ. '\ref')
     &         .OR. (WORD .EQ. '\ldots')
     &         .OR. (WORD .EQ. '\PROGRAM')
     &         .OR. (WORD .EQ. '\\')) THEN
               CONTINUE
            ELSE IF (WORD .EQ. '\item') THEN
               DONELAST = .TRUE.
               IF (IL .GE. IEND) DONETHIS = .TRUE.
            ELSE
               DONELAST = .TRUE.
               DONETHIS = .TRUE.
               IF ((WORD .EQ. '\begin') .OR. (WORD .EQ. '\end')) THEN
                  CALL SCANPARAM (THISLINE, IL, IEND, WORD2)
                  IF (WORD2 .EQ. 'verbatim') THEN
                     VERBATIM = (WORD .EQ. '\begin')
                     IF (.NOT. VERBATIM) ASIS = .FALSE.
                  END IF
               END IF
            END IF
         END IF

         IF (IEND .GT. 0) THEN
            IF (THISLINE(IEND:IEND) .EQ. '{') THEN
               DONETHIS = .TRUE.
            END IF
            IF (THISLINE(IEND:IEND) .EQ. '\') THEN
               DONETHIS = .TRUE.
            END IF
         END IF

         IF (ASIS) THEN
            DONELAST = .TRUE.
            DONETHIS = .TRUE.
         END IF

         IF (DONELAST) THEN
            IF (LASTEND .GE. 0) THEN
               WRITE (NFMT, 10000) LASTLINE(:LASTEND)
               LASTEND = -1
            END IF
         END IF

         IF (IOSTAT .LT. 0) GOTO 110

         IF (.NOT. ASIS) THEN
            DO WHILE (.TRUE.)
               I = INDEX (THISLINE, '.  ')
               IF ((I .LE. 0) .OR. (I .GE. IEND)) GOTO 100
               THISLINE(I+1:) = THISLINE(I+2:)
               IEND = IEND - 1
            END DO
  100       CONTINUE
         END IF

         DO WHILE ((.NOT. ASIS) .AND. (IEND .GT. 0))
            IF (LASTEND .LT. 0) LASTEND = 0
            IF (LASTEND .EQ. 0) THEN
               NCH = 0
            ELSE
               NCH = LASTEND + 1
            END IF
            IF ((NCH + IEND) .LE. MAXCOL) THEN
               IF (NCH .GT. 0) LASTLINE(NCH:NCH) = ' '
               LASTLINE(NCH+1:) = THISLINE(1:IEND)
               LASTEND = NCH + IEND
               IEND = -1
            ELSE
               ISTART = LIB$SKPC (' ', THISLINE)
               I = MAXCOL - NCH + 1
               IF (I .LT. ISTART) I = 0
               DO WHILE ((I .GT. 0) .AND. (THISLINE(I:I) .NE. ' '))
                  I = I - 1
               END DO
               IF ((LASTEND .EQ. 0) .AND. (I .LE. 0))
     &            I = INDEX (THISLINE, ' ')
               IF (I .GT. 0) THEN
                  IF (NCH .GT. 0) LASTLINE(NCH:NCH) = ' '
                  LASTLINE(NCH+1:) = THISLINE(1:I-1)
C#????                  LASTEND = NCH + I-1
                  IDUM = STR$TRIM (LASTLINE, LASTLINE, I2)
                  LASTEND = JZEXT (I2)
                  THISLINE = THISLINE(I+1:IEND)
                  IEND = IEND - I
               END IF
               IF (LASTEND .GE. 0) THEN
                  WRITE (NFMT, 10000) LASTLINE(:LASTEND)
                  LASTEND = -1
               END IF
            END IF
         END DO

         IF (DONETHIS) THEN
            IF (LASTEND .GE. 0) THEN
               WRITE (NFMT, 10000) LASTLINE(:LASTEND)
               LASTEND = -1
            END IF
            IF (IEND .GE. 0) THEN
               WRITE (NFMT, 10000) THISLINE(:IEND)
               IEND = -1
            END IF
         END IF
      END DO

  110 CONTINUE

      CLOSE (NTEX)
      CLOSE (NFMT)

  120 CONTINUE
10000  FORMAT (10A)
      END

C=======================================================================
      SUBROUTINE SCANSLASH (INLINE, ILINE, IEND, WORD)
C=======================================================================

C   --SCANSLASH scans the input line for the "\" word.  The "\" word
C   --may be a string of alphabetic characters or a single character.
C   --
C   --Parameters:
C   --   INLINE - IN - the input line
C   --   ILINE - IN/OUT - the start of the "\" word, returned end of
C   --      the "\" word
C   --   IEND - IN - the end of line
C   --   WORD - OUT - the "\" word (including "\")

      CHARACTER*(*) INLINE
      CHARACTER*(*) WORD

      CHARACTER CH

      ISTART = ILINE

      CH = INLINE(ILINE+1:ILINE+1)
      IF ((CH .GE. 'A') .AND. (CH .LE. 'Z'))
     &   CH = CHAR (ICHAR(CH) - ICHAR('A') + ICHAR('a'))
      ILINE = ILINE + 1

      IF ((CH .GE. 'a') .AND. (CH .LE. 'z')) THEN
         DO WHILE (ILINE .LT. IEND)
            CH = INLINE(ILINE+1:ILINE+1)
            IF ((CH .GE. 'A') .AND. (CH .LE. 'Z'))
     &         CH = CHAR (ICHAR(CH) - ICHAR('A') + ICHAR('a'))
            IF ((CH .LT. 'a') .OR. (CH .GT. 'z')) GOTO 100
            ILINE = ILINE + 1
         END DO
      END IF

  100 CONTINUE
      WORD = INLINE(ISTART:ILINE)

      RETURN
      END

C=======================================================================
      SUBROUTINE SCANPARAM (INLINE, ILINE, IEND, WORD)
C=======================================================================

C   --SCANPARAM scans the input line for the parameter.  The parameter
C   --starts with "{" and ends with "}" (stripped).
C   --
C   --Parameters:
C   --   INLINE - IN - the input line
C   --   ILINE - IN/OUT - the start of the parameter, returned end of
C   --      the parameter
C   --   IEND - IN - the end of line
C   --   WORD - OUT - the parameter without ("{" and "}")

      CHARACTER*(*) INLINE
      CHARACTER*(*) WORD

      CHARACTER CH

      ILINE = ILINE + 1
      CH = INLINE(ILINE:ILINE)
      DO WHILE (CH .EQ. ' ')
         ILINE = ILINE + 1
         CH = INLINE(ILINE:ILINE)
      END DO

      ISTART = ILINE
      IF (CH .NE. '{') GOTO 100

      L = LEN (INLINE)

      ILINE = ILINE + 1
      CH = INLINE(ILINE:ILINE)
      DO WHILE ((CH .NE. '}') .AND. (ILINE .LE. L))
         ILINE = ILINE + 1
         CH = INLINE(ILINE:ILINE)
      END DO

  100 CONTINUE
      WORD = INLINE(ISTART+1:ILINE-1)

      RETURN
      END
