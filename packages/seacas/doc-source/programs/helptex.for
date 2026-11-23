C=======================================================================
      PROGRAM HELPTEX
C=======================================================================

C                       Written by Amy Gilkey
C                       Last Revised 06/01/88

C   --HELPTEX changes the command section of a TEX file into a RUNOFF file
C   --to generate a help file for the commands.

      CHARACTER*120 INLINE
      CHARACTER*120 FILNAM
      CHARACTER*20 MODE
      CHARACTER*20 WORD
      LOGICAL ALPHA
      LOGICAL BREAK

      PARAMETER (MAXSECT = 20, MAXSSECT = 50, MAXCMD = 500)

      CHARACTER*80 SECTION, SUBSECT
      CHARACTER*10 CMD

      CHARACTER*80 SECTTBL(MAXSECT)
      INTEGER IXSUB(0:MAXSECT)
      CHARACTER*80 SUBTBL(MAXSSECT)
      INTEGER IXCMD(0:MAXSSECT)
      CHARACTER*10 CMDTBL(MAXCMD)

C   --Open the input and output files

      NTEX = 20

      OPEN (NTEX, FILE='INFILE', READONLY, STATUS='OLD', IOSTAT=IERR)
      IF (IERR .NE. 0) THEN
         PRINT *, 'ERROR - Cannot open input file'
         GOTO 120
      END IF

      NCMD = 30

      OPEN (NCMD, FILE='CMDFILE', CARRIAGECONTROL='LIST',
     &   STATUS='NEW')

      WRITE (NCMD, 10040) '.NOAUTOPARAGRAPH'
      WRITE (NCMD, 10040) '.PARAGRAPH 0,0'
      WRITE (NCMD, 10040) '.NO AUTOJUSTIFY .NO JUSTIFY'
      WRITE (NCMD, 10040) '.LEFT MARGIN 0 .RIGHT MARGIN 72'
      WRITE (NCMD, 10040) '.NO PAGING'
      WRITE (NCMD, 10040) '.SET DATE .FLAGS SUBSTITUTE'

C   --Set up the section, subsection, and command tables

      IXSUB(0) = 0
      IXCMD(0) = 0

      NSECTS = 1
      SECTTBL(1) = 'Commands'
      IXSUB(1) = 1
      SUBTBL(1) = ' '
      IXCMD(1) = 0

C   --Process the file

C   --Set mode to scan for section or command description
      MODE = 'SCANNING'

  100 CONTINUE
      DO WHILE (.TRUE.)
         READ (NTEX, '(Q, A)', END=110) IEND, INLINE
         IWRITE = 0
         BREAK = .FALSE.

         IF (MODE .EQ. 'SCANNING') THEN

C         --Scan for section or command description

            IF (INLINE(1:1) .EQ. '\') THEN
               IL = 1
               CALL SCANSLASH (INLINE, IL, IEND, WORD, ALPHA)

               IF ((WORD .EQ. '\input')
     &            .OR. (WORD .EQ. '\include')) THEN

C               --Change the input file (return to old file at end of file)

                  IS = INDEX (INLINE, '{') + 1
                  IE = INDEX (INLINE, '}') - 1
                  IDIR = INDEX (INLINE, ']') + 1
                  I = INDEX (INLINE(IDIR:), '.')
                  IF (I .LE. 0) THEN
                     FILNAM = INLINE(IS:IE) // '.TEX'
                  ELSE
                     FILNAM = INLINE(IS:IE)
                  END IF
                  IE = INDEX (INLINE, '}')
                  PRINT *, INLINE(:IE)

                  NTEX = NTEX + 1
                  OPEN (NTEX, FILE=FILNAM, READONLY, STATUS='OLD',
     &               IOSTAT=IERR)
                  IF (IERR .NE. 0) THEN
                     PRINT *, 'ERROR - Cannot open file ', FILNAM
                     NTEX = NTEX - 1
                  END IF

               ELSE IF ((WORD .EQ. '\section')
     &            .OR. (WORD .EQ. '\subsection')) THEN

C               --Change mode to processing start of section

                  NPAR = 0
                  MODE = 'SECTION HEADER'

               ELSE IF (WORD .EQ. '\cmddef') THEN

C               --Change mode to processing start of command

                  NPAR = 0
                  MODE = 'COMMAND HEADER'
               END IF
            END IF
         END IF

         IF ((MODE .EQ. 'SECTION HEADER')
     &      .OR. (MODE .EQ. 'SECTION REST')) THEN

C         --Process section

            CALL STRIPLINE (INLINE, IEND, NPAR, BREAK)

            IF (MODE .EQ. 'SECTION HEADER') THEN
               IS = LIB$SKPC (' ', INLINE)
               IDUM = STR$TRIM (INLINE, INLINE, IE)
               IF (INLINE(IS:IS) .EQ. '[') THEN
                  IS = IS + 1
                  I = INDEX (INLINE, ']')
                  IF (I .GE. IS) IE = I - 1
               END IF

               IF (WORD .EQ. '\section') THEN

C               --Start a new section (add dummy subsection)

                  SECTION = INLINE(IS:IE)
                  NSECTS = NSECTS + 1
                  SECTTBL(NSECTS) = SECTION
                  IXSUB(NSECTS) = IXSUB(NSECTS-1) + 1
                  NSUBS = IXSUB(NSECTS)
                  SUBTBL(NSUBS) = ' '
                  IXCMD(NSUBS) = IXCMD(NSUBS-1)

               ELSE

C               --Start a new subsection

                  SUBSECT = INLINE(IS:IE)
                  IXSUB(NSECTS) = IXSUB(NSECTS) + 1
                  NSUBS = IXSUB(NSECTS)
                  SUBTBL(NSUBS) = SUBSECT
                  IXCMD(NSUBS) = IXCMD(NSUBS-1)
                  WRITE (NCMD, 10040) '.P; !'
                  WRITE (NCMD, 10040) '.P; ! ', SECTTBL(NSECTS)
                  WRITE (NCMD, 10040) '.P; ! ', SUBSECT
               END IF

C            --Change mode to scanning rest of section header
               MODE = 'SECTION REST'
            END IF

            IF (NPAR .LE. 0) THEN
C            --Change mode to done with section header
               MODE = 'SCANNING'
            END IF

         ELSE IF (MODE .EQ. 'COMMAND HEADER') THEN

C         --Process start of command description

C         --Write command description header
            CALL STRIPLINE (INLINE, IEND, NPAR, BREAK)
            IWRITE = NCMD

            IS = LIB$SKPC (' ', INLINE)
            IE = INDEX (INLINE(IS:), ' ') + IS - 2
            CMD = INLINE(IS:IE)
            NSUBS = IXSUB(NSECTS)
            IXCMD(NSUBS) = IXCMD(NSUBS) + 1
            NCMDS = IXCMD(NSUBS)
            CMDTBL(NCMDS) = CMD
            WRITE (NCMD, 10040) '.P; !'
            WRITE (NCMD, 10040) '.P; 3 ', CMD
            WRITE (NCMD, 10040) '.P; !'
            WRITE (NCMD, 10040) '.LM +5 .P -5,0'

C         --Change mode to process command parameters
            MODE = 'COMMAND PARAM'

         ELSE IF (MODE .EQ. 'COMMAND PARAM') THEN

C         --Process command parameters until start of description

            IF (INLINE .EQ. '} {') THEN
C            --Change mode to process command description
               WRITE (NCMD, 10040) '.S .LM -5 .P 0,0'
               MODE = 'COMMAND'
            ELSE
C            --Process command parameters
               CALL STRIPLINE (INLINE, IEND, NPAR, BREAK)
               IWRITE = NCMD
            END IF

         ELSE IF (MODE .EQ. 'COMMAND') THEN

C         --Process command description until end of description

            IF (IEND .LE. 0) THEN
C            --Blank line starts new paragraph
               WRITE (NCMD, 10040) '.S .P'
            ELSE IF ((INLINE .EQ. '}') .AND. (NPAR .LE. 1)) THEN
C            --Change mode to scanning for section or command description
               MODE = 'SCANNING'
            ELSE
C            --Write command description
               CALL STRIPLINE (INLINE, IEND, NPAR, BREAK)
               IWRITE = NCMD
            END IF
         END IF

C      --Write processed line

         IF (IWRITE .GT. 0) THEN
            WRITE (IWRITE, 10040) INLINE(:IEND)
            IF (BREAK) THEN
               IF (MODE .EQ. 'COMMAND PARAM') THEN
                  WRITE (IWRITE, 10040) '.P'
               ELSE
                  WRITE (IWRITE, 10040) '.BR'
               END IF
            END IF
         END IF
      END DO

  110 CONTINUE

C   --Close TEX file, if input file resume processing with previous file

      IF (NTEX .GT. 20) THEN
         CLOSE (NTEX)
         NTEX = NTEX - 1
         GOTO 100
      END IF

  120 CONTINUE

      CLOSE (NTEX)
      CLOSE (NCMD)

C   --Write the command summary

      NSUM = 31

      OPEN (NSUM, FILE='SUMFILE', CARRIAGECONTROL='LIST',
     &   STATUS='NEW')

      WRITE (NSUM, 10040) '!'
      WRITE (NSUM, 10040) '2 Commands'
      WRITE (NSUM, 10040) '!'
      WRITE (*, *)

      DO ISECT = 1, NSECTS
         NC = IXCMD(IXSUB(ISECT)) - IXCMD(IXSUB(ISECT-1))
         IF (NC .GT. 0) THEN
            IF (SUBTBL(IXSUB(ISECT)) .NE. ' ') THEN
               WRITE (NSUM, *)
               WRITE (*, *)
               IDUM = STR$TRIM (INLINE, SECTTBL(ISECT), IEND)
               WRITE (NSUM, 10000) SECTTBL(ISECT)(:IEND)
               WRITE (*, 10010) SECTTBL(ISECT)(:IEND)
10000           FORMAT (25X, A)
10010           FORMAT (26X, A)
            END IF

            DO ISUB = IXSUB(ISECT-1)+1, IXSUB(ISECT)
               NC = IXCMD(ISUB) - IXCMD(ISUB-1)
               IF (NC .GT. 0) THEN
                  WRITE (NSUM, *)
                  WRITE (*, *)
                  IF (SUBTBL(ISUB) .NE. ' ') THEN
                     IDUM = STR$TRIM (INLINE, SUBTBL(ISUB), IEND)
                     WRITE (NSUM, 10040) SUBTBL(ISUB)(:IEND)
                     WRITE (*, *) SUBTBL(ISUB)(:IEND)
                  ELSE
                     IDUM = STR$TRIM (INLINE, SECTTBL(ISECT), IEND)
                     WRITE (NSUM, 10040) SECTTBL(ISECT)(:IEND)
                     WRITE (*, *) SECTTBL(ISECT)(:IEND)
                  END IF
                  WRITE (NSUM, 10020) (CMDTBL(ICMD),
     &               ICMD = IXCMD(ISUB-1)+1, IXCMD(ISUB))
                  WRITE (*, 10030) (CMDTBL(ICMD),
     &               ICMD = IXCMD(ISUB-1)+1, IXCMD(ISUB))
10020              FORMAT ((3X, 6 (A, :, 1X)))
10030              FORMAT ((4X, 6 (A, :, 1X)))
               END IF
            END DO
         END IF
      END DO

10040  FORMAT (10A)
      END

C=======================================================================
      SUBROUTINE SCANSLASH (INLINE, ILINE, IEND, WORD, ALPHA)
C=======================================================================

C   --SCANSLASH scans the input line for the "\" word.  The "\" word
C   --may be a string of lowercase alphabetic characters or a single
C   --character.
C   --
C   --Parameters:
C   --   INLINE - IN - the input line
C   --   ILINE - IN/OUT - the start of the "\" word, returned end of
C   --      the "\" word
C   --   IEND - IN - the end of line
C   --   WORD - OUT - the "\" word (including "\")
C   --   ALPHA - OUT - true iff word is alphabetic

      CHARACTER*(*) INLINE
      CHARACTER*(*) WORD
      LOGICAL ALPHA

      CHARACTER CH

      ISTART = ILINE

      CH = INLINE(ILINE+1:ILINE+1)
      ALPHA = ((CH .GE. 'a') .AND. (CH .LE. 'z'))
      IF (CH .NE. ' ') ILINE = ILINE + 1

      IF (ALPHA) THEN
         DO WHILE (ILINE .LT. IEND)
            CH = INLINE(ILINE+1:ILINE+1)
            IF ((CH .LT. 'a') .OR. (CH .GT. 'z')) GOTO 100
            ILINE = ILINE + 1
         END DO
      END IF

  100 CONTINUE
      WORD = INLINE(ISTART:ILINE)

      RETURN
      END

C=======================================================================
      SUBROUTINE STRIPLINE (INLINE, IEND, NPAR, BREAK)
C=======================================================================

C   --STRIPLINE strips the input line of TEX control sequences.
C   --
C   --Parameters:
C   --   INLINE - IN/OUT - the input line, returned stripped
C   --   IEND - IN/OUT - the end of line
C   --   NPAR - IN/OUT - the parenthesis (i.e., { and }) level
C   --   BREAK - OUT - returned true iff a new line should start after
C   --      the current line

      CHARACTER*(*) INLINE

      CHARACTER CH
      CHARACTER*20 WORD
      CHARACTER*20 TYPE
      LOGICAL ALPHA
      LOGICAL BREAK
      CHARACTER*20 TYPE_STACK(0:20)

      SAVE TYPE, TYPE_STACK

      DATA TYPE_STACK(0) / 'KNOWN' /

      IF (NPAR .LE. 0) TYPE = TYPE_STACK(NPAR)

      ILINE = 1
      BREAK = .FALSE.

      DO WHILE (ILINE .LE. IEND)
         CH = INLINE(ILINE:ILINE)
         ISTART = ILINE

         IF (CH .EQ. '{') THEN

C         --Process '{' according to current control sequence type

            NPAR = NPAR + 1
            TYPE_STACK(NPAR) = TYPE

            IF (TYPE .EQ. 'PARAM') THEN
               INLINE(ILINE:ILINE) = '{'
            ELSE IF (TYPE .EQ. 'OPTPARAM') THEN
               INLINE(ILINE:ILINE) = '['
            ELSE IF (TYPE .EQ. 'DEFAULT') THEN
               INLINE(ILINE:ILINE) = '<'
            ELSE IF (TYPE .EQ. 'STRIP') THEN
               CONTINUE
            ELSE IF (TYPE .EQ. 'CMDNOTE') THEN
               CALL SUBLINE (INLINE, ISTART, ILINE, IEND, 'x')
            ELSE IF (TYPE .EQ. 'CMDOPTION') THEN
               CALL SUBLINE (INLINE, ISTART, ILINE, IEND, 'x')
            ELSE IF (TYPE .EQ. 'CMDNOTE2') THEN
               CALL SUBLINE (INLINE, ISTART, ILINE, IEND,
     &            '.LM +2 .P 0')
            ELSE IF (TYPE .EQ. 'CMDOPTION2') THEN
               CALL SUBLINE (INLINE, ISTART, ILINE, IEND,
     &            '.LM +2 .P 0')
            ELSE
               CALL SUBLINE (INLINE, ISTART, ILINE, IEND, 'x')
            END IF

         ELSE IF (CH .EQ. '}') THEN

C         --Process '}' according to current control sequence type

            TYPE = TYPE_STACK(NPAR)
            NPAR = NPAR - 1

            IF (TYPE .EQ. 'PARAM') THEN
               INLINE(ILINE:ILINE) = '}'
            ELSE IF (TYPE .EQ. 'OPTPARAM') THEN
               INLINE(ILINE:ILINE) = ']'
            ELSE IF (TYPE .EQ. 'DEFAULT') THEN
               INLINE(ILINE:ILINE) = '>'
            ELSE IF (TYPE .EQ. 'STRIP') THEN
               CALL SUBLINE (INLINE, ISTRIP, ILINE, IEND, 'x')
            ELSE IF (TYPE .EQ. 'CMDNOTE') THEN
               TYPE = 'CMDNOTE2'
               CALL SUBLINE (INLINE, ISTART, ILINE, IEND, ':')
            ELSE IF (TYPE .EQ. 'CMDOPTION') THEN
               TYPE = 'CMDOPTION2'
               CALL SUBLINE (INLINE, ISTART, ILINE, IEND, 'x')
            ELSE IF (TYPE .EQ. 'CMDNOTE2') THEN
               CALL SUBLINE (INLINE, ISTART, ILINE, IEND,
     &            '.LM -4 .P')
            ELSE IF (TYPE .EQ. 'CMDOPTION2') THEN
               CALL SUBLINE (INLINE, ISTART, ILINE, IEND,
     &            '.LM -4 .P')
            ELSE
               CALL SUBLINE (INLINE, ISTART, ILINE, IEND, 'x')
            END IF

         ELSE IF (CH .EQ. '$') THEN

C         --Change '$_' to '__', delete other '$'

            CH = INLINE(ILINE+1:ILINE+1)
            IF (CH .EQ. '_') THEN
               INLINE(ILINE:ILINE) = '_'
               ILINE = ILINE + 1
            ELSE
               CALL SUBLINE (INLINE, ISTART, ILINE, IEND, 'x')
            END IF

         ELSE IF (CH .EQ. '`') THEN

C         --Change '``' to '"'

            CH = INLINE(ILINE+1:ILINE+1)
            IF (CH .EQ. '`') THEN
               ILINE = ILINE + 1
               CALL SUBLINE (INLINE, ISTART, ILINE, IEND, '"')
            END IF

         ELSE IF (CH .EQ. '''') THEN

C         --Change '''' to '"'

            CH = INLINE(ILINE+1:ILINE+1)
            IF (CH .EQ. '''') THEN
               ILINE = ILINE + 1
               CALL SUBLINE (INLINE, ISTART, ILINE, IEND, '"')
            END IF

         ELSE IF (CH .EQ. '-') THEN

C         --Change '--' to '-'

            CH = INLINE(ILINE+1:ILINE+1)
            IF (CH .EQ. '-') THEN
               CALL SUBLINE (INLINE, ISTART, ILINE, IEND, 'x')
            END IF

         ELSE IF (CH .EQ. '~') THEN

C         --Delete '~'

            CALL SUBLINE (INLINE, ISTART, ILINE, IEND, ' ')

         ELSE IF ((CH .EQ. '\')
     &      .AND. (TYPE_STACK(NPAR) .NE. 'STRIP')) THEN

C         --Process control sequence

            CALL SCANSLASH (INLINE, ILINE, IEND, WORD, ALPHA)

            IF (.NOT. ALPHA) THEN
               IF (WORD .EQ. '\\') THEN
C               --Delete new line indicator and set BREAK
                  BREAK = .TRUE.
                  CALL SUBLINE (INLINE, ISTART, ILINE, IEND, 'x')
               ELSE IF (WORD .EQ. '\') THEN
C               --Delete space indicator 
                  CALL SUBLINE (INLINE, ISTART, ILINE, IEND, 'x')
               ELSE
C               --Change '\' to '_'
                  INLINE(ISTART:ISTART) = '_'
               END IF

            ELSE
               TYPE = ' '
               IF (WORD .EQ. '\param') THEN
                  CALL SUBLINE (INLINE, ISTART, ILINE, IEND, 'x')
                  TYPE = 'PARAM'

               ELSE IF (WORD .EQ. '\optparam') THEN
                  CALL SUBLINE (INLINE, ISTART, ILINE, IEND, 'x')
                  TYPE = 'OPTPARAM'

               ELSE IF (WORD .EQ. '\label') THEN
                  ISTRIP = ISTART
                  TYPE = 'STRIP'

               ELSE IF (WORD .EQ. '\setlength') THEN
                  ISTRIP = ISTART
                  TYPE = 'STRIP'

               ELSE IF (WORD .EQ. '\notetome') THEN
                  ISTRIP = ISTART
                  TYPE = 'STRIP'

               ELSE IF (WORD .EQ. '\ref') THEN
                  PRINT *, '*** ERROR: invalid command ', WORD
                  PRINT *, '    ', INLINE(1:IEND)
                  CALL SUBLINE (INLINE, ISTART, ILINE, IEND,
     &               'in manual')
                  ISTRIP = ILINE + 1
                  TYPE = 'STRIP'

               ELSE IF (WORD .EQ. '\begin') THEN
                  IF (INLINE(ILINE+1:ILINE+9) .EQ. '{itemize}') THEN
                     ILINE = ILINE + 9
                     CALL SUBLINE (INLINE, ISTART, ILINE, IEND,
     &                  '.LIST 0,"o"')
                  END IF
               ELSE IF (WORD .EQ. '\end') THEN
                  IF (INLINE(ILINE+1:ILINE+9) .EQ. '{itemize}') THEN
                     ILINE = ILINE + 9
                     CALL SUBLINE (INLINE, ISTART, ILINE, IEND,
     &                  '.END LIST 0')
                  END IF
               ELSE IF (WORD .EQ. '\item') THEN
                  IF (ILINE .LT. IEND) THEN
                     CALL SUBLINE (INLINE, ISTART, ILINE, IEND,
     &                  '.LE;')
                  ELSE
                     CALL SUBLINE (INLINE, ISTART, ILINE, IEND,
     &                  '.S .LE')
                  END IF

               ELSE IF (WORD .EQ. '\cenlinesbegin') THEN
                  CALL SUBLINE (INLINE, ISTART, ILINE, IEND,
     &               '.LM +8 .S .P')
               ELSE IF (WORD .EQ. '\cenlinesend') THEN
                  CALL SUBLINE (INLINE, ISTART, ILINE, IEND,
     &               '.LM -8 .S .P')

               ELSE IF (WORD .EQ. '\nodefault') THEN
                  CALL SUBLINE (INLINE, ISTART, ILINE, IEND,
     &               '<no default>')
               ELSE IF (WORD .EQ. '\default') THEN
                  CALL SUBLINE (INLINE, ISTART, ILINE, IEND, 'x')
                  TYPE = 'DEFAULT'

               ELSE IF (WORD .EQ. '\cmdnote') THEN
                  CALL SUBLINE (INLINE, ISTART, ILINE, IEND,
     &               '.LM +2 .P -2; Issued from ')
                  TYPE = 'CMDNOTE'
               ELSE IF (WORD .EQ. '\cmdoption') THEN
                  CALL SUBLINE (INLINE, ISTART, ILINE, IEND, 
     &               '.LM +2 .P -2; ')
                  TYPE = 'CMDOPTION'

               ELSE IF ((WORD .EQ. '\cmd')
     &            .OR. (WORD .EQ. '\cmdverb')
     &            .OR. (WORD .EQ. '\caps')
     &            .OR. (WORD .EQ. '\bold')) THEN
                  CALL SUBLINE (INLINE, ISTART, ILINE, IEND, 'x')
               ELSE IF ((WORD .EQ. '\section')
     &            .OR. (WORD .EQ. '\subsection')) THEN
                  CALL SUBLINE (INLINE, ISTART, ILINE, IEND, 'x')
               ELSE IF ((WORD .EQ. '\cmddef')
     &            .OR. (WORD .EQ. '\cmdnext')) THEN
                  CALL SUBLINE (INLINE, ISTART, ILINE, IEND, 'x')
               ELSE IF (WORD .EQ. '\nth') THEN
                  CALL SUBLINE (INLINE, ISTART, ILINE, IEND, 'x')

               ELSE IF (WORD .EQ. '\ldots') THEN
                  CALL SUBLINE (INLINE, ISTART, ILINE, IEND, '_._._.')

               ELSE
                  PRINT *, '*** ERROR: unknown command ', WORD
                  PRINT *, '    ', INLINE(1:IEND)
               END IF
            END IF
         ELSE
            CONTINUE
         END IF
         ILINE = ILINE + 1
      END DO

      RETURN
      END

C=======================================================================
      SUBROUTINE SUBLINE (INLINE, ISTART, ILINE, IEND, SUBSTR)
C=======================================================================

C   --SUBLINE substitutes a string in the input line.
C   --
C   --Parameters:
C   --   INLINE - IN/OUT - the input line, returned with substituted string
C   --   ISTART - IN/OUT - the start of the string to be substituted;
C   --      set to one if INLINE(:ISTART-1) is empty
C   --   ILINE - IN/OUT - the end of the string to be substituted;
C   --      returned end of the substituted string
C   --   IEND - IN/OUT - the end of line
C   --   SUBSTR - OUT - the string to be substituted (including any leading
C   --      and trailing blanks; null string if 'x'

      CHARACTER*(*) INLINE
      CHARACTER*(*) SUBSTR

      IF (INLINE(:ISTART-1) .EQ. ' ') ISTART = 1
      LOUT = ILINE - ISTART + 1
      LSUB = LEN (SUBSTR)
      IF (SUBSTR .EQ. 'x') LSUB = 0
      INLINE(ISTART:) = SUBSTR(:LSUB) // INLINE(ILINE+1:IEND)
      IEND = IEND - LOUT + LSUB
      ILINE = ISTART + LSUB - 1
      IF ((SUBSTR(1:1) .EQ. '.') .AND. (ISTART .GT. 1)) THEN
         PRINT *, '*** ERROR: substituting ', SUBSTR
         PRINT *, '    ', INLINE(1:IEND)
      END IF

      RETURN
      END
