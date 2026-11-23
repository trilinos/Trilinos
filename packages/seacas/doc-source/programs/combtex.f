C=======================================================================
      PROGRAM COMBTEX
C=======================================================================

C                       Written by Amy Gilkey
C                       Last Revised 04/25/88

C   --COMBTEX combines a TEX file by including (inline) all \include and
C   --\input files.

$HP9000_800 INTRINSICS ON
      CHARACTER*120 INLINE
      CHARACTER*120 FILNAM
      CHARACTER*20 WORD
      LOGICAL INPUT, ABSTRACT, NXTABS

C   --Open the input and output files

      arg_count = iargc()
      if (arg_count .ge. 2) then
         call getarg(2, filnam)
      else
         write (*,*) 'ERROR - No filename specified'
         goto 120
      end if
      NTEX = 20

      OPEN (NTEX, FILE=FILNAM, READONLY, STATUS='OLD', IOSTAT=IERR)
      IF (IERR .NE. 0) THEN
         PRINT *, 'ERROR - Cannot open input file'
         GOTO 120
      END IF

      NCOMB = 19

      OPEN (NCOMB, FILE='OUTFILE', STATUS='NEW')

C   --Process the file

  100 CONTINUE
      ABSTRACT = .FALSE.

      DO WHILE (.TRUE.)
         READ (NTEX, '(Q, A)', END=110) IEND, INLINE
         NCOMB = NCOMB
         INPUT = .FALSE.
         NXTABS = .FALSE.

C      --Scan for \include or \input

         IF (INLINE(1:1) .EQ. '\') THEN
            IE = INDEX (INLINE, '{') - 1
            WORD = INLINE(1:IE)

            IF ((WORD .EQ. '\include') .OR. (WORD .EQ. '\input')
     &         .OR. (WORD .EQ. '\abstract')) THEN
               IF ((WORD .EQ. '\include')
     &            .OR. (WORD .EQ. '\input')) THEN
                  INPUT = .TRUE.
               ELSE IF (WORD .EQ. '\abstract') THEN
                  NXTABS = .TRUE.
               END IF

C            --Change the input file (return to old file at end of file)

               IS = INDEX (INLINE, '{') + 1
               IE = INDEX (INLINE, '}') - 1
               IDIR = INDEX (INLINE, '/') + 1
               I = INDEX (INLINE(IDIR:), '.')
               IF (I .LE. 0) THEN
                  FILNAM = INLINE(IS:IE) // '.tex'
               ELSE
                  FILNAM = INLINE(IS:IE)
               END IF
               IE = INDEX (INLINE, '}')
               PRINT *, INLINE(:IE)

               NTEX = NTEX + 1
               OPEN (NTEX, FILE=FILNAM, READONLY,
     $              STATUS='OLD', IOSTAT=IERR)
               IF (IERR .NE. 0) THEN
                  OPEN (NTEX, FILE=FILNAM//'.tex', READONLY,
     $                 STATUS='OLD', IOSTAT=IERR)
                  IF (IERR .NE. 0) THEN
                     PRINT *, 'ERROR - Cannot open file ', FILNAM
                     NTEX = NTEX - 1
                  end if
               END IF
            END IF
         END IF

C      --Write processed line

         IF (INPUT .OR. ABSTRACT) THEN
            WRITE (NCOMB, 10000) '%%% ', INLINE(:IEND)
         ELSE
            WRITE (NCOMB, 10000) INLINE(:IEND)
         END IF
10000     FORMAT (10A)

         IF (NXTABS) ABSTRACT = .TRUE.
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
      CLOSE (NCOMB)

      END
