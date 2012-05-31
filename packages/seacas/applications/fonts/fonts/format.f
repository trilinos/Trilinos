      DIMENSION IDEX(200), NVECT(200), XSIZE(200), YSIZE(200), X0(2300),
     * Y0(2300), X1(2300), Y1(2300)
      DIMENSION RVAL(4)
      CHARACTER*32FILE1, FILE2
      INTEGER OUTFONT
      LOGICAL FIRST, LAST
      INFONT = 12
      OUTFONT = 13
10    FORMAT(1X, A)
2000  CONTINUE
        WRITE(6, 10) '_Enter input font file name:'
        READ(5, '(a)') FILE1
        OPEN(UNIT = INFONT, FILE = FILE1, STATUS = 'old', IOSTAT = IOS)
        IF (IOS .NE. 0) THEN
          WRITE(6, *) ' Couldn''t open input font file; re-enter'
          STOP
        ENDIF
        GOTO 2020
2010  GOTO 2000
2020  CONTINUE
2030  CONTINUE
        WRITE(6, 10) '_Enter output font file name:'
        READ(5, '(a)') FILE2
        OPEN(UNIT = OUTFONT, FILE = FILE2, STATUS = 'new', FORM = 'unfor
     *matted', IOSTAT = IOS)
        IF (IOS .NE. 0) THEN
          WRITE(6, *) ' Couldn''t open output font file; re-enter'
          STOP
        ENDIF
        GOTO 2050
2040  GOTO 2030
2050  CONTINUE
      READ(INFONT, *, ERR=20) CHSP
      READ(INFONT, *, ERR=20) POINT
      DO 2060 IC = 1, 2
        LAST = .FALSE.
        NTOTVC = 0
        FIRST = .TRUE.
        NLOCVC = 0
        XMAX = 0.
        YMAX = 0.
        MAXVECT = 0
2080    CONTINUE
          READ(INFONT, *, END=100, ERR=20) (RVAL(I), I = 1, 4)
          IF (RVAL(1) .EQ. -99.) THEN
            IF (FIRST) THEN
              FIRST = .FALSE.
              ICH = RVAL(4)
              IDEX(ICH) = 1
              GOTO 2090
            ELSE
              IF (RVAL(4) .EQ. 999.) THEN
                LAST = .TRUE.
                GOTO 2100
              ENDIF
              IF (ICH .GT. 32) THEN
                XSIZE(ICH) = (XMAX + CHSP)/POINT
              ELSE
                XSIZE(ICH) = (XMAX)/POINT
              ENDIF
              YSIZE(ICH) = YMAX/POINT
              NVECT(ICH) = NLOCVC
              IF (MAXVECT .LT. NLOCVC) THEN
                MAXVECT = NLOCVC
                MAXCHR = ICH
              ENDIF
              ICH = RVAL(4)
              IDEX(ICH) = NTOTVC + 1
              XMAX = 0.
              YMAX = 0.
              NLOCVC = 0
            ENDIF
          ELSE
            NTOTVC = NTOTVC + 1
            NLOCVC = NLOCVC + 1
            X0(NTOTVC) = RVAL(1)/POINT
            Y0(NTOTVC) = RVAL(2)/POINT
            X1(NTOTVC) = RVAL(3)/POINT
            Y1(NTOTVC) = RVAL(4)/POINT
            IF (RVAL(1) .GT. XMAX) THEN
              XMAX = RVAL(1)
            ENDIF
            IF (RVAL(3) .GT. XMAX) THEN
              XMAX = RVAL(3)
            ENDIF
            IF (RVAL(2) .GT. YMAX) THEN
              YMAX = RVAL(2)
            ENDIF
            IF (RVAL(4) .GT. YMAX) THEN
              YMAX = RVAL(4)
            ENDIF
          ENDIF
          IF (LAST) THEN
            GOTO 2100
          ENDIF
2090    GOTO 2080
2100    CONTINUE
100     CONTINUE
        IF (ICH .GT. 32) THEN
          XSIZE(ICH) = (XMAX + CHSP)/POINT
        ELSE
          XSIZE(ICH) = (XMAX)/POINT
        ENDIF
        YSIZE(ICH) = YMAX/POINT
        NVECT(ICH) = NLOCVC
        DO 2110 I = 48, 57
          XSIZE(I) = 1.
2110    CONTINUE
        XSIZE(43) = 1.
        XSIZE(45) = 1.
        WRITE(6, *) ' Total number of vectors in font is ', NTOTVC
        WRITE(6, *) ' The most number of vectors per character is ', MAX
     *VECT
        WRITE(6, *) '        for character ', MAXCHR
        NVECT(32) = 0
        WRITE(6, *) ' Number of vectors for character 32 (blank) set to 
     *zero.'
        WRITE(OUTFONT) NTOTVC, CHSP, POINT
        WRITE(OUTFONT) (IDEX(I), I=1,200)
        WRITE(OUTFONT) (NVECT(I),I=1,200)
        WRITE(OUTFONT) (XSIZE(I),I=1,200)
        WRITE(OUTFONT) (YSIZE(I),I=1,200)
        WRITE(OUTFONT) (X0(I),I=1,NTOTVC)
        WRITE(OUTFONT) (Y0(I),I=1,NTOTVC)
        WRITE(OUTFONT) (X1(I),I=1,NTOTVC)
        WRITE(OUTFONT) (Y1(I),I=1,NTOTVC)
2060  CONTINUE
      CLOSE(INFONT)
      CLOSE(OUTFONT)
      STOP
20    WRITE(6, *) ' Bad format for input'
      STOP
21    WRITE(6, *) 'Error in opening font file.'
      END
