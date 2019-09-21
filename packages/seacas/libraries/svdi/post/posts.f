C Copyright (C) 2009-2017 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.
C
C     * Neither the name of NTESS nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C

      SUBROUTINE POST
c
c read frame skip parameter from unit 21 (if it is there)
c     READ(21,*,ERR=10,END=10) IFRAM
c     CALL FRMSKP(IFRAM)
C
c call postprocessor loop routine
10    CALL PPLOOP
      RETURN
      END
C
C
C
      SUBROUTINE PPLOOP
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C
C
C LOOP             - LOOP THROUGH BGP FILE, CALLING VDI EQUIVALENTS
C
C P. WATTERBERG    - 25 OCT 81
C P. L. CROTTY     -  4 APR 85 - UPDATETO INCLUDE POLYGON FILL (PLC)
C K. Cole          - 05 oct 90 - added SAVE
C
C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C
C ENTRY CONDITIONS - VDI has been initialized by calling VDINIT.
C
C CALLS            - Most of the VDI routines.
C
C EXIT CONDITIONS  -
C
C NARRATIVE        - SCAN A BGP FILE, PICKING OFF OPCODES AND OPERANDS.
C                    MAKE A CALL TO THE EQUIVALENT VDI ROUTINE.
C                    THE OPCODE IS USED AS AN INDEX INTO A JUMP TABLE OF
C                    A COMPUTED GOTO.
C
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C
C
C          VARIABLES THAT DON'T HAVE ANYTHING TO DO WITH COLOR
C
      SAVE
      INTEGER OPCODE, COUNT, TEXT(132), ESNUMB(100)
      INTEGER NUMESC, IFRAM
      REAL ARGS(2000)
      LOGICAL EOFOK, NOEOF, MARKER, TWODIM, FICHE, START, DOFISH
C (PLC)
      REAL XARRAY(508),YARRAY(508)
C
C          DESCRIPTION OF THE ABOVE VARIABLES
C
C OPCODE      - When a BGP opcode is expected, the next 8 bits from the BGP
C               file are read into OPCODE.  OPCODE-127 is then used as an index
C               into a jump table.  Also, may be the first eight bits of an
C               X coordinate.
C COUNT       - The next 8 bits after the opcode which is the word count for
C               the opcode.  May be the second eight bits of an X coordinate.
C TEXT(132)   - When a text string command is seen, the text is read
C               into TEXT one character or 8 bits at a time.  This means that
C               TEXT will contain ascii characters, one per word, right
C               justified and zero filled, just what VDTEXT needs.
C ESNUMB(100) - Is a list of known VDI escape codes with numeric arguments.
C NUMESC      - The number of elements in
C IFRAM       - Draws every Ith frame.
C ARGS(2000)  - A real array to hold the arguments of the escape code.
C EOFOK       - Always .TRUE.  When requesting the next few bits from the file,
C               it is passed to the retrieval routine, PPBTR, to let it know
C               that an end of file is ok at this point.
C NOEOF       - Always .FALSE. When requesting the next few bits from the file,
C               it is passed to the retrieval routine, PPBTR, to let it know
C               that an end of file here is illegal.
C MARKER      - True when marker plot mode has been selected, false otherwise.
C TWODIM      - True when two dimensional coordinate mode is in effect and
C               false when in three dimensional mode.
C FICHE       - True when the device that is being post processed to will
C               handle the VDI fiche break escape.
C START       - True for the first file id command then false.
C DOFISH      - True if the VDI fiche break escape should be called rather than
C               a VDNWPG.
C (PLC)
C XARRAY,YARRAY-Real arrays to pass x,y coordinates to polygon routine (VDPOLY)
C
C          VARIABLES TO HANDLE COLOR
C
      INTEGER NEXTFC, NEXTBC
      INTEGER CSPOT, INDICES(256), NXTDEX
      INTEGER NUMCLR, COLMAP(0:255), INVMAP(0:255)
      REAL RGB(3,256), TABLE(3,0:255), VECTOR(7)
      LOGICAL BATCH, DRAWN, FIRST
C
C          DESCRIPTION OF THE ABOVE VARIABLES
C
C TABLE(3,0:255) - This is used for holding the color table as the meta file
C                  thinks it is.  The device may not have all 256 colors that
C                  the meta file has, so we can't just ship the color calls
C                  down to the device.
C COLMAP(0:255)  - This is the mapping from TABLE to the devices real color
C                  table.  COLMAP(n) is -1 if color n in TABLE has not been
C                  sent to the device yet.  If color n in TABLE has been
C                  defined at the device, then COLMAP(n) is the color number
C                  in the devices real table.  The basic idea is that when the
C                  device only has 8 colors, the user will get the first eight
C                  colors that are used in each frame.
C INVMAP(0:255)  - This is the inverse map from the devices real color table
C                  to TABLE.  Like COLMAP, INVMAP(n) is -1 if the link has not
C                  been established.
C NXTDEX         - The next position in the devices color table that has not
C                  been used.  That is, the next slot in INVMAP that is -1.
C NUMCLR         - The number of colors that the device can support.
C NEXTFC         - The color that the next primitive should be drawn in.
C NEXTBC         - The background color that should be used for the next VDNWPG
C BATCH          - True when in the middle of a batch update
C                - of the color table.
C RGB(3,256)     - When doing a batch update of the color table, this array is
C                  used to hold the color entries which can be sent to the
C                  device, i.e. those entries for which COLMAP is not -1.
C INDICES(256)    - The color table indexes associated with the entries in RGB.
C CSPOT          - The next open spot in RGB and INDICES.
C VECTOR(7)      - Used to find out the default foregound and background
C                  colors from VDIQOS and to reinitialize the system at a
C                  new file id.
C DRAWN          - True if something has been drawn since inititalization or
C                  the background color has changed from default.
C FIRST          - True until after the first newpage since initialization
C
C
C          NON COLOR DATA STATEMENTS
C
      DATA NUMESC   /     58  /
      DATA (ESNUMB(I),I=1,58)
     +  / 99, 100, 101, 200, 201, 202, 203, 204, 205, 206,
     +   206, 208, 209, 210, 211, 212, 213, 214, 215, 216,
     +   302, 400, 401, 402, 501, 502, 600, 601, 700, 701,
     +   702, 703, 704, 705, 706, 707, 708, 900, 1020, 1021,
     +   1022, 1023, 1024, 1025, 1026, 1400, 1505, 1506,1507,
     +   1508, 1509, 1510, 1700, 2100, 2101, 2300, 2305, 2306 /
      DATA DOFISH   / .FALSE. /
      DATA EOFOK    /  .TRUE. /
      DATA NOEOF    / .FALSE. /
      DATA FICHE    / .FALSE. /
      DATA START    /  .TRUE. /
      DATA IFRAM    /      1  /
      DATA IFCOUN   /      0  /
      DATA COUNT    /      0  /
      GOTO 80000
C
      ENTRY FRMSKP(IFRMS)
      IFRAM = IFRMS
      RETURN
C
80000 CONTINUE
C
C          MAIN LOOP TO PROCESS FILE
C
C          THROW AWAY UNUSED WORDS OF THE PREVIOUS INSTRUCTION
C            (NOTE: first time through count is 0)
C
    1 DO 2 I=1,COUNT
         CALL PPBTR(16,OPCODE,NOEOF)
    2 CONTINUE
C
C          GET NEXT INSTRUCTION OR MOVE/DRAW
C
    3 CALL PPBTR(8,OPCODE,EOFOK)
      CALL PPBTR(8,COUNT,EOFOK)
c
c          If just starting up, check for file id (84HEX=132)
C            or escape (82HEX=130)
C
      IF (START) THEN
         IF (OPCODE.NE.132.AND.OPCODE.NE.130) THEN
            PRINT 10010
10010       FORMAT(1X,
     +        'NO FILE ID PRESENT, POST PROCESSING TERMINATING.')
            RETURN
         ENDIF
      ENDIF
C
C          CHECK FOR MOVE OR DRAW
C
      IF(OPCODE.GE.128) GOTO 20
          XCOORD = FLOAT(256*OPCODE+COUNT)*SCALE
          CALL PPBTR(8,OPCODE,NOEOF)
          CALL PPBTR(8,COUNT,NOEOF)
          IF(.NOT.TWODIM) CALL PPBTR(16,OPCODE,NOEOF)
C
C          IF WE ARE SKIPPING THIS FRAME THEN DON'T DO ANYTHING
C
          IF(IFCOUN.NE.0) GOTO 3
C
C          IF OPCODE < 128 THEN MOVE, ELSE DRAW A POINT OR LINE
C
          IF(OPCODE.LT.128) THEN
              YCOORD = FLOAT(256*OPCODE+COUNT)*SCALE
              CALL VDMOVA(XCOORD,YCOORD)
            ELSE
C
C          MAKE SURE THE COLOR WE WANT HAS BEEN SET UP
C
              IF(COLMAP(NEXTFC).EQ.-1) THEN
                  INVMAP(NXTDEX) = NEXTFC
                  COLMAP(NEXTFC) = NXTDEX
                  CALL VDSTFC(NXTDEX)
                  CALL VDSTCO(1,NXTDEX,TABLE(1,NEXTFC),0)
C
C          SEARCH FOR THE NEXT UNUSED ENTRY IN INVMAP
C
    4             NXTDEX = NXTDEX + 1
                  IF(INVMAP(NXTDEX).NE.-1) GOTO 4
                  ENDIF
C
C          OK, MAKE A SMUDGE OF SOME SORT
C
              DRAWN = .TRUE.
              YCOORD = FLOAT(256*(OPCODE-128)+COUNT)*SCALE
              IF(MARKER) THEN
                  CALL VDPNTA(XCOORD,YCOORD)
                ELSE
                  CALL VDLINA(XCOORD,YCOORD)
                ENDIF
            ENDIF
          GOTO 3
C
C          WE HAVE AN OPCODE, SO BRANCH TO APPROPRIATE CODE
C
   20 IF(OPCODE.EQ.133) THEN
         IFCOUN = IFCOUN + 1
         IF(IFCOUN.EQ.IFRAM) IFCOUN = 0
         ENDIF
      IF(IFCOUN.NE.0.AND.OPCODE.NE.134) GOTO 1
      GOTO (
     +     1,    1,  300,  400,  500,  600,  700,    8,    8, 1000,
C         80    81    82    83    84    85    86    87    88    89
C
     +     8,    8,    8,    8,    8,    8, 1700, 1800, 1900, 2000,
C         8A    8B    8C    8D    8E    8F    90    91    92    93
C
     +     8,    8,    8,    8,    8,    8,    8,    8,    8,    8,
C         94    95    96    97    98    99    9A    9B    9C    9D
C
     +     8,    8, 3300, 3400, 3500, 3600, 3700,    1,    3,    3,
C         9E    9F    A0    A1    A2    A3    A4    A5    A6    A7
C
     +  4100, 4200, 4400,    8,    8,    8,    8,    8, 4900,    1,
C         A8    A9    AA    AB    AC    AD    AE    AF    B0    B1
C
     +  5100,    1,    1,    1,    8,    8,    8,    8,    8,    8,
C         B2    B3    B4    B5    B6    B7    B8    B9    BA    BB
C
     +     8,    8,    8,    8, 6500, 6600, 6700, 6800, 6900, 7000,
C         BC    BD    BE    BF    C0    C1    C2    C3    C4    C5
C
     +     1, 7200, 7300,    8,    8,    8,    8,    8,    8,    8,
C         C6    C7    C8    C9    CA    CB    CC    CD    CE    CF
C
     +     8,    8,    8,    8,    8,    8,    8,    8,    8,    8,
C         D0    D1    D2    D3    D4    D5    D6    D7    D8    D9
C
     +     8,    8,    8,    8,    8,    8,    8,    8,    8,    8,
C         DA    DB    DC    DD    DE    DF    E0    E1    E2    E3
C
     +     8,    8,    8,    8,    8,    8,    8,    8,    8,    8,
C         E4    E5    E6    E7    E8    E9    EA    EB    EC    ED
C
     +     8,    8,    8,    8,    8,    8,    8,    8,    8,    8,
C         EE    EF    F0    F1    F2    F3    F4    F5    F6    F7
C
     +     8,    8,    8,    8,    8,    8,    8,    8), OPCODE-127
C         F8    F9    FA    FB    FC    FD    FE    FF
C
C
C          82 -- ESCAPE CODES
C
  300 CALL PPBTR(8,OPCODE,NOEOF)
      CALL PPBTR(8,ICOUNT, NOEOF)
      COUNT = COUNT - 1
      IF(OPCODE.NE.1) GOTO 1
C
C          IT'S A VDI ESCAPE SO GET THE ESCAPE CODE
C
      CALL PPBTR(16,OPCODE,NOEOF)
      COUNT = COUNT - 1
C
C          CHECK TO SEE IF IT IS A KNOWN NUMERIC ESCAPE
C
      DO 310 I=1, NUMESC
         IF(OPCODE.EQ.ESNUMB(I)) GOTO 320
  310 CONTINUE
      GOTO 340
C
C          GET THE NUMERIC ARGUMENTS
C
  320 ISTOP = MIN0(COUNT,2000)
      DO 330 I=1, ISTOP/2
         CALL PPBTR(16,INTEGR,NOEOF)
         CALL PPBTR(16,IFRACT,NOEOF)
         IF(INTEGR.GE.32768) THEN
               ARGS(I) = -(FLOAT(INTEGR-32768) + IFRACT/32768.)
            ELSE
               ARGS(I) = FLOAT(INTEGR) + IFRACT/32768.
            ENDIF
  330 CONTINUE
      GOTO 350
C
C          IS IT A KNOWN ALHPA ESCAPE?
C
  340 CONTINUE
      GOTO 1
C
  350 CALL VDESCP(OPCODE,ISTOP/2,ARGS)
      COUNT = COUNT - ISTOP
      GOTO 1
C
C          83 -- ASPECT RATIO DEFINITION
C
  400 CALL PPBTR(16,IXDIM,NOEOF)
      CALL PPBTR(16,IYDIM,NOEOF)
      CALL PPBTR(16, JUNK,NOEOF)
      SCALE = AMIN1(XMAX/IXDIM,YMAX/IYDIM)
      GOTO 3
C
C          84 -- NEW FILE ID SO REINITIALIZE EVERYTHING
C          IF FIRST FILE ID, THEN DO ALL THE ONE TIME INITIALIZATION STUFF
C
  500 IF (START) THEN
c         CALL VBPKG('POST    ')
         CALL VDINIT(0.,5)
         CALL VDFRAM(0)
C
C          FIND OUT MAXIMUM X AND Y
C          USED TO SCALE BGP COORDINATES TO NDC COORDINATES
C
         CALL VDIQND(XMAX,YMAX)
C
C          FIND OUT DEFAULT FOREGROUND AND BACKGROUND COLORS
C
         CALL VDIQOS(VECTOR)
         NEXTFC = VECTOR(1)
         NEXTBC = VECTOR(2)
C
C          FIND OUT THE NUMBER OF COLORS IN THE COLOR TABLE
C
         INDEXA = 4
         CALL VDIQDC(4,DUMMY)
         NUMCLR = DUMMY
C
C          FIND OUT IF WE CAN USE FICHE BREAKS
C
         CALL VDIQES(211,IANS)
         IF(IANS.EQ.2) FICHE = .TRUE.
      ENDIF
      CSPOT = 1
      DRAWN = .FALSE.
      FIRST = .TRUE.
      MARKER = .FALSE.
      TWODIM = .TRUE.
      BATCH  = .FALSE.
      SCALE = AMIN1(XMAX/32767.,YMAX/32767.)
C
C          RE-ESTABLISH DEFAULT COLOR TABLE
C
      DO 510 I=0,255
         COLMAP(I) = I
         TABLE(1,I) = 1.
         TABLE(2,I) = 1.
         TABLE(3,I) = 1.
  510 CONTINUE
      TABLE(1,0) = 0.
      TABLE(2,0) = 0.
      TABLE(3,0) = 0.
      TABLE(2,1) = 0.
      TABLE(3,1) = 0.
      TABLE(1,2) = 0.
      TABLE(3,2) = 0.
      TABLE(3,3) = 0.
      TABLE(1,4) = 0.
      TABLE(2,4) = 0.
      TABLE(2,5) = 0.
      TABLE(1,6) = 0.
      IF(.NOT.START) THEN
         CALL VDWAIT
C
C          DELETE ALL SEGMENTS
C
         CALL VDESCP(203,1,0.)
C
C          REESTABLISH DEFAULT SETTINGS
C
         CALL VDSTOS(VECTOR)
         CALL VDSTCO(255,COLMAP,TABLE,0)
         CALL VDNWPG
         CALL VDMOVA(0.,0.)
         ENDIF
      START = .FALSE.
      NXTDEX = 256
C
C          IF THE DEVICE HAS LESS THAN 250 COLORS, THEN
C          SET UP THE MAPPINGS SO THAT THE USER WILL GET THE FIRST
C          FEW COLORS THEY ASK FOR.
C
      IF(NUMCLR.LE.250) THEN
         DO 520 I=0,255
            INVMAP(I) = -1
            COLMAP(I) = -1
  520    CONTINUE
         NXTDEX = 0
         IF(NEXTBC.EQ.0) NXTDEX = 1
         COLMAP(NEXTBC) = NEXTBC
         ENDIF
      NEXTFC = VECTOR(1)
      NEXTBC = VECTOR(2)
      GOTO 1
C
C          85 -- NEW PAGE
C
C          IF NOTHING HAS BEEN DRAWN, THEN DON'T CALL FOR A NEW PAGE.
C
  600 IF(FIRST.AND..NOT.DRAWN) GOTO 1
      FIRST = .FALSE.
      CALL VDWAIT
C
C          RESET THE COLOR MAPS SO THE USER CAN GET A NEW SET FOR NEXT FRAM
C
      IF(NUMCLR.LT.250) THEN
         DO 610 I= 0, 255
            COLMAP(I) = -1
            INVMAP(I) = -1
  610    CONTINUE
         COLMAP(NEXTBC) = 0
         CALL VDSTCO(1,0,TABLE(1,NEXTBC),0)
         NXTDEX = 1
         ENDIF
      CALL VDSTBC(COLMAP(NEXTBC))
      IF(DOFISH) THEN
            CALL VDESCP(211,0,0.)
         ELSE
            CALL VDNWPG
         ENDIF
      DOFISH = .FALSE.
      GOTO 1
C
C          86 -- END OF DATA FILE
C
  700 CALL VDWAIT
      CALL VDFRAM(1)
      CALL VDTERM
      RETURN
C
C          89 -- FICHE BREAK
C
 1000 DOFISH = FICHE
      GOTO 600
C
C          90 -- DEFINE COLOR TABLE INDEX
C
 1700 CALL PPBTR(16,INDEX,NOEOF)
      CALL PPBTR(16,IR,NOEOF)
      CALL PPBTR(8,IG,NOEOF)
      CALL PPBTR(8,IY,NOEOF)
      CALL PPBTR(8,IB,NOEOF)
      CALL PPBTR(8,IM,NOEOF)
      CALL PPBTR(8,IC,NOEOF)
      CALL PPBTR(8,IW,NOEOF)
      TABLE(1,INDEX) = AMIN1(IR/256.+IY/256.+IM/256.+IW/256.,1.)
      TABLE(2,INDEX) = AMIN1(IG/256.+IY/256.+IC/256.+IW/256.,1.)
      TABLE(3,INDEX) = AMIN1(IB/256.+IC/256.+IM/256.+IW/256.,1.)
C
C          IF THE LINK HAS ALREADY BEEN ESTABLISHED, THEN WE MAY WANT TO
C          CHANGE THE COLOR DYNAMICALLY.
C
      IF(COLMAP(INDEX).NE.-1) THEN
         IF(.NOT.BATCH) THEN
               CALL VDSTCO(1,COLMAP(INDEX),TABLE(1,INDEX),0)
            ELSE
               INDICES(CSPOT) = COLMAP(INDEX)
               RGB(1,CSPOT) = TABLE(1,INDEX)
               RGB(2,CSPOT) = TABLE(2,INDEX)
               RGB(3,CSPOT) = TABLE(3,INDEX)
               CSPOT = MIN0(CSPOT+1,256)
            ENDIF
         ENDIF
      GOTO 3
C
C          91 -- SELECT COLOR
C
 1800 CALL PPBTR(16,INDEX,NOEOF)
      IF(COLMAP(INDEX).GE.0) CALL VDSTFC(COLMAP(INDEX))
      NEXTFC = INDEX
      GOTO 3
C
C          92 -- INTENSITY
C
 1900 CALL PPBTR(16,INTEN,NOEOF)
      CALL VDSTIN(INTEN/32767.)
      GOTO 3
C
C          93 -- SET BACKGROUND COLOR
C
 2000 CALL PPBTR(16,INDEX,NOEOF)
      NEXTBC = INDEX
      IF(NEXTBC.NE.INT(VECTOR(2))) DRAWN = .TRUE.
      GOTO 3
C
C          A0 -- SET 2D MODE / CLEAR 3D MODE
C
 3300 TWODIM = .TRUE.
      GOTO 3
C
C          A1 -- SET 3D MODE / CLEAR 2D MODE
C
 3400 TWODIM = .FALSE.
      GOTO 3
C
C          A2 -- SET DRAW MODE / CLEAR MARK MODE
C
 3500 MARKER = .FALSE.
      GOTO 3
C
C          A3 -- SET MARK MODE / CLEAR DRAW MODE
C
 3600 MARKER = .TRUE.
      GOTO 3
C
C          A4 -- PLOT MARKER AT CURRENT POSITION
C
 3700 CALL VDIQCP(XPOS,YPOS)
      CALL VDPNTA(XPOS,YPOS)
      GOTO 3
C
C          A8 -- SET LINE STYLE
C
 4100 CALL PPBTR(16,LSTYL,NOEOF)
      IF(LSTYL.EQ.32767) THEN
          CALL VDSTLS(0)
      ELSE IF(LSTYL.EQ.0) THEN
          CALL VDSTLS(1)
      ELSE IF(LSTYL.EQ.16382) THEN
          CALL VDSTLS(2)
      ELSE IF(LSTYL.EQ.5461) THEN
          CALL VDSTLS(3)
      ELSE IF(LSTYL.EQ.27305) THEN
          CALL VDSTLS(4)
      ELSE IF(LSTYL.EQ.16383) THEN
          CALL VDSTLS(5)
      ENDIF
      GOTO 3
C
C          A9 -- SET LINE WIDTH
C
 4200 CALL PPBTR(16,LWID,NOEOF)
      CALL VDSTLW(LWID/32767.)
      GOTO 3
C
C (PLC)
C
C          AA -- POLYGONS
C
 4400 NPTS =0
 4401 CALL PPBTR(16,OPCODE,NOEOF)
      IF(OPCODE.GE.32768)GOTO 4410
        NPTS = NPTS + 1
        XARRAY(NPTS) = FLOAT(OPCODE) * SCALE
        CALL PPBTR(16,OPCODE,NOEOF)
         IF(OPCODE.GE.32768)GOTO 4410
         YARRAY(NPTS) = FLOAT(OPCODE) * SCALE
         GOTO 4401
C
C          IF IT'S AN OPCODE BUT ISN'T = AB (END OF POLYGON COMMAND),
C          THEN IT MUST BE AN OUT OF RANGE COORDINATE.
C
 4410 IF(OPCODE.NE.43776) THEN
         PRINT 10030, OPCODE
10030    FORMAT(1X,'ILLEGAL X,Y COORDINATE -- ',I5)
         GOTO 8
         ENDIF
C
C          MAKE SURE THE COLOR WE WANT HAS BEEN SET UP
C
      IF(COLMAP(NEXTFC).EQ.-1) THEN
         INVMAP(NXTDEX) = NEXTFC
         COLMAP(NEXTFC) = NXTDEX
         CALL VDSTFC(NXTDEX)
         CALL VDSTCO(1,NXTDEX,TABLE(1,NEXTFC),0)
 4411    NXTDEX = NXTDEX + 1
         IF(INVMAP(NXTDEX).NE.-1) GOTO 4411
         ENDIF
C
      CALL VDPOLY(XARRAY,YARRAY,NPTS)
      GOTO 3
C
C          B0 -- TEXT STRING
C
 4900 NCHARS = MIN0(132,COUNT*2)
      IF(NCHARS.EQ.0) GOTO 3
      DO 4910 I=1,NCHARS
         CALL PPBTR(8,TEXT(I),NOEOF)
 4910 CONTINUE
C
C          MAKE SURE THE COLOR WE WANT HAS BEEN SET UP
C
      IF(COLMAP(NEXTFC).EQ.-1) THEN
         INVMAP(NXTDEX) = NEXTFC
         COLMAP(NEXTFC) = NXTDEX
         CALL VDSTFC(NXTDEX)
         CALL VDSTCO(1,NXTDEX,TABLE(1,NEXTFC),0)
 4920    NXTDEX = NXTDEX + 1
         IF(INVMAP(NXTDEX).NE.-1) GOTO 4920
         ENDIF
C
      IF(TEXT(NCHARS).EQ.0) NCHARS = NCHARS - 1
      DRAWN = .TRUE.
      CALL VDTEXT(NCHARS,TEXT)
      COUNT = COUNT - (NCHARS+1)/2
      GOTO 4900
C
C          B2 -- SET CHARACTER SIZE
C
 5100 CALL PPBTR(16,ICHAR,NOEOF)
      CALL VDSTCS(ICHAR/32767./.65)
      CALL PPBTR(16,JUNK,NOEOF)
      GOTO 3
C
C          C0 -- START SEGMENT
C
 6500 CALL PPBTR(16,NAME,NOEOF)
      CALL PPBTR(16,ITYPE,NOEOF)
      ARGS(1) = NAME
      ARGS(2) = ITYPE
      CALL VDESCP(200,2,ARGS)
      GOTO 3
C
C          C1 -- END OF SEGMENT
C
 6600 CALL VDESCP(201,0,0.)
      GOTO 3
C
C          C2 -- DELETE SEGMENT
C
 6700 CALL PPBTR(16,NAME,NOEOF)
      ARGS(1) = NAME
      CALL VDESCP(203,1,ARGS)
      GOTO 3
C
C          C3 -- DELETE ALL SEGMENTS
C
 6800 CALL VDESCP(203,1,0.)
      GOTO 3
C
C          C4 -- RENAME SEGMENT
C
 6900 CALL PPBTR(16,NAME1,NOEOF)
      CALL PPBTR(16,NAME2,NOEOF)
      ARGS(1) = NAME1
      ARGS(2) = NAME2
      CALL VDESCP(202,2,ARGS)
      GOTO 3
C
C          C5 -- SEGMENT ATTRIBUTES
C
 7000 CALL PPBTR(16,NAME,NOEOF)
      CALL PPBTR(16,IVALU,NOEOF)
      ARGS(1) = NAME
      ARGS(2) = IVALU
      CALL VDESCP(204,2,ARGS)
      GOTO 3
C
C          C7 -- BATCH COLOR TABLE UPDATE
C
 7200 BATCH = .TRUE.
      GOTO 3
C
C          C8 -- SEND COLOR TABLE / END BATCH UPDATE
C
 7300 BATCH = .FALSE.
      IF(CSPOT.GT.1) CALL VDSTCO(CSPOT-1,INDICES,RGB,0)
      CSPOT = 1
      GOTO 3
C
C          ILLEGAL OP CODE
C
    8 PRINT 10020,OPCODE
10020 FORMAT(1X,'ILLEGAL OP CODE -- ',I5)
C
C          SCAN FOR NEW FRAME
C
    9 CALL PPBTR(8,OPCODE,NOEOF)
      CALL PPBTR(8, COUNT,NOEOF)
      IF(OPCODE.NE.255) GOTO 9
      IF( COUNT.NE.255) GOTO 9
      CALL PPBTR(8,OPCODE,NOEOF)
      CALL PPBTR(8, COUNT,NOEOF)
      IF(OPCODE.NE.255) GOTO 9
      IF( COUNT.NE.255) GOTO 9
      CALL PPBTR(8,OPCODE,NOEOF)
      CALL PPBTR(8, COUNT,NOEOF)
      IF(OPCODE.NE.255) GOTO 9
      IF( COUNT.NE.255) GOTO 9
      COUNT = 2
      CALL VDNWPG
      GOTO 1
C
      END
C
C
C
      SUBROUTINE PPBTR(IWIDTH,RESULT,EOFOK)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C
C
C PPBTR            - FETCH IWIDTH BITS FROM THE INPUT FILE
C
C P. WATTERBERG    - 25 OCT 81
C K. Cole          - 05 oct 90 added SAVE
C
C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C
C ENTRY CONDITIONS - IWIDTH IS THE NUMBER OF BITS TO GET FROM THE FILE.
C                    EOFOK IS TRUE IF AN END OF FILE CAN BE ENCOUNTERED
C                          AND FALSE OTHERWISE.
C
C CALLS            - CDRUPK, CDRRFS
C
C EXIT CONDITIONS  - RESULT CONTAINS THE NEXT IWIDTH BITS OF THE FILE, RIGHT
C                    JUSTIFIED, ZERO FILLED.  OR, IF AN END OF FILE IS OK,
C                    AND ENCOUNTERED, RESULT CONTAINS AN HEX 86 TO SIGNIFY
C                    END OF FILE.
C
C NARRATIVE        -
C
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C
C
      SAVE
      INTEGER IWIDTH, RESULT
      LOGICAL EOFOK
C
c TLC - made buffer one bigger so cdrupk doesn't go out of memory
      INTEGER BUFFER(10001), POINTR, BUFSIZ, PPREAD
C
      DATA POINTR /    10001 /
      DATA BUFSIZ /        0 /
      DATA PPREAD /       55 /
C
      RESULT = 0
      CALL CDRUPK(BUFFER,POINTR,IBITLC,IWIDTH,RESULT)
      IF(POINTR.GT.BUFSIZ) THEN
c          call cdrpst("vdimet.dat",10)
          CALL CDRRFS(PPREAD,10000,BUFFER,IOSTAT)
          IF(IOSTAT.NE.0) THEN
              BUFSIZ = IABS(IOSTAT)
              POINTR = 1
              IBITLC = 0
              CALL CDRUPK(BUFFER,POINTR,IBITLC,IWIDTH,RESULT)
            ELSE
              IF(EOFOK) THEN
                  RESULT = 134
                ELSE
                  PRINT 10010
10010             FORMAT(1X,'UNEXPECTED END OF FILE, PROGRAM ABORT.')
                  CALL VDFRAM(1)
                  CALL VDTERM
                  STOP
                ENDIF
            ENDIF
        ENDIF
      RETURN
      END
