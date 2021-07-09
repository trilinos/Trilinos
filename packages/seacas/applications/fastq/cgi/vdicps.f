C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

      SUBROUTINE VIINIT(ASPECT,JUSTIF)

C     VDI-PostScript driver - B&W and COLOR versions
C     Adapted for all systems by S.L.Thompson
C     Original code from D.Campbell and J.LONG

C     vdi device numbers are
C                                                   device number
C   black & white, batch, no poly fill                  799.1
C   black & white, interactive, no poly                 799.2
C   black & white, batch, poly fill                     799.3
C   black & white, interactive, poly fill               799.4
C   color, batch                                        799.5
C   color, interactive                                  799.6
C   color, batch, black-white interchange               799.7
C   color, interactive, black-white interchange         799.8
C   color, batch, black background                      799.9
C   color, interactive, black background                799.11

C                                                 last mod 6/20/90 slt

C     Note that there are several parameters to set depending on how
C     the package is to be used. Most are in routine pstsel routine
C     which is called at the first of this routine (viinit.) Two other
c     parameters (xinch,yinch) are set in this routine and vdiqd9.

C     This code is for BOTH color and black & white systems.
C     Flag is set for mode in pstsel.

C     Device can be set with escape call before call to vdinit.
C     Otherwise, code will interactively ask for device type.
C     There is also an escape flag for landscape or portrait format.

C     This deck was generated from a qms driver and still has the
C     qms comments in places.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VIINIT           -Initialize SVDI.  postscript device

C D.L. CAMPBELL    -1-DEC-1986
C J.P. LONG        -9-NOV-1987

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   (postscript)

C ENTRY CONDITIONS -ASPECT = real ratio of X dimension to Y dimension.
C                   Range >0.0.  Default: 0. (device dependent).
C                   JUSTIF = integer justification of NDC space on the
C                   device.  Range 0-9.  Default: 0 (device dependent.)

C CALLS            -VBERRH,VDSTCS,VDSTLW,VIMOVA

C EXIT CONDITIONS  -XNDCMX,YNDCMX = real NDC maximum valid values(as
C                   constrained by ASPECT).
C                   VECTOR = real array of attribute values(all device
C                   dependent except VECTOR(4)=0.0).

C NARRATIVE        -This must be the first SVDI call made.  All
C                   attribute values, the color table, and current
C                   position are set to appropriate defaults for the
C                   device.  All necessary input device initialization
C                   is done.  The screen is cleared or paper advanced
C                   if necessary to guarantee a blank view surface for
C                   drawing on.

C                   ASPECT specifies the ratio of the X dimension to the
C                   Y dimension .  Maximum NDC values (at least one of
C                   which will be 1.0) are computed to give the ASPECT
C                   specified.  The default ASPECT (0.0) is device
C                   dependent and equals the aspect ratio of the
C                   physical device, except for variable aspect devices
C                   (such as drum plotters) which are assigned a default
C                   aspect of 1.0.  The NDC rectangle is scaled until
C                   one dimension fills the corresponding dimension of
C                   the device.

C                   JUSTIF determines where the rectangle is located on
C                   the device as diagrammed below:
C                            ---------
C                           | 7  8  9 |
C                           | 4  5  6 |
C                           | 1  2  3 |
C                            ---------
C                   For example, JUSTIF = 7 indicates NDC space will be
C                   upper left justified on the device.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Set parameters for type of usage.
C     Two settings are coded - one for square field of view
C     and one for full field of view.

C     If VDIQDC is called before vdinit, full field of view is selected.
C     Otherwise, square is used.

C     size of full view
      PARAMETER (XINCHO=10.0)
      PARAMETER (YINCHO=7.5)

C     size of square view
C     PARAMETER (XINCHO=7.5)
C     PARAMETER (YINCHO=7.5)
*- INCLUDE PSTSQUAR
C     size of square view window
C     parameters set to get same size plot as imagen and qms b&w.
C      PARAMETER (XINCHO=7.4412525)
C      PARAMETER (YINCHO=7.4412525)
*-
      COMMON /VCMODR/ XINCH, YINCH

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      REAL ASPECT
      INTEGER JUSTIF
      COMMON /VCVEC1/ IVECT
      INTEGER IVECT

C MAXIMUM VALID NDC VALUES. (DEVICE-INDEPENDENT)
      REAL XNDCMX,YNDCMX
      COMMON /VCNDCM/ XNDCMX,YNDCMX
C CURRENT ATTRIBUTE VALUES. (DEVICE-INDEPENDENT)
      REAL VECTOR(7)
      COMMON /VCATTR/ VECTOR
C     VECTOR(1)=FOREGROUND COLOR
C           (2)=BACKGROUND COLOR
C           (3)=INTENSITY
C           (4)=LINE STYLE
C           (5)=LINE WIDTH
C           (6)=CHARACTER BOX Y
C           (7)=CHARACTER BOX X
C CURRENT POSITION.
      REAL XCP,YCP
      COMMON /VCCRPS/ XCP,YCP

C COMPUTER DEPENDENT COMMON VARIABLES AND CONSTANTS.
      include 'vcpstc.blk'
C SCALE FACTORS FOR NDC TO DC MAPPING. (LXY,HC1)
      REAL XSCALE,YSCALE
      COMMON /VCSCAL/ XSCALE,YSCALE
C DC DIMENSIONS OF OFFSETS AND PICTURE. (LXY,HC1)
      INTEGER XPAD,YPAD,XDEVIC,YDEVIC
      COMMON /VCDDIM/ XPAD,YPAD,XDEVIC,YDEVIC
C JOB ID INFORMATION. (HC1, DIC)
      include 'vcjob.blk'

      COMMON /DEVCAP/ DEV(33)
C ESCAPE FLAGS
C PATNO AND BORDER USED BY VIPOLY FOR FILL PATTERN AND BORDER ON/OFF;
C DEFAULT COMPLETE FILL WITH BORDER. PLC.
      COMMON /VCESCP/ PGFORM,PATNO,BORDER
      INTEGER PGFORM,PATNO,BORDER
      CHARACTER COORD*20,XCOORD*4,YCOORD*4

C     pstmlw controls minimum line width
C     kpstbg controls background coloring
C            = 0,  not colored (white ground from paper)
C            = 1,  colored black
C     kpstci controls black-white interchange (colors 0 & 7 only)
C            = 0,  no interchange
C            = 1,  colors interchanged
      COMMON /VCPSTA/ PSTMLW, KPSTBG, KPSTCI

C     mopoly controls polygon fill =0, on ;  =1, off
C     mocolr controls color =0, on ;  =1, off
      COMMON /VCPSTB/ MOPOLY, MOCOLR

      DIMENSION COLDEF(3)

C     vcpstd variables control what to do with empty frames with
C     command is received to dump data to output
C        kempty=0,  frame is void - do not draw
C              >0,  frame has data - draw it
      COMMON /VCPSTD/ KEMPTY

      integer dummy(1)

      DEV(1) = 0.0
      dev(2) = 1.0
      dev(3) = 1.0
      dev(4) = 1.0
      dev(5) = 15.0
      dev(6) = 2.0
      dev(7) = 0.0
      dev(8) = 0.0
      dev(9) = 0.0
      dev(10) = 0.0

      dev(11) = 0.0
      dev(12) = 0.0
      dev(13) = 0.0
      dev(14) = 0.0
      dev(15) = 7230.0
      dev(16) = 5040.0
      dev(17) = 254.0
      dev(18) = 178.0
      dev(19) = 4.0
      dev(20) = 10.0

      dev(21) = 84.0
      dev(22) = 0.0
      dev(23) = 0.0
      dev(24) = 3.0
      dev(25) = 99999.
      dev(26) = 0.0
      dev(27) = 1.0
      dev(28) = 0.0
      dev(29) = 0.0
      dev(30) = 5000.

      dev(31) = 750.
      dev(32) = 0.0
      dev(33) = 1.0

C SET DEFAULT ATTRIBUTE VALUES.  ALL ARE DEVICE DEPENDENT EXCEPT
C VECTOR(4)=0.0.
C .. following removed since should be in block data...
c      DATA VECTOR /0.,7.,1.,0.,.06255,.01,.0/
C     VECTOR(1)=FOREGROUND COLOR - BLACK
C           (2)=BACKGROUND COLOR - WHITE
C           (3)=INTENSITY        - MAXIMUM
C           (4)=LINE STYLE       - SOLID
C           (5)=LINE WIDTH       - ABOUT 1/72 INCHES
C           (6)=CHARACTER BOX Y  - ABOUT 1/10 INCHES
C           (7)=CHARACTER BOX X  - 5/7 OF BOX-Y

      vector(1) = 0.0
      vector(2) = 7.0
      vector(3) = 1.0
      vector(4) = 0.0
      vector(5) = 0.06255
      vector(6) = 0.01
      vector(7) = 0.0

C PROTECT INPUT PARAMETERS FROM BEING CHANGED.
      ASPEC1=ASPECT
      JUSTI1=JUSTIF
      KEMPTY=0

      PGFORM = 0
      PATNO = 20
      BORDER = 1
      XCP = 0.0
      YCP = 0.0

C CHECK FOR VALID ASPECT.  IF(ASPECT.LT.0.0) THEN CALL VBERRH(721,5),
C AND USE DEFAULT ASPECT.
      IF(ASPECT.LT.0.0) THEN
         CALL VBERRH(721,5)
         ASPEC1=0.0
      END IF

C CHECK FOR VALID JUSTIF.  IF(JUSTIF.LT.0 .OR. JUSTIF.GT.9) THEN
C CALL VBERRH(720,5), AND USE DEFAULT JUSTIF.
      IF(JUSTIF.LT.0.OR.JUSTIF.GT.9) THEN
         CALL VBERRH(720,5)
         JUSTI1=0
      END IF

C SCALE NDC UNITS TO DEVICE UNITS.
C FOR QMS, THE PHYSICAL PLOT SURFACE IS XINCH X YINCH (10.x7.5).
C    DEVICE COORDINATES ARE KEPT IN 1/723 INCH TO GAIN SIMPLICITY
C    IN ASSEMBLING CHARACTER STRINGS WITH FLOATING NUMBERS
C    ACCURATE TO TENTHS OF DEVICE UNITS
C THE DEVICE ASPECT IS XINCH/YINCH
C BUT WE MUST MAP THE ASPECT SPECIFIED IN DC INTO THIS
C ADDRESSABILITY,USING AS MUCH OF THE SPACE AS POSSIBLE.
      XINCH=XINCHO
      YINCH=YINCHO

C     test for rscors post or direct mode. Use 7.5x7.5 for direct
C     and 10.0x7.5 for post

C     if VDIQDC has already been called, we are in post mode;
C     otherwise in direct mode
      CALL VDIQD9(XINCH,YINCH)

C                                  CHECK PAGE FORMAT - IF PORTRAIT,
C                                     THEN SWITCH THINGS AROUND
      IF (PGFORM.EQ.1) THEN
         TEMP=XINCH
         XINCH=YINCH
         YINCH=TEMP
         TEMP=DEV(15)
         DEV(15)=DEV(16)
         DEV(16)=TEMP
         TEMP=DEV(17)
         DEV(17)=DEV(18)
         DEV(18)=TEMP
      ENDIF
      XUNITS=XINCH*723.
      YUNITS=YINCH*723.
      DASPEC=XUNITS/YUNITS

C DEFAULT ASPECT = 1., DEFAULT JUSTIF = 1.
      IF(ASPEC1.EQ.0.) ASPEC1=DASPEC
      IF(JUSTI1.EQ.0) JUSTI1=1

      IF(ASPEC1.GE.DASPEC) THEN

C THEN X DIMENSION IS FILLED.
         XDEVIC=XUNITS
         YDEVIC=XUNITS/ASPEC1
         XPAD=0
C FIGURE HOW MUCH Y PADDING IS NEEDED.
         IF(JUSTI1.LE.3) THEN
            YPAD=0
         ELSE IF(JUSTI1.LE.6) THEN
            YPAD=(YUNITS-YDEVIC)/2
         ELSE
            YPAD=YUNITS-YDEVIC
         END IF
      ELSE

C ELSE Y DIMENSION IS FILLED.
         XDEVIC=YUNITS*ASPEC1
         YDEVIC=YUNITS
         YPAD=0
C FIGURE HOW MUCH X PADDING IS NEEDED.
         IF(JUSTI1.EQ.3.OR.JUSTI1.EQ.6.OR.JUSTI1.EQ.9) THEN
            XPAD=XUNITS-XDEVIC
         ELSE IF(JUSTI1.EQ.2.OR.JUSTI1.EQ.5.OR.JUSTI1.EQ.8) THEN
            XPAD=(XUNITS-XDEVIC)/2
         ELSE
            XPAD=0
         END IF
      END IF

C FIGURE MAXIMUM NDC VALUES XNDCMX AND YNDCMX.
      IF(ASPEC1.GE.DASPEC) THEN
         XNDCMX=MIN(1.,ASPEC1)
         YNDCMX=XNDCMX/ASPEC1
      ELSE
         XNDCMX=ASPEC1
         YNDCMX=1.
      END IF

C SET SCALE FACTORS FOR NDC-TO-DEVICE MAPPING.
      XSCALE=DBLE(XDEVIC)/XNDCMX
      YSCALE=DBLE(YDEVIC)/YNDCMX
      IF (PGFORM .GT. 0) THEN
         YPAD = YPAD+280.
         XPAD = XPAD+360.
      ELSE
         XPAD = XPAD+280.
         YPAD = YPAD-180.
      ENDIF

      CALL PSTSEL(' ')

C  SET UP MONITORING INFORMATION
      CALL VBDEV('V PST   ')
      CALL VDMONI(0)
      IVECT=0

C OPEN OUTPUT FILE
      CALL PSTOFS(KOUTFL)

C INITIALIZE the printer

      CALL PSTINI

      CALL PSTBUF(38, '%%Title: Graphics SVDI PostScript File')
      CALL PSTBUF(0,' ')
      CALL PSTBUF(40, '%%Creator: SNL SEACAS SVDI Driver -- cps')
      CALL PSTBUF(0,' ')
      CALL PSTBUF(24, '%%Orientation: Landscape')
      CALL PSTBUF(0,' ')
      CALL PSTBUF(16, '%%Pages: (atend)')
      CALL PSTBUF(0,' ')
      CALL PSTBUF(28, '%%BoundingBox: 36 30 574 750')
      CALL PSTBUF(0,' ')
      CALL PSTBUF(24, '%%DocumentFonts: Courier')
      CALL PSTBUF(0,' ')
      CALL PSTBUF(13, '%%EndComments')
      CALL PSTBUF(0,' ')
      IF(MOCOLR.EQ.0) THEN
        CALL PSTBUF(35,'% this file contains color commands')
        CALL PSTBUF(0,' ')
C       default user dictionary is too small to contain color
C       definition commands - make a larger one.
        CALL PSTBUF(14,'300 dict begin')
        CALL PSTBUF(0,' ')
      END IF
      CALL PSTBUF(27,'/y {/Courier findfont} def ')
      CALL PSTBUF(27,'/x {scalefont setfont} def ')
      CALL PSTBUF(32,'/m {moveto} def /l {lineto} def ')
      CALL PSTBUF
     *     (50,'/c {closepath} def /v {save} def /r {restore} def ')
      CALL PSTBUF
     *    (54,'/f {eofill} def /s {stroke} def /w {setlinewidth} def ')
      CALL PSTBUF(31,'/h {setdash} def /t {show} def ')
      CALL PSTBUF(33,'/d {gsave} def /u {grestore} def ')
      if (dev(23) .ge. 799.75 .and. dev(23) .le. 799.85) then
         CALL PSTBUF(39,'/q {exch pop exch dup 2 exp 1 exch sub ')
         CALL PSTBUF(39,'setgray add /brsum exch def brsum 0 le ')
         CALL PSTBUF(40,'{ 0 setgray } if brsum 2 ge{ 1 setgray }')
         CALL PSTBUF(10,' if } def ')
      else if (dev(23) .ge. 799.05 .and. dev(23) .le. 799.15) then
         CALL PSTBUF(28,'/q {exch pop exch dup 2 exp ')
         CALL PSTBUF(39,'setgray add /brsum exch def brsum 0 le ')
         CALL PSTBUF(40,'{ 0 setgray } if brsum 2 ge{ 1 setgray }')
         CALL PSTBUF(10,' if } def ')
      else
         CALL PSTBUF(21,'/q {setrgbcolor} def ')
      end if
      CALL PSTBUF(14,'1 setlinejoin ')
      CALL PSTBUF(0,' ')
C                                       SET PAGE FORMAT (LANDSCAPE/PORTRAIT)
       IF (PGFORM.EQ.0) THEN
          CALL PSTBUF(4,'/o {')
          CALL PSTBUF(10,'90 rotate ')
          CALL PSTI2C(0,4,XCOORD)
          CALL PSTI2C(INT(YDEVIC+YDEVIC*3./32.),4,YCOORD)
          COORD = ' '//XCOORD(1:3)//'.'//XCOORD(4:4)//' -'//
     1    YCOORD(1:3)//'.'//YCOORD(4:4)
          CALL PSTBUF( 13,COORD)
          CALL PSTBUF(11,' translate ')
          CALL PSTBUF(6,'} def ')
          YPAD = -YPAD
       ELSE
          CALL PSTBUF(17,'/o {newpath} def ')
       ENDIF
      CALL PSTBUF(35,'/p {showpage} def 1 setlinecap v o ')

C     check for color or black & white mode

      IF(MOCOLR.EQ.0) THEN

C       color is on

C       define some kind of color table

        DO 120 IC=0,7
        COLDEF(1)=0.
        COLDEF(2)=0.
        COLDEF(3)=0.
        IF(IC.EQ.1) THEN
          COLDEF(1)=1.
        ELSEIF(IC.EQ.2) THEN
          COLDEF(2)=1.
        ELSEIF(IC.EQ.3) THEN
          COLDEF(1)=1.
          COLDEF(2)=1.
        ELSEIF(IC.EQ.4) THEN
          COLDEF(3)=1.
        ELSEIF(IC.EQ.5) THEN
          COLDEF(1)=1.
          COLDEF(3)=1.
        ELSEIF(IC.EQ.6) THEN
          COLDEF(2)=1.
          COLDEF(3)=1.
        ELSEIF(IC.EQ.7) THEN
          COLDEF(1)=1.
          COLDEF(2)=1.
          COLDEF(3)=1.
        END IF
          DO 115 IK=0,255,8
          DUMMY(1) = IC+IK
          CALL VDSTCO(1,DUMMY,COLDEF,0)
          IF(IC.EQ.0) THEN
            COLDEF(1)=0.2
            COLDEF(2)=0.2
            COLDEF(3)=0.2
          END IF
  115     CONTINUE
  120   CONTINUE
      END IF
      VECTOR(1)=7.
      VECTOR(2)=0.

C     define the postscript current position
      CALL VBVECT(0,XCP,YCP)

C     shade background if appropriate

      IF(KPSTBG.NE.0) THEN
        CALL PSTBBG
        KEMPTY=0
      END IF

C INIT LINE WIDTH,CHARACTER SIZE
      CALL VDSTLW(VECTOR(5))
c      CALL VDSTCS(VECTOR(6))
      CALL VDSTFC(NINT(VECTOR(1)))
      CALL PSTBUF(0,' ')
      RETURN
      END
      SUBROUTINE VDIQDC(INDEX,VALUE)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDIQDC           -Inquire Device Capabilities.

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77

C ENTRY CONDITIONS -INDEX = integer capability number.  Range 1-33.

C CALLS            -

C EXIT CONDITIONS  -VALUE = real value of the capability indicated by
C                   INDEX.

C NARRATIVE        -Return values of various device capabilities.  INDEX
C                   is the integer capability number (as given below),
C                   and the real value is returned in VALUE.
C     1. Erasability
C        0.  None (hard copy)
C        1.  Screen (Tektronix 4010)
C        2.  Point or Line SLOW (Plasma)
C        3.  Point or Line MEDIUM (Refresh -Serially connected)
C        4.  Point or Line FAST (Refresh -Direct connected)
C        5.  (1) and some (3) (Tektronix 4014 with write-thru mode)
C     2. Scan Type
C        0.  Vector
C        1.  Raster
C        2.  Matrix (Plasma)
C     3. Intensities (1-N)
C     4. Colors (1-N)  This is the number of colors that can be
C                      displayed at one time and may be less than the
C                      total number of colors the device can produce.
C     5. Line Widths (1-N)
C     6. Line Styles (0-N) A bit pattern indicating which of the 5
C                      non-solid line styles are supported in the
C                      device.  Bits 4,3,2,1, and 0 correspond to line
C                      styles medium dash, long dash, short dash, dot
C                      dash, and dot. (0 - device has no hardware line
C                      styles - simulate.)
C     7. Character Sizes (0-N)  (0 - device has no hardware - simulate)
C     8. Number of Locator Devices
C     9. Number of Valuator Devices
C    10. Number of Button Devices
C    11. Number of Keyboard Devices
C    12. Number of Stroke Devices
C    13. Input
C        0. none
C        1. synchronous only - program requests input, then the user
C           supplies it.
C        2. synchronous and asynchronous - synchronous is the same
C           as in (1) above.  Asynchronous means the user can provide
C           input at any time; this input is then saved by the system
C           in an event queue until the program calls for it.
C    14. Input Timing
C        0. no timeout supported
C        1. unreliable timing
C        2. timeout with reliable timing
C    15. X Dimension of View Surface in Device Coordinates
C    16. Y Dimension of View Surface in Device Coordinates
C    17. X Dimension of View Surface in Physical Units (mm) (0 if
C        undefined).  If this dimension is variable (as for drum
C        plotters), it should be set equal to the Y dimension to
C        guarantee a device aspect ratio of 1.0.
C    18. Y Dimension of View Surface in Physical Units (mm) (0 if
C        undefined).
C    19. Smallest Line Width (DC) at default intensity
C    20. Smallest Point (DC) at default intensity
C    21. Smallest Character Size (DC)
C    22. Header and Trailer Frames Required (0=no,1=yes)
C    23. Device Identifier
C        1.  TK4 - Tektronix 4014
C        1.1 TK6 - Tektronix 4016
C        1.2 TEK - Tektronix 4010, 4012
C        1.3 TP2 - Tektronix 4662
C        1.4 TP8 - Tektronix 4662 with 8 pen option
C        1.5 T14 - Tektronix 4114
C        1.6 T13 - Tektronix 4113
C        1.7 T05 - TEKTRONIX 4105
C        1.8 T07 - TEKTRONIX 4107
C        1.9 T15 - TEKTRONIX 4115
C        2.1 16C - Dicomed 16mm color Movies
C        2.2 16B - Dicomed 16mm black and white movies
C        2.3 35C - Dicomed 35mm color slides
C        2.31 3MC - Dicomed 35mm movie color
C        2.4 35B - Dicomed 35mm black and white slides
C        2.41 3MB - Dicomed 35mm movie black and white
C        2.5 35A - Dicomed 35mm Aperture Card
C        2.6 24L - Dicomed 24X Fiche
C        2.7 48L - Dicomed 48X Fiche
C        2.8 CSQ - Dicomed Color Full Frame(square aspect ratio)
C        2.9 BSQ - Dicomed Black and White Full Frame(square aspect)
C        3   R94 - Ramtek 9400
C        4.  T27 - Tektronix 4027
C        4.1 T25 - Tektronix 4025
C        5.  ALP - Alphanumeric Terminal
C        6.  HC1 - Remote Hard Copy
C        7.  LXY - Printronix
C        8.  TST - Test Driver, Print VDI calls made.
C        9.  V25 - Digital VT125
C       10.  AED - AED 512
C       10.1 AE7 - AED 767
C       10.2 AE1 - AED 1024
C       11.  MET - SVDI Metafile
C       12.  HPP - Hewlett Packard Printer/Plotter 2671G
C       12.1 H75 - HP 7580
C       12.2 H72 - HP 7221C OR T
C       12.3 H74 - HP 7475A
C       14.  RET - Retrographics
C       15.  AP5 - Aps 5 Phototypesetter
C       16.  JP7 - JUPITER 7
C       16.1 JP1 - Jupiter 1024
C       17.  GER - Gerber 4400 Photoplotter
C       18.  XYN - XYNETICS
C       20.  PS3 - E & S Picture System 300
C       21.  QMS - QMS LASER PRINTER
C       22.  C51 - CALCOMP 1051 DRUM PLOTTER
C       23.  R25 - RASTER TECHNOLOGIES MODEL ONE/25
C       24.  QLF - QCR large format (8 x 10)
C       24.1 Q35 - QCR 35mm format
C       25.  T45 - Tektronix 4510 Rasterizer
C    24. Polygon support level
C        0.  no support
C        1.  fills convex polygons
C        2.  fills simple polygons (may be concave but not
C            self-crossing)
C        3.  full complex polygon fill support
C    25. Maximum number of points in a polygon (99999. if infinite)
C    26. Setable color table
C        0.  no
C        1.  yes
C    27. Device color palette size (1-N) - the number of different
C            colors a device can produce (may be more than the device
C            can display simultaneously)
C    28. Direct color space size (0-N) - the number of colors
C            available via direct RGB color specification
C            (displayable simultaneously)
C    29. Vector verses Raster VDI
C        0.  SVDI
C        1.  SVDI+Raster
C    30. Maximum character height (DC)
C    31. Maximum line width (DC)
C    32. Color verses monochrome (greyscale) device
C        0.  monochrome
C        1.  color
C    33. Device pixel aspect - the ratio of the spacing of device
C            pixels in x divided by the spacing in y (1 for square
C            pixels)

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      SAVE
      INTEGER INDEX
      REAL VALUE
C ESCAPE FLAGS
C PATNO AND BORDER USED BY VIPOLY FOR FILL PATTERN AND BORDER ON/OFF;
C DEFAULT COMPLETE FILL WITH BORDER. PLC.
      COMMON /VCESCP/ PGFORM,PATNO,BORDER
      INTEGER PGFORM,PATNO,BORDER

C     mopoly controls polygon fill =0, on ;  =1, off
C     mocolr controls color =0, on ;  =1, off
      COMMON /VCPSTB/ MOPOLY, MOCOLR

*- INCLUDE PSTFULL
C     size of full view window
C     parameters set to get same size plot as imagen and qms b&w.
      PARAMETER (XINCHF=9.92167)
      PARAMETER (YINCHF=7.4412525)
*-
C INITIALIZE THE DEVICE CAPABILITIES VECTOR.
      COMMON /DEVCAP/ DEV(33)

      DATA NOCALL /0/

C     If device is 0, call to reset

      IF(NINT(DEV(23)).EQ.0) THEN
        CALL PSTSEL(' ')
      END IF

C CHECK FOR VALID INDEX.
      IF(INDEX.LT.1.OR.INDEX.GT.33) THEN
         CALL VBERRH(726,5)
         GOTO 999
      END IF

C RETURN INDEXED VALUE.
      VALUE=DEV(INDEX)
      IF(INDEX.EQ.23) NOCALL=1

  999 RETURN

C**********************************************************************
      ENTRY VDIQD9(XINCH,YINCH)

C     This is an added entry for rscors version of pst driver to
C     tell if direct or post mode operation. If post mode, vdiqdc
C     is called before vdinit to get terminal type. In direct mode
C     it is never called to get type.

      IF(NOCALL.NE.0) THEN
C       XINCH=10.0
C       YINCH=7.5
        XINCH=XINCHF
        YINCH=YINCHF
      END IF
      RETURN
      END
      SUBROUTINE VBERRH(ERRNUM,ERRSEV)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VBERRH           -Error Handler.

C R.W.Simons       -08APR81
C                   30SEP81

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   Device-independent.

C ENTRY CONDITIONS -ERRNUM = integer error number.
C                   ERRSEV = integer severity code.  If > 12, error is
C                   fatal.

C CALLS            -VDLOGE.

C EXIT CONDITIONS  -

C NARRATIVE        -An error will normally cause an error message to
C                   be printed on the error output device and possible
C                   termination of the program, unless a routine VBERRH
C                   is supplied by the user.  This routine will replace
C                   the default VBERRH provided by the system.  The
C                   system supplied VBERRH calls VDLOGE before
C                   returning.  All versions of VBERRH, whether user-
C                   supplied or default, must STOP on any error severity
C                   greater than 12.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      INTEGER ERRNUM,ERRSEV

C REPORT THE ERROR USING VDLOGE.
      CALL VDLOGE(ERRNUM,ERRSEV)

C CHECK FOR FATAL ERROR.
      IF(ERRSEV.GT.12) STOP

      RETURN
      END
      SUBROUTINE VDGNAM(NAME)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDGNAM           -Name the graphics output file

C P.L.Crotty       -OCT88

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   Device-independent.

C ENTRY CONDITIONS -NAME = character string; < 80 characters

C CALLS

C EXIT CONDITIONS  -output graphics file is assigned the name NAME

C NARRATIVE        -This subroutine associates a file name with
C                   the graphics output file (KOUTFL). If this
C                   routine is not called, a system dependent
C                   default name is used.  VDGNAM must be called
C                   before VDINIT.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      SAVE
      CHARACTER*(*) NAME
      CHARACTER*132 PSTNAM
      INTEGER LENGTH,ISTART,IEND,I
      integer*4 koutff, koutfl

      DATA PSTNAM /'vdicps.ps'/

      DATA ISTAT /0/
      LENGTH = MIN(LEN(NAME),132)

C Strip off any leading blanks
      ISTART = 0
      DO 10 I=1,LENGTH
       IF(NAME(I:I) .NE. ' ')THEN
         ISTART = I
         GOTO 11
       ENDIF
10    CONTINUE
11    CONTINUE

C Strip off trailing blanks
      IEND = 0
      IF(ISTART.GT.0)THEN
        DO 20 I=LENGTH,1,-1
         IF(NAME(I:I) .NE. ' ')THEN
           IEND = I
           GOTO 21
         ENDIF
20      CONTINUE
      ENDIF
21    CONTINUE
      PSTNAM=NAME(ISTART:IEND)
      RETURN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      ENTRY PSTOFS(KOUTFL)
      IF(ISTAT.EQ.0) THEN
        OPEN(KOUTFL,FILE=PSTNAM,FORM='FORMATTED',STATUS='UNKNOWN',
     &  ERR=202,IOSTAT=ISTAT)
        ISTAT=1
      END IF
      GO TO 210
  202 WRITE(6,203) ISTAT,PSTNAM(1:128)
  203 FORMAT(//,' ERROR OPENING PST OUTPUT FILE UNIT =',I8,/,
     &1X,A)
      STOP 'NOOPEN'
  210 CONTINUE
      RETURN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      ENTRY PSTCFS(KOUTFF,KK)
      IF(ISTAT.NE.0) THEN
        CLOSE(KOUTFF,ERR=303)
  303   ISTAT=0
      END IF
      RETURN
      END
      SUBROUTINE VDINIT(ASPECT,JUSTIF)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDINIT           -Initialize SVDI.

C R.W.Simons       -08APR81
C                   30SEP81

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   Device-independent.

C ENTRY CONDITIONS -ASPECT = real ratio of X dimension to Y dimension.
C                   Range >0.0.  Default: 0.0 (device dependent).
C                   JUSTIF = integer justification of NDC space on the
C                   device.  Range 0-9.  Default: 0 (device dependent).

C CALLS            -PSTJOB, VBERRH, VIINIT.

C EXIT CONDITIONS  -XNDCMX,YNDCMX = real NDC maximum valid values.
C                   VECTOR = real array of default attribute values (all
C                   device-dependent except VECTOR(4)=0.0).

C NARRATIVE        -This must be the first SVDI call made.  All
C                   attribute values, the color table, and current
C                   position are set to appropriate defaults for the
C                   device.  All necessary input device initialization
C                   is done.  The screen is cleared or paper advanced
C                   if necessary to guarantee a blank view surface for
C                   drawing.

C                   ASPECT specifies the ratio of the X dimension to the
C                   Y dimension.  Maximum NDC values (at least one of
C                   which will be 1.0) are computed to give the ASPECT
C                   specified.  The default ASPECT (0.0) is device
C                   dependent and equals the aspect ratio of the
C                   physical device, except for variable aspect devices
C                   (such as drum plotters) which are assigned a default
C                   aspect of 1.0.  The NDC rectangle is scaled until
C                   one dimension fills the corresponding dimension of
C                   the device.

C                   JUSTIF determines where the rectangle is located on
C                   the device as diagrammed below:
C                            ---------
C                           | 7  8  9 |
C                           | 4  5  6 |
C                           | 1  2  3 |
C                            ---------
C                   For example, JUSTIF = 7 indicates NDC space will be
C                   upper left justified on the device.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      REAL ASPECT
      INTEGER JUSTIF

C JOB ID INFORMATION. (HC1, DIC)
      include 'vcjob.blk'

C COMPUTER DEPENDENT COMMON VARIABLES AND CONSTANTS.
      include 'vcpstc.blk'

C     set common variables
      KWRTFL=6
      KRDFL=0
      KOUTFL=77
      KINFL=5
      KWRDSZ=0
      KBYTEL=0
      KCPW=0
      KBAUD=0
      KCOMTP=0

C CHECK FOR VALID CLASSIFICATION. Because of output format ignore.
      CALL PSTJOB

C     IF(KSECUR.NE.0) THEN
C        CALL VBERRH(957,13)
C     END IF

C THIS IS JUST A DUMMY ROUTINE WHICH DOES NOTHING BUT CALL VIINIT.
C THIS ORGANIZATION FACILITATES ADDING SECURITY MARKINGS TO SVDI.
      CALL VIINIT(ASPECT,JUSTIF)

      RETURN
      END
      SUBROUTINE VDIQND(XNDC,YNDC)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDIQND           -Inquire NDC Space.

C R.W.Simons       -08APR81

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   Device-independent.

C ENTRY CONDITIONS -XNDCMX,YNDCMX = real maximum valid NDC values.

C CALLS            -

C EXIT CONDITIONS  -XNDC,YNDC = real maximum valid NDC values (XNDCMX,
C                   YNDCMX).

C NARRATIVE        -Return the maximum NDC values as set to realize the
C                   aspect defined by VDINIT.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      REAL XNDC,YNDC

C MAXIMUM VALID NDC VALUES. (DEVICE-INDEPENDENT)
      REAL XNDCMX,YNDCMX
      COMMON /VCNDCM/ XNDCMX,YNDCMX

C RETURN THE MAXIMUM VALID NDC VALUES.
      XNDC=XNDCMX
      YNDC=YNDCMX

      RETURN
      END
      SUBROUTINE VDIQOS(ATTARR)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDIQOS           -Inquire Output Status (of Attributes).

C K.M. ERICKSON    -14 NOV 80

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   Device-independent.

C ENTRY CONDITIONS -VECTOR = real array of current attribute values.

C CALLS            -

C EXIT CONDITIONS  -ATTARR = real array of current attribute value
C                   (VECTOR).

C NARRATIVE        -Return the current attribute values in ATTARR as
C                   given below.
C                   ATTARR(1)=Foreground Color
C                         (2)=Background Color
C                         (3)=Intensity
C                         (4)=Line Style
C                         (5)=Line Width
C                         (6)=Character Box Y
C                         (7)=Character Box X

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      REAL ATTARR(7)

C CURRENT ATTRIBUTE VALUES. (DEVICE-INDEPENDENT)
      REAL VECTOR(7)
      COMMON /VCATTR/ VECTOR
C     VECTOR(1)=FOREGROUND COLOR
C           (2)=BACKGROUND COLOR
C           (3)=INTENSITY
C           (4)=LINE STYLE
C           (5)=LINE WIDTH
C           (6)=CHARACTER BOX Y
C           (7)=CHARACTER BOX X

      INTEGER I

      DO 100 I=1,7
         ATTARR(I)=VECTOR(I)
  100 CONTINUE

      RETURN
      END
      SUBROUTINE VDLINA(X,Y)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDLINA           -Line Absolute.

C R.W.Simons       -08APR81

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   Device-independent.

C ENTRY CONDITIONS -X,Y = real NDC position.

C CALLS            -VILINA.

C EXIT CONDITIONS  -

C NARRATIVE        -Draw a line from current position to absolute NDC
C                   position X,Y and update current position.
C                   Attributes foreground color, intensity, line style,
C                   and line width apply.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      REAL X,Y

C THIS IS JUST A DUMMY ROUTINE WHICH DOES NOTHING BUT CALL VILINA.
C THIS ORGANIZATION FACILITATES ADDING SECURITY MARKINGS TO SVDI.
      CALL VILINA(X,Y)

      RETURN
      END
      SUBROUTINE VDLOGE(ERRNUM,ERRSEV)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDLOGE           -Log Error.

C R.W.Simons       -08APR81
C K.M.Erickson     -8OCT84 - add buffer flush

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   Device-independent.

C ENTRY CONDITIONS -ERRNUM = integer error number.
C                   ERRSEV = integer error severity.

C CALLS            -PSTTBK, VDBUFL

C EXIT CONDITIONS  -

C NARRATIVE        -Report error with message to user and possibly
C                   terminate job depending on severity.  Notice that
C                   by judicious use of the error routines (see VBERRH)
C                   it is possible to write very "nice" error routines
C                   that, for example, only report the first two
C                   occurrences of a particular error, or terminate
C                   if more than 10 errors of a particular severity
C                   occur.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      INTEGER ERRNUM,ERRSEV

C COMPUTER DEPENDENT COMMON VARIABLES AND CONSTANTS.
      include 'vcpstc.blk'

C flush buffer before we do a write
      CALL VDBUFL

C WRITE THE ERROR TO THE LISTING.
      WRITE(KWRTFL,10)ERRNUM,ERRSEV
   10 FORMAT(' SVDI ERROR NUMBER ',I5,'   SEVERITY CODE ',I5)

C TRACEBACK.
csam  CALL PSTTBK

      RETURN
      END
      SUBROUTINE VDMONI(ISTATE)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDMONI           -Logs Usage Information..

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   Device-independent.

C ENTRY CONDITIONS -ISTATE = 0 - initialization
C                            1 - new page
C                            2 - terminate

C CALLS

C EXIT CONDITIONS  -

C NARRATIVE        -For ISTATE=0, job information is initialized, and
C                   timers are initialized called by VIINIT.
C                   ISTATE=1 will increment a common block page
C                   counter called by VINWPG.
C                   ISTATE=2 is called by VITERM and will cause
C                   the collected usage monitoring information to
C                   be written to a file.
C                   Contains entry points VBPKG which will has
C                   an 8 character parameter which will set a common
C                   block variable specifying the package being used.
C                   Entry point VBDEV has an 8 character parameter
C                   which will set a common block variable specifying
C                   the device being used.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C     dummy routine

      CHARACTER*(*) C1,C2

      RETURN
C Usage Monitoring Information

      ENTRY VBPKG (C1)
      RETURN
      ENTRY VBDEV (C2)
      RETURN
      ENTRY VBIQPK(C1)
      C1=' '
      RETURN
      ENTRY VBIQDV(C2)
      C2=' '
      RETURN
      END
      SUBROUTINE VDMOVA(X,Y)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDMOVA           -Move Absolute.

C R.W.Simons       -08APR81

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   Device-independent.

C ENTRY CONDITIONS -X,Y = real NDC position.

C CALLS            -VIMOVA

C EXIT CONDITIONS  -

C NARRATIVE        -Set current position to absolute NDC position X,Y.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      REAL X,Y

C THIS IS JUST A DUMMY ROUTINE WHICH DOES NOTHING BUT CALL VIMOVA.
C THIS ORGANIZATION FACILITATES ADDING SECURITY MARKINGS TO SVDI.
      CALL VIMOVA(X,Y)

      RETURN
      END
      SUBROUTINE VDNWPG
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDNWPG           -New Page.

C R.W.Simons       -08APR81

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   Device-independent.

C ENTRY CONDITIONS -

C CALLS            -VINWPG.

C EXIT CONDITIONS  -

C NARRATIVE        -Physically advance the medium or clear the screen,
C                   whichever is appropriate.  Also flood the screen
C                   with the background color on devices that support
C                   this function. The current position is not changed.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C THIS IS JUST A DUMMY ROUTINE WHICH DOES NOTHING BUT CALL VINWPG.
C THIS ORGANIZATION FACILITATES ADDING SECURITY MARKINGS TO SVDI.
      CALL VINWPG

      RETURN
      END
      SUBROUTINE VDPNTA(X,Y)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDPNTA           -Point Absolute.

C R.W.Simons       -08APR81

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   Device-independent.

C ENTRY CONDITIONS -X,Y = real NDC position.

C CALLS            -VIPNTA.

C EXIT CONDITIONS  -

C NARRATIVE        -Set current position to absolute NDC position X,Y
C                   and put a visible point there.  Attributes
C                   foreground color and intensity apply.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      REAL X,Y

C THIS IS JUST A DUMMY ROUTINE WHICH DOES NOTHING BUT CALL VIPNTA.
C THIS ORGANIZATION FACILITATES ADDING SECURITY MARKINGS TO SVDI.
      CALL VIPNTA(X,Y)

      RETURN
      END
      SUBROUTINE VDPOLY(XARRAY,YARRAY,NPTS)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDPOLY           -POLYGON FILL ROUTINE

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77

C ENTRY CONDITIONS -XARRAY-ARRAY OF X VALUES OF THE POLYGON
C                   YARRAY-CORRESPONDING ARRAY OF Y VALUES
C                   NPTS- NUMBER OF POINTS IN (XARRAY,YARRAY)

C CALLS            -VIPOLY

C EXIT CONDITIONS  -

C NARRATIVE        -The polygon defined by XARRAY,YARRAY will be drawn
C                   and filled (constrained by any limitations of the
C                   physical device -- see below).  No checking will be
C                   done -- all points will be passed to the device.
C                   Current foreground color is used and the polygon
C                   boundary is drawn using the solid line style.
C                   VDI will close the polygon (i.e. the last point
C                   will be connected to the first).

C                   The level of support for this primitive is device-
C                   dependent.  The level of support is categorized
C                   as follows:

C                     Level 0 -- no polygon fill.  Only the polygon
C                        boundary is drawn.
C                     Level 1 -- the device fills convex polygons.
C                     Level 2 -- the device fills simple polygons (may
C                        be concave but not self-crossing)
C                     Level 3 -- full support for complex polygons (may
C                        be self-crossing). In general, the interior of
C                        a complex polygon is defined by the set of points
C                        such that, for each point, when an imaginary line
C                        is drawn to that point from a point far outside
C                        the polygon, that line intersects the polygon
C                        boundary an odd number of times.

C                   Note that the level of support for a particular device
C                   can be inquired using the function VDIQDC.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      INTEGER NPTS
      REAL XARRAY(NPTS),YARRAY(NPTS)

C THIS IS JUST A DUMMY ROUTINE WHICH DOES NOTHING BUT CALL VIPOLY.
C THIS ORGANIZATION FACILITATES ADDING SECURITY MARKINGS TO SVDI.
C     mopoly controls polygon fill =0, on ;  =1, off
C     mocolr controls color =0, on ;  =1, off
      COMMON /VCPSTB/ MOPOLY, MOCOLR

      IF(MOPOLY.EQ.0) THEN
        CALL VIPOLY(XARRAY,YARRAY,NPTS)
      END IF

      RETURN
      END
      SUBROUTINE VDSTOS(ATTARR)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDSTOS           -Set Output Status (of Attributes).

C K.M. ERICKSON    -14 NOV 80

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   Device-independent.

C ENTRY CONDITIONS -ATTARR = real array of attribute values.

C CALLS            -VDSTBC,VDSTCS,VDSTFC,VDSTIN,VDSTLS,VDSTLW

C EXIT CONDITIONS  -VECTOR = real updated attribute values (ATTARR).

C NARRATIVE        -Set the attribute values from ATTARR as given below.
C                   ATTARR(1)=Foreground Color
C                         (2)=Background Color
C                         (3)=Intensity
C                         (4)=Line Style
C                         (5)=Line Width
C                         (6)=Character Box Y

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      REAL ATTARR(6)

C CURRENT ATTRIBUTE VALUES. (DEVICE-INDEPENDENT)
      REAL VECTOR(7)
      COMMON /VCATTR/ VECTOR
C     VECTOR(1)=FOREGROUND COLOR
C           (2)=BACKGROUND COLOR
C           (3)=INTENSITY
C           (4)=LINE STYLE
C           (5)=LINE WIDTH
C           (6)=CHARACTER BOX Y
C           (7)=CHARACTER BOX X

C CALL EACH OF THE INDIVIDUAL ATTRIBUTE SETTING ROUTINES.
C CHECK FOR VALIDITY OF INPUT VALUES WILL BE DONE IN EACH INDIVIDUAL
C ROUTINE.
      CALL VDSTFC(NINT(ATTARR(1)))
      CALL VDSTBC(NINT(ATTARR(2)))
      CALL VDSTIN(ATTARR(3))
      CALL VDSTLS(NINT(ATTARR(4)))
      CALL VDSTLW(ATTARR(5))
c      CALL VDSTCS(ATTARR(6))

      RETURN
      END
      SUBROUTINE VDTERM
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDTERM           -Terminate SVDI.

C R.W.Simons       -08APR81

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   Device-independent.

C ENTRY CONDITIONS -

C CALLS            -VITERM.

C EXIT CONDITIONS  -

C NARRATIVE        -Terminate the SVDI by flushing buffers, etc.  This
C                   should be the last SVDI call made.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C THIS IS JUST A DUMMY ROUTINE WHICH DOES NOTHING BUT CALL VITERM.
C THIS ORGANIZATION FACILITATES ADDING SECURITY MARKINGS TO SVDI.
      CALL VITERM

      RETURN
      END
      SUBROUTINE VDTEXT(LENGTH,CHARS)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDTEXT           -Text from Array.

C R.W.Simons       -08APR81

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   Device-independent.

C ENTRY CONDITIONS -LENGTH = integer number of characters in CHARS.
C                   Range 1-136.
C                   CHARS = integer array of ASCII characters, one
C                   character per array element, right justified.
C                   Range 8,10,32-126.

C CALLS            -VITEXT.

C EXIT CONDITIONS  -

C NARRATIVE        -Draw LENGTH characters in CHARS, starting at
C                   current position and update current position to
C                   the point after the last character box where the
C                   next character would begin.  Current position
C                   indicates the lower left corner of the first
C                   character box.  Only printable characters (32-126
C                   decimal) and backspace and linefeed are allowed.
C                   All values in this range must produce "reasonable"
C                   output; mapping lower case to upper case letters is
C                   considered reasonable.  Attributes foreground color,
C                   background color, intensity, and  character size
C                   apply.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      INTEGER LENGTH,CHARS(136)

C THIS IS JUST A DUMMY ROUTINE WHICH DOES NOTHING BUT CALL VITEXT.
C THIS ORGANIZATION FACILITATES ADDING SECURITY NARKINGS TO SVDI.
      CALL VITEXT(LENGTH,CHARS)

      RETURN
      END
      SUBROUTINE VDFRAM(ITYPE)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDFRAM           - Draw header or trailer frame

C P. Watterberg    - 27 Aug 81

C ENVIRONMENT      -Computer-independent, system-independent, FORTRAN 77

C ENTRY CONDITIONS - ITYPE = 0   for header frame
C                          = 1   for trailer frame

C CALLS            - VIFRAM

C EXIT CONDITIONS  -

C NARRATIVE        - Calls vifram to get time and date from the
C                    system via the computer-dependent routine PSTTOD(entry
C                    point in PSTJOB) and writes it on an identification frame.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      INTEGER ITYPE

      CALL VIFRAM(ITYPE)
      RETURN
      END
      SUBROUTINE VIFRAM(ITYPE)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VIFRAM           - Draw header or trailer frame

C P. Watterberg    - 27 Aug 81

C ENVIRONMENT      -Computer-independent, system-independent, FORTRAN 77

C ENTRY CONDITIONS - ITYPE = 0   for header frame
C                          = 1   for trailer frame

C CALLS            -

C EXIT CONDITIONS  -

C NARRATIVE        -NULL ROUTINE

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      INTEGER ITYPE

      RETURN
      END
      SUBROUTINE VDAABU(BTNNUM)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDAABU           -Await Any Button.

C R.W.Simons       -08APR81

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   All Batch Devices.

C ENTRY CONDITIONS -

C CALLS            -

C EXIT CONDITIONS  -BTNNUM = integer number of the button pressed.
C                   Range 1 to a device dependent maximum which must be
C                   at least 8.

C NARRATIVE        -When a button has been pressed, its integer button
C                   number is returned in BTNNUM.  This function flushes
C                   the button buffer, if any.  This function flushes
C                   the output buffers before doing input.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      INTEGER BTNNUM

C COMPUTER DEPENDENT COMMON VARIABLES AND CONSTANTS.
      include 'vcpstc.blk'

C BATCH DEVICES DON'T NEED TO FLUSH BUFFERS.

      BTNNUM=32

      RETURN
      END
      SUBROUTINE VDABGL(BTNNUM,X,Y)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDABGL           -Await Button, Get Locator.

C R.W.Simons       -08APR81

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   All Batch Devices.

C ENTRY CONDITIONS -

C CALLS            -

C EXIT CONDITIONS  -BTNNUM = integer number of the button pressed.
C                   Range 1 to a device dependent maximum that must be
C                   at least 8.
C                   X,Y = real NDC position of the locator.

C NARRATIVE        -Wait until a button is hit, then return the number
C                   of the button in BTNNUM and the NDC value of the
C                   locator in X,Y.  This function flushes the output
C                   buffers before doing input.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      REAL X,Y
      INTEGER BTNNUM

C COMPUTER DEPENDENT COMMON VARIABLES AND CONSTANTS.
      include 'vcpstc.blk'

C BATCH DEVICES DON'T NEED TO FLUSH BUFFERS.

      BTNNUM=32
      X=0
      Y=0

      RETURN
      END
      SUBROUTINE VDAKGL(CHAR,X,Y)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDAKGL           -Await Keyboard, Get Locator.

C R.W.SIMONS       -02DEC80

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   All Batch Devices.

C ENTRY CONDITIONS -

C CALLS            -

C EXIT CONDITIONS  -CHAR = integer ASCII character input from the
C                   keyboard, right-justified, zero fill.  Range 32-126.
C                   X,Y = real NDC position of the locator.

C NARRATIVE        -Wait until a key is hit, then return the character
C                   entered in CHAR and the NDC value of the locator
C                   in X,Y.  If the character entered does not fall in
C                   the range 32-126, a blank(32) is returned in CHAR.
C                   This function flushes the output buffers before
C                   doing input.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      REAL X,Y
      INTEGER CHAR

C     dummy routine

      CHAR=32
      X=0.
      Y=0.
      RETURN
      END
      SUBROUTINE VDALOC(X,Y)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDALOC           -Await Locator.

C R.W.Simons       -08APR81

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   All Batch Devices.

C ENTRY CONDITIONS -

C CALLS            -

C EXIT CONDITIONS  -X,Y = real NDC position of the locator.

C NARRATIVE        -Wait until the locator is positioned, then return
C                   the NDC value of the locator in X,Y.  The fact that
C                   the locator is positioned can be signaled in a
C                   variety of device dependent ways, such as clicking
C                   the switch in a tablet pen, hitting a button, or
C                   hitting a key on the keyboard.  Any information
C                   contained in this signal is ignored by this
C                   function, as only the locator position is returned.
C                   This function flushes the output buffers before
C                   doing input.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      REAL X,Y

C COMPUTER DEPENDENT COMMON VARIABLES AND CONSTANTS.
      include 'vcpstc.blk'

C BATCH DEVICES DON'T NEED TO FLUSH BUFFERS.

      X=0
      Y=0

      RETURN
      END
      SUBROUTINE VDBELL
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDBELL           -Ring Bell

C R.W.SIMONS       -02DEC80

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   All Batch Devices.

C ENTRY CONDITIONS -

C CALLS            -

C EXIT CONDITIONS  -

C NARRATIVE        -Ring user's bell to get his attention.  This
C                   function is ignored by batch devices.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C THIS FUNCTION IS IGNORED BY BATCH DEVICES.

      RETURN
      END
      SUBROUTINE VDBUFL
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDBUFL           -Buffer Flush.

C R.W.Simons       -19DEC80

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   All Batch Devices.

C ENTRY CONDITIONS -

C CALLS            -

C EXIT CONDITIONS  -

C NARRATIVE        -Assure that the picture is up-to-date by flushing
C                   buffers if necessary.  Also prepare the device to
C                   operate in alphanumeric (as opposed to graphic)
C                   mode.  This is necessary on some devices so that
C                   alphanumeric data from FORTRAN I/O won't be
C                   misinterpreted as graphic data.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C THIS FUNCTION IS IGNORED BY BATCH DEVICES.

      RETURN
      END
      SUBROUTINE VDSTLA(LOCX,LOCY)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDSTLA           -Set Initial Locator Position.

C R.W.Simons       -08APR81

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   All Batch Devices.

C ENTRY CONDITIONS -LOCX,LOCY = real NDC position that the locator is
C                   initilaized to.

C CALLS            -

C EXIT CONDITIONS  -

C NARRATIVE        -Set the initial locator position (light pen tracking
C                   cross, for example) each time this function is
C                   called.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      REAL LOCX,LOCY

C BATCH DEVICES IGNORE THIS FUNCTION.

      RETURN
      END
      SUBROUTINE VDWAIT
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDWAIT           -Wait for User.

C R.W.SIMONS       -02DEC80

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   All Batch Devices.

C ENTRY CONDITIONS -

C CALLS            -

C EXIT CONDITIONS  -

C NARRATIVE        -Wait for the user to view the screen and signal he
C                   is done, normally by hitting any key.  This function
C                   flushes the output buffers before doing input.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C BATCH DEVICES IGNORE THIS COMMAND.

      RETURN
      END
      SUBROUTINE VDIQCO(NUM,INDEX,CLRARY,CLRMOD)
C     C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C
C
C     VDIQCO           -Inquire Color Table.
C
C     R.W.Simons       -08APR81
C     H. S. LAUSON      29MAY86 - changed for current HLS interpretation
C
C     ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C     All Black and White Devices. (LXY, HC1, ALP)
C
C     ENTRY CONDITIONS -NUM = integer number of color indexes to inquire.
C     Range 1-256.
C     INDEX = integer array of indexes to inquire.  Range
C     0-255.
C     CLRMOD = integer color model to be used.  Range 0,1.
C
C     CALLS            -VBERRH
C
C     EXIT CONDITIONS  -CLRARY = real array of 3 by NUM elements returning
C     the values of the components of the indexes inquired.
C     Range for RGB: red 0.0-1.0
C     green 0.0-1.0
C     blue 0.0-1.0
C     Range for HLS: hue 0.0-360.0
C     lightness 0.0-1.0
C     saturation 0.0-1.0
C
C     NARRATIVE        -Inquire one or more color table entries.  NUM and
C     INDEX specify how many and which indexes are being
C     inquired.  CLRMOD specifies which color model
C     (0=RGB, 1=HLS) should be used in constructing values
C     to return in CLRARY.  A device which does not
C     support a color table index specified will
C     return -1.0 in the first element of the CLRARY value
C     for that index.
C
C     C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C
C
      INTEGER NUM,INDEX(NUM),CLRMOD
      REAL CLRARY(3,NUM)
C
      COMMON /PCOLST/ PCOLS(3,256)
C
C     CHECK FOR VALID NUM.
      IF(NUM.LT.1.OR.NUM.GT.256) THEN
         CALL VBERRH(723,5)
         GOTO 999
      END IF
C
C     CHECK FOR VALID CLRMOD.
      IF(CLRMOD.NE.0.AND.CLRMOD.NE.1) THEN
         CALL VBERRH(725,5)
         GOTO 999
      END IF
C
      IF(CLRMOD.NE.0) STOP 'HLS COLORS NOT SUPPORTED'
C
C     CHECK FOR VALID INDEXES.
      DO I=1,NUM
         INDEXN=INDEX(I)
         IF(INDEXN.LT.0.OR.INDEXN.GT.255) THEN
            CALL VBERRH(724,5)
            GOTO 100
         END IF
         CLRARY(1,I)=PCOLS(1,INDEXN)
         CLRARY(2,I)=PCOLS(2,INDEXN)
         CLRARY(3,I)=PCOLS(3,INDEXN)
 100     continue
      end do
C
 999  RETURN
      END
      SUBROUTINE VDIQCP(X,Y)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDIQCP           -Inquire Where Current Position Is.

C R.W.Simons       -02DEC80

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   All Devices that support a software CP.
C                   (AP5,GER,H50,HC1,HCB,HPP,I10,I30,LXY,QCR,QMS,XYN)

C ENTRY CONDITIONS -

C CALLS            -

C EXIT CONDITIONS  -X,Y = real NDC position.

C NARRATIVE        -Return the value of current position.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C CURRENT POSITION.
      REAL XCP,YCP
      COMMON /VCCRPS/ XCP,YCP

C ASSIGN THE CP TO X,Y.
      X=XCP
      Y=YCP

      RETURN
      END
      SUBROUTINE VDSTBC(COLOR)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDSTBC           -Set Background Color.

C R.W.Simons       -05DEC80

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   All Devices with a constant white background. (LXY,
C                   HC1, ALP)

C ENTRY CONDITIONS -COLOR = integer color table index. Range 0-255.
C                   Default: device dependent, in range 0-7.

C CALLS            -VBERRH

C EXIT CONDITIONS  -VECTOR(2) = real updated background color (COLOR).

C NARRATIVE        -Set the background color for following VDNWPG or
C                   TEXT primitives for devices supporting these
C                   features.  For example, many raster devices support
C                   both an overall background color and a background
C                   for character drawing(e.g.,red letters on a green
C                   box).
C                   All devices must support at least a single device
C                   dependent value in the range 0-7.
C                   If an unsupported value is specified, set to
C                   default.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      INTEGER COLOR

C CURRENT ATTRIBUTE VALUES. (DEVICE-INDEPENDENT)
      REAL VECTOR(7)
      COMMON /VCATTR/ VECTOR
C     VECTOR(1)=FOREGROUND COLOR
C           (2)=BACKGROUND COLOR
C           (3)=INTENSITY
C           (4)=LINE STYLE
C           (5)=LINE WIDTH
C           (6)=CHARACTER BOX Y
C           (7)=CHARACTER BOX X

C CHECK FOR VALID COLOR.
      IF(COLOR.LT.0.OR.COLOR.GT.255) THEN
         CALL VBERRH(724,5)
         GOTO 999
      END IF

C ONLY THE SINGLE BACKGROUND COLOR 7 (WHITE) IS SUPPORTED,
C SO NO ACTION IS NECESSARY.

      vector(2) = color
  999 RETURN
      END
      SUBROUTINE VDSTCO(NUM,INDEX,CLRARY,CLRMOD)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDSTCO           -Set Color Table.

C R.W.SIMONS       -02DEC80

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   All Black and White Devices. (LXY, HC1, ALP)

C ENTRY CONDITIONS -NUM = integer number of color indexes to be set.
C                   Range 1-256.
C                   INDEX = integer array of indexes to be set.  Range
C                   0-255.
C                   CLRARY = real array of 3 by NUM elements specifying
C                   the values of the components of the index to be
C                   set.
C                   Range for RGB: red 0.0-1.0
C                                  green 0.0-1.0
C                                  blue 0.0-1.0
C                   Range for HLS: hue 0.0-360.0
C                                  lightness 0.0-1.0
C                                  saturation 0.0-1.0
C                   Default:
C                   Index  Color  RGB Values
C                     0    black   0.,0.,0.
C                     1    red     1.,0.,0.
C                     2    green   0.,1.,0.
C                     3    yellow  1.,1.,0.
C                     4    blue    0.,0.,1.
C                     5    magenta 1.,0.,1.
C                     6    cyan    0.,1.,1.
C                     7    white   1.,1.,1.
C                   CLRMOD = integer color model being used.  Range 0,1.
C                   Default: 0 (RGB).

C CALLS            -VBERRH

C EXIT CONDITIONS  -

C NARRATIVE        -Set one or more color table entries.  This is a
C                   dynamic setting, if the device will support it.
C                   "Dynamic" neans that primitives which have already
C                   been drawn will change their appearance when a
C                   dynamic setting is changed.  INDEX is the
C                   position (or array of positions) in the table
C                   (0-255).  CLRARY is a three-element vector (or 3 by
C                   NUM array) with the fractions (0.-1.) of RGB or
C                   values (0.0-360.0, 0.0-1.0, 0.0-1.0) of HLS.
C                   A translator for HLS to RGB will be used from
C                   GSPC79.  CLRMOD specifies which color model is being
C                   used (0=RGB, 1=HLS).
C                   All devices must support at least a single device
C                   dependent value for each of red, green, and blue and
C                   the corresponding values for hue, lightness, and
C                   saturation.  If unsupported values are specified,
C                   set to the closest supported values.
C                   All devices must support both RGB and HLS color
C                   models.
C                   All devices must support at least a single device
C                   dependent INDEX value in the range 0-7.  If an
C                   unsupported value is specified, it should be ignored.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      INTEGER NUM,INDEX(NUM),CLRMOD
      REAL CLRARY(3,NUM)
      CHARACTER*6 KOLIND
      CHARACTER*20 KOLCOM
      COMMON /VCVEC1/ IVECT
      INTEGER IVECT

C     ARRAY TO CONTAIN COMPLETE COLOR TABLE

      COMMON /PCOLST/ PCOLS(3,256)

C     mopoly controls polygon fill =0, on ;  =1, off
C     mocolr controls color =0, on ;  =1, off
      COMMON /VCPSTB/ MOPOLY, MOCOLR

C CHECK FOR VALID NUM.
      IF(NUM.LT.1.OR.NUM.GT.256) THEN
         CALL VBERRH(723,5)
         GOTO 999
      END IF

C CHECK FOR VALID CLRMOD.
      IF(CLRMOD.NE.0.AND.CLRMOD.NE.1) THEN
         CALL VBERRH(725,5)
         GOTO 999
      END IF

C CHECK FOR VALID INDEXES.
      DO 100 I=1,NUM
         INDEXN=INDEX(I)
         IF(INDEXN.LT.0.OR.INDEXN.GT.255) THEN
            CALL VBERRH(724,5)
            GOTO 100
         END IF
C CHECK FOR VALID CLRARY.
         CLRAR1=CLRARY(1,I)
         CLRAR2=CLRARY(2,I)
         CLRAR3=CLRARY(3,I)
         IF(CLRMOD.EQ.0) THEN
            IF(CLRAR1.LT.0..OR.CLRAR1.GT.1.
     X      .OR.CLRAR2.LT.0..OR.CLRAR2.GT.1.
     X      .OR.CLRAR3.LT.0..OR.CLRAR3.GT.1.) THEN
               CALL VBERRH(727,5)
               GOTO 100
            END IF

C 256 INDEXES ARE SUPPORTED:
              DO IC=1,3
                 PCOLS(IC,INDEXN+1)=CLRARY(IC,I)
              end do

C           define symbol for color reference

            IF(MOCOLR.NE.0) GO TO 390

C           if a set of vectors was in process, issue stroke command
C           to draw them - then start a new path.

            IF(IVECT.NE.0) THEN
              CALL PSTBUF(2,'s ')
              IVECT=0
            END IF
            CALL PSTBUF(0,' ')
            CALL PSTBUF(2,'r ')
            KOLIND='/c'
            IF(INDEXN.LE.9) THEN
              WRITE(KOLIND(3:3),'(I1)',ERR=310) INDEXN
              NNN=4
            ELSEIF(INDEXN.LE.99) THEN
              WRITE(KOLIND(3:4),'(I2)',ERR=310) INDEXN
              NNN=5
            ELSE
              WRITE(KOLIND(3:5),'(I3)',ERR=310) INDEXN
              NNN=6
            END IF
            WRITE(KOLCOM,300,ERR=310) (PCOLS(IC,INDEXN+1),IC=1,3)
  300       FORMAT(F5.3,2F6.3,' q}')
            CALL PSTBUF(NNN+26,KOLIND(1:NNN)//'{'//KOLCOM//' def ')
  310       CONTINUE
C           save and restore can not be in same line - why?
            CALL PSTBUF(0,' ')
            CALL PSTBUF(1,'v')
            CALL PSTBUF(0,' ')
  390       CONTINUE
         ELSE
            IF(CLRAR1.LT.0..OR.CLRAR1.GT.360.
     X      .OR.CLRAR2.LT.0..OR.CLRAR2.GT.1.
     X      .OR.CLRAR3.LT.0..OR.CLRAR3.GT.1.) THEN
               CALL VBERRH(727,5)
               GOTO 100
            END IF

C 256 INDEXES ARE SUPPORTED:
           STOP 'HLS COLORS NOT AVAILABLE'
         END IF
  100 CONTINUE

  999 RETURN
      END
      SUBROUTINE VDSTFC(COLOR)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDSTFC           -Set Foreground Color.

C R.W.Simons       -05DEC80

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   All Devices with a constant black foreground. (LXY,
C                   HC1)

C ENTRY CONDITIONS -COLOR = integer color table index . Range 0-255.
C                   Default is device dependent, in range 0-7.

C CALLS            -VBERRH

C EXIT CONDITIONS  -VECTOR(1) = real updated foreground color (COLOR).

C NARRATIVE        -Set the foreground color index, i.e., the color
C                   table index used for drawing future primitives.
C                   Color is an integer from 0-255 which is used as an
C                   index into the color table (see VDSTCO).
C                   All devices must support at least a single device
C                   dependent value in the range 0-7.
C                   If an unsupported value is specified, set to
C                   default.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      INTEGER COLOR
      CHARACTER*5 KOLIND

C     ARRAY TO CONTAIN COMPLETE COLOR TABLE

      COMMON /PCOLST/ PCOLS(3,256)

C     pstmlw controls minimum line width
C     kpstbg controls background coloring
C            = 0,  not colored (white ground from paper)
C            = 1,  colored black
C     kpstci controls black-white interchange (colors 0 & 7 only)
C            = 0,  no interchange
C            = 1,  colors interchanged
      COMMON /VCPSTA/ PSTMLW, KPSTBG, KPSTCI

C     mopoly controls polygon fill =0, on ;  =1, off
C     mocolr controls color =0, on ;  =1, off
      COMMON /VCPSTB/ MOPOLY, MOCOLR

      COMMON /VCVEC1/ IVECT
      INTEGER IVECT
C CURRENT ATTRIBUTE VALUES. (DEVICE-INDEPENDENT)
      REAL VECTOR(7)
      COMMON /VCATTR/ VECTOR
C     VECTOR(1)=FOREGROUND COLOR
C           (2)=BACKGROUND COLOR
C           (3)=INTENSITY
C           (4)=LINE STYLE
C           (5)=LINE WIDTH
C           (6)=CHARACTER BOX Y
C           (7)=CHARACTER BOX X

C CHECK FOR VALID COLOR.
      IF(COLOR.LT.0.OR.COLOR.GT.255) THEN
         CALL VBERRH(724,5)
         GO TO 999
      END IF

      VECTOR(1)=COLOR
      IF(MOCOLR.EQ.0) THEN

C       draw any vectors in stack before changing colors
        IF(IVECT.NE.0) THEN
          CALL PSTBUF(4,'s r ')
          CALL PSTBUF(0,' ')
          CALL PSTBUF(4,'v o ')
          CALL PSTBUF(0,' ')
          IVECT=0
        END IF

C       code using symbols
        KOLIND='c'
        IF(COLOR.LE.9) THEN
          KOLOR=COLOR
C         test for interchange of colors 0 and 7
          IF(KPSTCI.NE.0) THEN
            IF(KOLOR.EQ.7) THEN
              KOLOR=0
            ELSEIF(KOLOR.EQ.0) THEN
              KOLOR=7
            END IF
          END IF
          WRITE(KOLIND(2:2),'(I1)',ERR=999) KOLOR
          NNN=3
        ELSEIF(COLOR.LE.99) THEN
          WRITE(KOLIND(2:3),'(I2)',ERR=999) COLOR
          NNN=4
        ELSE
          WRITE(KOLIND(2:4),'(I3)',ERR=999) COLOR
          NNN=5
        END IF
        CALL PSTBUF(NNN,KOLIND(1:NNN))

      END IF
  999 RETURN
      END
      SUBROUTINE VDSTIN(INTEN)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDSTIN           -Set Intensity.

C R.W.Simons       -05DEC80

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   All Single Intensity Devices. (LXY, HC1)

C ENTRY CONDITIONS -INTEN = real intensity of the image of an output
C                   primitive.  Range 0.-1.  Default: device dependent.

C CALLS            -

C EXIT CONDITIONS  -VECTOR(3) = real updated intensity (INTEN).

C NARRATIVE        -Set the intensity value indicated for future
C                   primitives.  Intensity is a real value between 0
C                   (not visible) and 1 (maximum).  Intensities are
C                   monotonically increasing within this range.
C                   All devices must support at least a single value:
C                   1.0.  If an unsupported value is specified, set to
C                   the closest supported intensity.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      REAL INTEN

C CURRENT ATTRIBUTE VALUES. (DEVICE-INDEPENDENT)
      REAL VECTOR(7)
      COMMON /VCATTR/ VECTOR
C     VECTOR(1)=FOREGROUND COLOR
C           (2)=BACKGROUND COLOR
C           (3)=INTENSITY
C           (4)=LINE STYLE
C           (5)=LINE WIDTH
C           (6)=CHARACTER BOX Y
C           (7)=CHARACTER BOX X

C CHECK FOR VALID INTEN.
      IF(INTEN.LT.0.0.OR.INTEN.GT.1.0) THEN
         CALL VBERRH(401,5)
         GOTO 999
      END IF

C ONLY THE SINGLE INTENSITY 1.0 (MAXIMUM) IS SUPPORTED,
C SO NO ACTION IS NECESSARY.

  999 RETURN
      END
      SUBROUTINE VITERM
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VITERM           -TERMINATE.

C D.L. CAMPBELL    -1-DEC-1986
C J.P. LONG        -9-NOV-1987

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77

C ENTRY CONDITIONS -

C CALLS            -

C EXIT CONDITIONS  -

C NARRATIVE        -Terminate graphics device.  Close output file.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C COMPUTER DEPENDENT COMMON VARIABLES AND CONSTANTS.
      include 'vcpstc.blk'
      COMMON /VCPAGE/ TOTPAG
      INTEGER TOTPAG
      CHARACTER*10 KPAGE
C     mopoly controls polygon fill =0, on ;  =1, off
C     mocolr controls color =0, on ;  =1, off
      COMMON /VCPSTB/ MOPOLY, MOCOLR

C  put out the last page and restore postscript environment so
C                            nothing is left on the stack
      CALL VINWPG
      CALL PSTBUF(2,'r ')
C FLUSH BUFFER
      CALL PSTBUF(0,' ')
C     write end of data message

      WRITE(KPAGE,'(I10)',ERR=345) TOTPAG
      GO TO 349
  345    KPAGE=' ???'
 349  CONTINUE
      CALL PSTBUF(9, '%%Trailer')
      CALL PSTBUF(0,' ')
      IF(MOCOLR.EQ.0) THEN
        CALL PSTBUF(3,'end')
        CALL PSTBUF(0,' ')
      END IF
      CALL PSTBUF(19,'%%Pages: '//KPAGE)
      CALL PSTBUF(0,' ')
      CALL PSTBUF(5, '%%EOF')
      CALL PSTBUF(0,' ')
C CLOSE OUTPUT FILE
      CALL PSTCFS(KOUTFL,1)
      CALL VDMONI(2)

      RETURN
      END
      SUBROUTINE VIMOVA(X,Y)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VIMOVA           -Move Absolute.

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77

C ENTRY CONDITIONS -X,Y = real NDC position.

C CALLS            -

C EXIT CONDITIONS  -XCP,YCP = real updated current position. (X,Y)

C NARRATIVE        -Set current position to absolute NDC position X,Y.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      REAL X,Y

C move
      CALL VBVECT(0,X,Y)

      RETURN
      END
      SUBROUTINE VIPNTA(X,Y)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VIPNTA           -Point Absolute.

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77

C ENTRY CONDITIONS -X,Y = real NDC position.

C CALLS            -VIMOVA,VILINA

C EXIT CONDITIONS  -

C NARRATIVE        -Set current position to absolute NDC position X,Y
C                   and put a visible point there.  Attributes
C                   foreground color and intensity apply.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      REAL X,Y

      CALL VIMOVA(X,Y)
      CALL VILINA(X,Y)

      RETURN
      END
      SUBROUTINE VIPOLY(XARRAY,YARRAY,NPTS)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VIPOLY           -POLYGON FILL ROUTINE

C D.L. CAMPBELL    -1-DEC-1986
C J.P. LONG        -9-NOV-1987

C ENVIRONMENT      -Fortran77, QMS

C ENTRY CONDITIONS -XARRAY-ARRAY OF X VALUES OF THE POLYGON
C                   YARRAY-CORRESPONDING ARRAY OF Y VALUES
C                   NPTS- NUMBER OF POINTS IN (XARRAY,YARRAY)

C CALLS            -

C EXIT CONDITIONS  -

C NARRATIVE        -The polygon defined by XARRAY,YARRAY will be drawn
C                   and filled (constrained by any limitations of the
C                   physical device -- see below).  No checking will be
C                   done -- all points will be passed to the device.
C                   Current foreground color is used and the polygon
C                   boundary is drawn using the solid line style.
C                   VDI will close the polygon (i.e. the last point
C                   will be connected to the first).

C                   The level of support for this primitive is device-
C                   dependent.  The level of support is categorized
C                   as follows:

C                     Level 0 -- no polygon fill.  Only the polygon
C                        boundary is drawn.
C                     Level 1 -- the device fills convex polygons.
C                     Level 2 -- the device fills simple polygons (may
C                        be concave but not self-crossing)
C                     Level 3 -- full support for complex polygons (may
C                        be self-crossing). In general, the interior of
C                        a complex polygon is defined by the set of points
C                        such that, for each point, when an imaginary line
C                        is drawn to that point from a point far outside
C                        the polygon, that line intersects the polygon
C                        boundary an odd number of times.

C                   Note that the level of support for a particular device
C                   can be inquired using the function VDIQDC.

********************************************************************************

C                   The level for this device is level 2.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      REAL XARRAY(NPTS),YARRAY(NPTS)

C ESCAPE FLAGS
C PATNO AND BORDER USED BY VIPOLY FOR PATTERN FILL AND BORDER ON/OFF. DEFAULT
C COMPLETE FILL AND BORDER ON
C CURRENT ATTRIBUTE VALUES. (DEVICE-INDEPENDENT)
      REAL VECTOR(7)
      COMMON /VCATTR/ VECTOR
C     VECTOR(1)=FOREGROUND COLOR
C           (2)=BACKGROUND COLOR
C           (3)=INTENSITY
C           (4)=LINE STYLE
C           (5)=LINE WIDTH
C           (6)=CHARACTER BOX Y
C           (7)=CHARACTER BOX X
      COMMON /VCESCP/ PGFORM,PATNO,BORDER
      INTEGER PGFORM,PATNO,BORDER
      COMMON /VCVEC1/ IVECT
      COMMON /VCVEC2/ COORD, LSTCRD
      CHARACTER COORD*20, LSTCRD*20
      INTEGER IVECT

C     mopoly controls polygon fill =0, on ;  =1, off
C     mocolr controls color =0, on ;  =1, off
      COMMON /VCPSTB/ MOPOLY, MOCOLR

C CHECK FOR VALID N
      IF (NPTS.LT.1 .OR. NPTS.GT.1490) THEN
         CALL VBERRH(802,5)
         GO TO 999
      END IF

C IF A SET OF VECTORS WAS IN PROCESS, ISSUE STROKE COMMAND TO DRAW THEM
C Start a new path.

      IF(IVECT.NE.0) THEN
        CALL PSTBUF(2,'s ')
        IVECT=0
      END IF

      CALL PSTBUF(2,'r ')
      CALL PSTBUF(0,' ')
      CALL PSTBUF(4,'v o ')
      WSAVE=VECTOR(5)
      CALL PSTBUF(4,'0 w ')
C     CALL VDSTLW(0.)
      IF(MOCOLR.EQ.0) THEN
        CALL VDSTFC(NINT(VECTOR(1)))
      END IF
      CALL PSTBUF(0,' ')

C DRAW POLYGON VECTORS

C MOVE TO FIRST POINT
      CALL VIMOVA(XARRAY(1),YARRAY(1))

C CALL VDLINA TO DRAW POINTS FROM 1ST POINT TO NTH POINT
      DO 100 I=2,NPTS
        CALL VILINA(XARRAY(I),YARRAY(I))
100   CONTINUE

C THEN DRAW A LINE TO THE FIRST POINT TO CLOSE THE POLYGON
      CALL VILINA(XARRAY(1),YARRAY(1))

C CLOSE THE POLYGON, GRAPHICS SAVE, FILL IT, GRAPHICS RESTORE, STROKE
C    TO PROVIDE THE SAME FILLED AREA AS IF IT WERE FILLED WITH VECTORS
C    THEN RESTORE AND SAVE POSTSCRIPT ENVIRONMENT TO AVOID INPUT BUFFER OVERFLOW
      CALL PSTBUF(12,'c d f u s r ')
      IVECT=0
      CALL PSTBUF(0,' ')
      CALL PSTBUF(4,'v o ')
      CALL VDSTLW(WSAVE)
C ... if color is on (mocolr = 0), then font is set in vdstfc
      IF(MOCOLR.EQ.0) THEN
        CALL VDSTFC(NINT(VECTOR(1)))
      ELSE
c         CALL VDSTCS(VECTOR(6))
      END IF
      CALL PSTBUF(0,' ')

C INIT THE CURRENT POSITION WITHIN POSTSCRIPT
      CALL VDMOVA(XARRAY(NPTS),YARRAY(NPTS))
      IVECT=0

999   RETURN
      END
      SUBROUTINE VINWPG
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VINWPG           -New Page.

C D.L. CAMPBELL    -1-DEC-1986
C J.P. LONG        -9-NOV-1987

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77

C ENTRY CONDITIONS -

C CALLS            -

C EXIT CONDITIONS  -

C NARRATIVE        -Physically advance the medium or clear the screen,
C                   whichever is appropriate.  Also flood the screen
C                   with the background color on devices that support
C                   this.  The current position is not changed.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      COMMON /VCVEC1/ IVECT
      COMMON /VCVEC2/ COORD,LSTCRD

C     vcpstd variables control what to do with empty frames with
C     command is received to dump data to output
C        kempty=0,  frame is void - do not draw
C              >0,  frame has data - draw it
      COMMON /VCPSTD/ KEMPTY

      CHARACTER COORD*20, LSTCRD*20
      CHARACTER*10 KPAGE

C     pstmlw controls minimum line width
C     kpstbg controls background coloring
C            = 0,  not colored (white ground from paper)
C            = 1,  colored black
C     kpstci controls black-white interchange (colors 0 & 7 only)
C            = 0,  no interchange
C            = 1,  colors interchanged
      COMMON /VCPSTA/ PSTMLW, KPSTBG, KPSTCI
      INTEGER IVECT
C CURRENT ATTRIBUTE VALUES. (DEVICE-INDEPENDENT)
      REAL VECTOR(7)
      COMMON /VCATTR/ VECTOR
C     VECTOR(1)=FOREGROUND COLOR
C           (2)=BACKGROUND COLOR
C           (3)=INTENSITY
C           (4)=LINE STYLE
C           (5)=LINE WIDTH
C           (6)=CHARACTER BOX Y
C           (7)=CHARACTER BOX X

      COMMON /VCPAGE/ TOTPAG
      INTEGER TOTPAG
      DATA NPAGE /0/

C        check for void page draw request
C        if nothing is on page, skip request

      NPAGE=NPAGE+1
      TOTPAG = NPAGE
      WRITE(KPAGE,'(I10)',ERR=345) NPAGE
      GO TO 349
  345    KPAGE=' ???'
  349 IF(KEMPTY.EQ.0) GO TO 350

C  stroke the path in case there are any vectors and show text
      CALL PSTBUF(2,'s ')
      IVECT=0

C   showpage and restore postscript environment to avoid buffer overflow
C            flush buffer because save and restore won't work back-to-back

      CALL PSTBUF(4,'p r ')
      CALL PSTBUF(0,' ')

C     comment frame number in output file

      CALL PSTBUF(31,'%%Page: "'//KPAGE//'" '//KPAGE)
      CALL PSTBUF(0,' ')
      CALL PSTBUF(28, '%%PageOrientation: Landscape')
      CALL PSTBUF(0,' ')
      CALL PSTBUF(32, '%%PageBoundingBox: 36 30 574 750')
      CALL PSTBUF(0,' ')

      CALL PSTBUF(4,'v o ')
      CALL VDMONI(1)

C     shade background if appropriate

      IF(KPSTBG.NE.0) THEN
        CALL PSTBBG
      END IF
      GO TO 370

C     void frame -- First Page

  350 CALL PSTBUF(2, 'r ')
      CALL PSTBUF(0,' ')
      CALL PSTBUF(31,'%%Page: "'//KPAGE//'" '//KPAGE)
      CALL PSTBUF(0,' ')
      CALL PSTBUF(28, '%%PageOrientation: Landscape')
      CALL PSTBUF(0,' ')
      CALL PSTBUF(32, '%%PageBoundingBox: 36 30 574 750')
      CALL PSTBUF(0,' ')
      CALL PSTBUF(4, 'v o ')

  370 CALL VDSTLW(VECTOR(5))
c      CALL VDSTCS(VECTOR(6))
      CALL VDSTFC(NINT(VECTOR(1)))
      CALL PSTBUF(0,' ')
      KEMPTY=0

      RETURN
      END
      SUBROUTINE VDESCP(ESCPCD,N,ARGS)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDESCP           -Escape Code Routine.

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77

C ENTRY CONDITIONS -ESCPCD = integer escape function code.
C                   N = integer number of arguments in ARG.  RANGE >=0.
C                   ARGS = real array of arguments for the escape
C                   function specified.

C CALLS            -

C EXIT CONDITIONS  -

C NARRATIVE        -Invoke the nonstandard, device-dependent function
C                   ESCPCD.  N is the number of arguments used by this
C                   function and ARGS is a real array containing those
C                   arguments.  Unsupported values of ESCPCD are
C                   ignored, not causing an error.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      INTEGER ESCPCD,N
      REAL ARGS(*)

C COMPUTER DEPENDENT COMMON VARIABLES AND CONSTANTS.
      include 'vcpstc.blk'
C SCALE FACTORS FOR NDC TO DC MAPPING. (LXY,HC1)
      REAL XSCALE,YSCALE
      COMMON /VCSCAL/ XSCALE,YSCALE
C DC DIMENSIONS OF OFFSETS AND PICTURE. (LXY,HC1)
      INTEGER XPAD,YPAD,XDEVIC,YDEVIC
      COMMON /VCDDIM/ XPAD,YPAD,XDEVIC,YDEVIC

C USED BY VIPOLY FOR PATTERN FILL AND BORDER ON/OFF. DEFAULT COMPLETE FILL
C AND BORDER ON. PLC.
      COMMON/VCESCP/PGFORM,PATNO,BORDER
      INTEGER PGFORM,PATNO,BORDER

C CHECK FOR VALID N.
      IF(N.LT.0) THEN
         CALL VBERRH(802,5)
         GOTO 999
      END IF

C 2100 - PAGE FORMAT (0=LANDSCAPE,1=PORTRAIT)
      IF (ESCPCD.EQ.2100) THEN
         IF (ARGS(1).EQ.0) THEN
            PGFORM=0
         ELSE
            PGFORM=1
         ENDIF

C     set output format

      ELSEIF (ESCPCD.EQ.2101) THEN
         CALL PSTSEL('1')
      ELSEIF (ESCPCD.EQ.2102) THEN
         CALL PSTSEL('2')
      ELSEIF (ESCPCD.EQ.2103) THEN
         CALL PSTSEL('3')
      ELSEIF (ESCPCD.EQ.2104) THEN
         CALL PSTSEL('4')
      ELSEIF (ESCPCD.EQ.2105) THEN
         CALL PSTSEL('5')
      ELSEIF (ESCPCD.EQ.2106) THEN
         CALL PSTSEL('6')
      ELSEIF (ESCPCD.EQ.2107) THEN
         CALL PSTSEL('7')
      ELSEIF (ESCPCD.EQ.2108) THEN
         CALL PSTSEL('8')
      ELSEIF (ESCPCD.EQ.2109) THEN
         CALL PSTSEL('9')
      ELSEIF (ESCPCD.EQ.2110) THEN
         CALL PSTSEL('10')
      ENDIF

 999  RETURN
      END
      SUBROUTINE VILINA (X,Y)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VILINA

C D.L. CAMPBELL    -1-DEC-1986
C J.P. LONG        -9-NOV-1987

C ENVIRONMENT      -DEVICE DEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77

C ENTRY CONDITIONS -X,Y REAL NDC COORDINATES

C CALLS            -

C EXIT CONDITIONS  -CURRENT POSITION IS SET

C NARRATIVE
C                   LINE-DRAW A LINE FROM CP TO ABSOLUTE NDC POSITION X,Y
C                        AND UPDATE CP . ATTRIBUTES COLOR,INTEN,LINSTY AND
C                        LINWTH APPLY.

C        OTHER VARIABLES:
C                    XCP,YCP-NDC COORDINATES
C***************************************************************************

      REAL X,Y

C     vcpstd variables control what to do with empty frames with
C     command is received to dump data to output
C        kempty=0,  frame is void - do not draw
C              >0,  frame has data - draw it
      COMMON /VCPSTD/ KEMPTY

C draw
      ENTRY VBLINA(X,Y)
      CALL VBVECT(1,X,Y)
      KEMPTY=1

      RETURN
      END
      SUBROUTINE VBVECT(IPEN,X,Y)
C****************************************************
C vbvect - do move or draw to x,y (depending on ipen)

C     ipen = 0 for move, 1 for draw
C     x,y = NDC coordinates to be moved/drawn to

C******************************************************

      REAL X,Y,XOFF,YOFF
      CHARACTER CTEMP*20,XCOORD*4,YCOORD*4

C     mopoly controls polygon fill =0, on ;  =1, off
C     mocolr controls color =0, on ;  =1, off
      COMMON /VCPSTB/ MOPOLY, MOCOLR

C SCALE FACTORS FOR NDC TO DC MAPPING. (LXY,HC1)
      REAL XSCALE,YSCALE
      COMMON /VCSCAL/ XSCALE,YSCALE
C DC DIMENSIONS OF OFFSETS AND PICTURE. (LXY,HC1)
      INTEGER XPAD,YPAD,XDEVIC,YDEVIC
      COMMON /VCDDIM/ XPAD,YPAD,XDEVIC,YDEVIC
C CURRENT POSITION.
      REAL XCP,YCP
      COMMON /VCCRPS/ XCP,YCP
C CURRENT ATTRIBUTE VALUES. (DEVICE-INDEPENDENT)
      REAL VECTOR(7)
      COMMON /VCATTR/ VECTOR
C     VECTOR(1)=FOREGROUND COLOR
C           (2)=BACKGROUND COLOR
C           (3)=INTENSITY
C           (4)=LINE STYLE
C           (5)=LINE WIDTH
C           (6)=CHARACTER BOX Y
C           (7)=CHARACTER BOX X
      COMMON /VCVEC1/ IVECT
      COMMON /VCVEC2/ COORD, LSTCRD
      CHARACTER COORD*20, LSTCRD*20
      INTEGER IVECT

C compute new point in dev. coord.
C     convert to floating offsets
      XOFF=XPAD
      YOFF=YPAD

      IXDC=X*XSCALE+XOFF
      IYDC=Y*YSCALE+YOFF

C        write(xcoord,'(i5)')ixdc
C        write(ycoord,'(i5)')iydc
C                                ...include both x,y
      CALL PSTI2C(IXDC,4,XCOORD)
      CALL PSTI2C(IYDC,4,YCOORD)
      COORD = XCOORD(1:3)//'.'//XCOORD(4:4)//' '//
     1 YCOORD(1:3)//'.'//YCOORD(4:4)

C pack up move/draw command, send it down
C      if (lstcrd(1:11) .ne. coord(1:11)) then
         IF (IPEN.EQ.0) THEN
            CTEMP= COORD(1:11) // ' m '
         ELSE
            CTEMP= COORD(1:11) // ' l '
         ENDIF
         CALL PSTBUF(14,CTEMP)
C     ...count the coordinate pair
         IVECT=IVECT+1
C      end if
      lstcrd(1:11) = coord(1:11)

C  stroke the path if we are approaching the 1500-coord pair limit
C                also restore and save postscript environment to avoid
C                input buffer overflow (must have a c/r between restore
C                and save)
      IF(IVECT.GT.1400) THEN
         lstcrd(1:11) = '           '
         CALL PSTBUF(4,'s r ')
         CALL PSTBUF(0,' ')
         CALL PSTBUF(4,'v o ')
         CALL VDSTLW(VECTOR(5))
         IF(MOCOLR.EQ.0) THEN
           CALL VDSTFC(NINT(VECTOR(1)))
         END IF
C        ...reset the vector count - vdstls (called by vdstlw)
C                 reinitted the current posn
         IVECT=1
      ENDIF

C UPDATE CURRENT POSITION
      XCP=X
      YCP=Y

      RETURN
      END
      SUBROUTINE VITEXT(LENGT1,CHARS)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VITEXT           - Text from Array.

C P. Watterberg    - 24 MAR 81
C J. P. LONG       -  3 DEC 87

C ENVIRONMENT      - COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77

C ENTRY CONDITIONS - LENGT1 = integer number of characters in CHARS.
C                    Range 1-136.
C                    CHARS = integer array of ASCII characters, one
C                    character per array element, right justified.
C                    Range 8,10,32-126.

C CALLS            - vbout

C EXIT CONDITIONS  - XCP,YCP = integer updated current position (at the end
C                    of the string).

C NARRATIVE        - Draw LENGT1 characters in CHARS, starting at
C                    current position and update current position to
C                    the point after the last character box where the
C                    next character would begin.  Current position
C                    indicates the lower left corner of the first
C                    character box.  Only printable characters (32-126
C                    decimal) and backspace and linefeed are allowed.
C                    All values in this range must produce "reasonable"
C                    output; mapping lower case to upper case letters is
C                    considered reasonable.  Attributes foreground color,
C                    background color, intensity, and character size
C                    apply.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      INTEGER LENGT1, CHARS(136), LENGTH

      CHARACTER CTEMP*150,STR*3
C CURRENT POSITION.
      REAL XCP,YCP
      COMMON /VCCRPS/ XCP,YCP
C CURRENT ATTRIBUTE VALUES. (DEVICE-INDEPENDENT)
      REAL VECTOR(7)
      COMMON /VCATTR/ VECTOR
C     VECTOR(1)=FOREGROUND COLOR
C           (2)=BACKGROUND COLOR
C           (3)=INTENSITY
C           (4)=LINE STYLE
C           (5)=LINE WIDTH
C           (6)=CHARACTER BOX Y
C           (7)=CHARACTER BOX X

C SCALE FACTORS FOR NDC TO DC MAPPING. (LXY,HC1)
      REAL XSCALE,YSCALE
      COMMON /VCSCAL/ XSCALE,YSCALE
C DC DIMENSIONS OF OFFSETS AND PICTURE. (LXY,HC1)
      INTEGER XPAD,YPAD,XDEVIC,YDEVIC
      COMMON /VCDDIM/ XPAD,YPAD,XDEVIC,YDEVIC

C     vcpstd variables control what to do with empty frames with
C     command is received to dump data to output
C        kempty=0,  frame is void - do not draw
C              >0,  frame has data - draw it
      COMMON /VCPSTD/ KEMPTY

C          check for valid length.

      call vdstcs(vector(6))
      KEMPTY=1
      LENGTH = LENGT1
      IF(LENGTH.LT.1) THEN
         CALL VBERRH(212,5)
         GO TO 999
      END IF

C          if(length.gt.136) then call vberrh(213,5), and use the
C          maximum length of 136.

      IF(LENGTH.GT.136) THEN
         CALL VBERRH(213,5)
         LENGTH = 136
      ENDIF

      CTEMP='('
      LENOUT=1

C          loop through length characters.

      DO 100 I=1,LENGTH

C          check for valid chars.

C          ignore control characters, except for:
C          8 is backspace
C          10 is linefeed
C          13 is carriage return

         IF(CHARS(I).LT.32 .OR. CHARS(I).GT.126) THEN

            IF(CHARS(I).EQ.8) THEN
               DX=-VECTOR(7)
               DY=0.
            ELSE IF(CHARS(I).EQ.10) THEN
               DX=0.
               DY= -VECTOR(6)
            ELSE IF(CHARS(I).EQ.13) THEN
               DX=XPAD-XCP
               DY=0.
            ELSE
               DX=0.
               DY=0.
               CALL VBERRH(208,5)
               GOTO 100
            ENDIF

C           finish the string, emulate the control char, and start a new one

C           send the buffered chars to the printer if there are any
            IF(LENOUT.NE.1) THEN
               CTEMP(LENOUT+1:150)=') t '
               LENOUT=LENOUT+4
               CALL PSTBUF(LENOUT,CTEMP)
C              reset the cp from the characters
               XCP=XCP+(LENOUT-5)*VECTOR(7)
            ENDIF

C           calculate the new current position after the control char
            XCP=XCP+DX
            YCP=YCP+DY
            CALL VBVECT(0,XCP,YCP)

C           start a new string
            CTEMP='('
            LENOUT=1

         ELSE

C           Char value is 32-126 inclusive.  Put \ before these:
C              92 is \
C              40 is (
C              41 is )

            IF(CHARS(I).EQ.40.OR.CHARS(I).EQ.41.OR.CHARS(I).EQ.92) THEN
               CTEMP(LENOUT+1:150)='\\'
               LENOUT=LENOUT+1
            ENDIF

C           now pack the chars into the buffer

            CALL PSTA2C(CHARS(I),STR)
            CTEMP(LENOUT+1:150)=STR(1:1)
            LENOUT=LENOUT+1
         ENDIF

  100 CONTINUE

C          send the chars to the printer

      CTEMP(LENOUT+1:150)=') t '
      LENOUT=LENOUT+4
      CALL PSTBUF(LENOUT,CTEMP)

C          reset the cp from the characters

      XCP=XCP+(LENOUT-5)*VECTOR(7)

  999 RETURN
      END
      SUBROUTINE VDSTLS(LINSTY)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDSTLS           -Set Line Style.

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77

C ENTRY CONDITIONS -LINSTY = integer linestyle of line drawing output
C                   primitives.  Range 0-5.  Default:0.

C CALLS            -

C EXIT CONDITIONS  -VECTOR(4) = real updated line style (LINSTY).

C NARRATIVE        -Set the style of line as below.  This applies only
C                   to line drawing primitives.  The line styles are:
C                          0 - solid
C                          1 - dotted
C                          2 - dot dash
C                          3 - short dash
C                          4 - long dash
C                          5 - medium dash
C                   All devices must support at least the values 0 and
C                   5.  If an unsupported value is specified, set to 5.
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      REAL LW
      INTEGER LINSTY,ILL,JLL
      COMMON /VCVEC1/ IVECT
      COMMON /VCVEC2/ COORD,LSTCRD
      CHARACTER COORD*20, LSTCRD*20
      INTEGER IVECT
      CHARACTER CTEMP*30,STRL*3,STRS*3,STRG*3

C CURRENT ATTRIBUTE VALUES. (DEVICE-INDEPENDENT)
      REAL VECTOR(7)
      COMMON /VCATTR/ VECTOR
C     VECTOR(1)=FOREGROUND COLOR
C           (2)=BACKGROUND COLOR
C           (3)=INTENSITY
C           (4)=LINE STYLE
C           (5)=LINE WIDTH
C           (6)=CHARACTER BOX Y
C           (7)=CHARACTER BOX X
C SCALE FACTORS FOR NDC TO DC MAPPING. (LXY,HC1)
      REAL XSCALE,YSCALE
      COMMON /VCSCAL/ XSCALE,YSCALE
C CURRENT POSITION.
      REAL XCP,YCP
      COMMON /VCCRPS/ XCP,YCP
C CURRENT ATTRIBUTE VALUES. (DEVICE-INDEPENDENT)
C      REAL VECTOR(7)
C     COMMON /VCATTR/ VECTOR

C     VECTOR(1)=FOREGROUND COLOR
C           (2)=BACKGROUND COLOR
C           (3)=INTENSITY
C           (4)=LINE STYLE
C           (5)=LINE WIDTH
C           (6)=CHARACTER BOX Y
C           (7)=CHARACTER BOX X

      ENTRY VBSTLS(LINSTY)

C CHECK FOR VALID LINSTY.
      IF(LINSTY.LT.0.OR.LINSTY.GT.5) THEN
         CALL VBERRH(401,5)
         VECTOR(4) = 0
         GOTO 999
      END IF

      IF(IVECT.NE.0) THEN
        CALL PSTBUF(2,'s ')
        IVECT=0
      END IF
C GENERATE THE LINESTYLE COMMANDS
      IF(LINSTY.EQ.0) THEN
        CALL PSTBUF(7,'[] 0 h ')
      ENDIF

C calculate the linewidth -- it's needed below in every case

C        actual xscale is xscale*.1; linewidth=1 => .01 in NDC
         LW=VECTOR(5)
         LW=XSCALE*VECTOR(5)*.001
C        a linewidth of zero isn't good with postscript
         IF(LW.LT.1.) LW=1.

C     from here on, set up patterns that depend on the linewidth and
C          the extra length added to the line segment
C          by the hemispherical end cap

      IF(LINSTY.EQ.1) THEN
         ILL=NINT(0.5*LW)
         IGAP=NINT(3.*LW)
         CALL PSTI2C(ILL,3,STRL)
         CALL PSTI2C(IGAP,3,STRG)
         CTEMP='['//STRL(1:3)//' '//STRG(1:3)//'] 0 h '
         CALL PSTBUF(14,CTEMP)

      ELSE IF(LINSTY.EQ.2) THEN
         ILL=NINT(18.*LW)
         JLL=NINT(1.5*LW)
         IGAP=NINT(3.*LW)
         CALL PSTI2C(ILL,3,STRL)
         CALL PSTI2C(JLL,3,STRS)
         CALL PSTI2C(IGAP,3,STRG)
         CTEMP='['//STRS(1:3)//' '//STRG(1:3)//' '//STRL(1:3)
     *            //' '//STRG(1:3)//'] 0 h '
         CALL PSTBUF(22,CTEMP)
C         call pstbuf(14,'[2 2 6 2] 0 h ')

      ELSE IF(LINSTY.EQ.3) THEN
         ILL=NINT(6.*LW)
         IGAP=NINT(7.*LW)
         CALL PSTI2C(ILL,3,STRL)
         CALL PSTI2C(IGAP,3,STRG)
         CTEMP='['//STRL(1:3)//' '//STRG(1:3)//'] 0 h '
         CALL PSTBUF(14,CTEMP)
C         call pstbuf(8,'[4] 0 h ')

      ELSE IF(LINSTY.EQ.4) THEN
         ILL=NINT(24.*LW)
         IGAP=NINT(18.*LW)
         CALL PSTI2C(ILL,3,STRL)
         CALL PSTI2C(IGAP,3,STRG)
         CTEMP='['//STRL(1:3)//' '//STRG(1:3)//'] 0 h '
         CALL PSTBUF(14,CTEMP)
C         call pstbuf(8,'[8] 0 h ')

      ELSE IF(LINSTY.EQ.5) THEN
         ILL=NINT(12.*LW)
         IGAP=NINT(10.*LW)
         CALL PSTI2C(ILL,3,STRL)
         CALL PSTI2C(IGAP,3,STRG)
         CTEMP='['//STRL(1:3)//' '//STRG(1:3)//'] 0 h '
         CALL PSTBUF(14,CTEMP)

      ENDIF

C     redefine the postscript current position

C     the code below is equivalent to
C      call vbvect(0,xcp,ycp)
C     but can't do it because vbvect calls vdstlw which calls this routine

      CTEMP=COORD(1:11)//' m '
      CALL PSTBUF(14,CTEMP)

      VECTOR(4)=LINSTY

  999 RETURN
      END
      SUBROUTINE VDSTCS(YSIZE)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDSTCS           -Set Character Size.

C R.W.Simons       -05DEC80
C J. P. LONG       -03 DEC 87

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   All Devices that support only software characters.
C                   (LXY, HC1)

C ENTRY CONDITIONS -YSIZE = real Y dimension of the character box in NDC
C                   space.  Range 0.-1.  Default: device dependent,
C                   typically the smallest hardware size.

C CALLS            -

C EXIT CONDITIONS  -VECTOR(6) = real updated character box Y (YSIZE).
C                   VECTOR(7) = real updated character box X.

C NARRATIVE        -Set the character size for text primitives.  Size
C                   is given by YSIZE as the Y dimension of the
C                   character box.  The SVDI will assign the X dimension
C                   of the character box and X and Y character size
C                   within the box according to the font used.  Applies
C                   only to text primitives.
C                   All devices must support at least a single device
C                   dependent value that is the default.  If an
C                   unsupported value is specified, set to the largest
C                   supported character size that does not exceed the
C                   specified size.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      REAL YSIZE
      CHARACTER STR*4,CTEMP*10

C CURRENT ATTRIBUTE VALUES. (DEVICE-INDEPENDENT)
      REAL VECTOR(7)
      COMMON /VCATTR/ VECTOR
C     VECTOR(1)=FOREGROUND COLOR
C           (2)=BACKGROUND COLOR
C           (3)=INTENSITY
C           (4)=LINE STYLE
C           (5)=LINE WIDTH
C           (6)=CHARACTER BOX Y
C           (7)=CHARACTER BOX X
C SCALE FACTORS FOR NDC TO DC MAPPING. (LXY,HC1)
      REAL XSCALE,YSCALE
      COMMON /VCSCAL/ XSCALE,YSCALE

C CHECK FOR VALID YSIZE.
      IF(YSIZE.LT.0.0.OR.YSIZE.GT.1.0) THEN
         CALL VBERRH(401,5)
         GOTO 999
      END IF

C PROTECT INPUT PARAMETER FROM BEING CHANGED.
      YSIZE1=YSIZE

C DON'T ALLOW VALUES BELOW THE MINIMUM "HARDWARE" SIZE.
      IF(YSIZE1.LT.0.01) YSIZE1=0.01

C VALUES ESTABLISHED HERE ARE USED BY VBSIM IN SIMULATING CHARACTERS.
C ALWAYS USE A CHARACTER ASPECT RATIO OF 5/7.
      VECTOR(6)=YSIZE1
      VECTOR(7)=YSIZE1*5./7.

C convert the character size into device coords

      IYSIZE=NINT(XSCALE*YSIZE1)

C output the postscript command

      CALL PSTI2C(IYSIZE,4,STR)
C     iysize is in tenths of device units
      CTEMP='y '//STR(1:3)//' x '
      CALL PSTBUF(8,CTEMP)

  999 RETURN
      END
      SUBROUTINE VDSTLW(LINWTH)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDSTLW           -Set Line Width.

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77

C ENTRY CONDITIONS -LINWTH = real line width of line drawing output
C                   primitives.  Range 0.-1.  Default: device dependent.

C CALLS            -

C EXIT CONDITIONS  -VECTOR(5) = real updated line width (LINWTH).

C NARRATIVE        -Set the relative width of an output line.  Values
C                   are 0.-1. with 1. being .01 in NDC space.
C                   All devices must support at least a single device
C                   dependent value that is the default.  If an
C                   unsupported value is specified, set to the closest
C                   supported line width.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      REAL LINWTH,LW
      CHARACTER CTEMP*19,STR*5

      COMMON /VCVEC1/ IVECT
      INTEGER IVECT

C COMPUTER DEPENDENT COMMON VARIABLES AND CONSTANTS.
      include 'vcpstc.blk'
C CURRENT ATTRIBUTE VALUES. (DEVICE-INDEPENDENT)
      REAL VECTOR(7)
      COMMON /VCATTR/ VECTOR
C     VECTOR(1)=FOREGROUND COLOR
C           (2)=BACKGROUND COLOR
C           (3)=INTENSITY
C           (4)=LINE STYLE
C           (5)=LINE WIDTH
C           (6)=CHARACTER BOX Y
C           (7)=CHARACTER BOX X
C SCALE FACTORS FOR NDC TO DC MAPPING. (LXY,HC1)
      REAL XSCALE,YSCALE
      COMMON /VCSCAL/ XSCALE,YSCALE

C     pstmlw controls minimum line width
C     kpstbg controls background coloring
C            = 0,  not colored (white ground from paper)
C            = 1,  colored black
C     kpstci controls black-white interchange (colors 0 & 7 only)
C            = 0,  no interchange
C            = 1,  colors interchanged
      COMMON /VCPSTA/ PSTMLW, KPSTBG, KPSTCI

C CHECK FOR VALID LINWTH.
      IF(LINWTH.LT.0.0.OR.LINWTH.GT.1.) THEN
         CALL VBERRH(401,5)
         GOTO 999
      END IF

C     test user define minimum

      WIDTH=MAX(PSTMLW,LINWTH)

C CONVERT LINE-WIDTH TO NDC
      LW=WIDTH*.005

C CONVERT WIDTH TO DEVICE COORDINATES AND ADD A DIGIT; NEED IT TO HUNDREDTHS
      ILW=NINT(XSCALE*LW*10.)
C     A LINEWIDTH OF ZERO WORKS ONLY PART OF THE TIME
      IF(ILW.LT.10) ILW=10

C SET LINE WIDTH
      CALL PSTI2C(ILW,5,STR)
      IF(IVECT.NE.0) THEN
        CTEMP='s '//STR(1:3)//'.'//STR(4:5)//' w '
        CALL PSTBUF(11,CTEMP)
        IVECT=0
      ELSE
        CTEMP=STR(1:3)//'.'//STR(4:5)//' w '
        CALL PSTBUF(9,CTEMP)
      END IF

      VECTOR(5)=WIDTH

C     since linestyle uses the linewidth in setting the pattern, call it

      LINSTY=VECTOR(4)
      CALL VBSTLS(LINSTY)
  999 RETURN
      END
      SUBROUTINE VDIQES(ESCPCD,SUPPORT)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDIQES           -Inquire Escape.

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77

C ENTRY CONDITIONS -ESCPCD = integer escape function code.

C CALLS            -

C EXIT CONDITIONS  -SUPPRT = integer level of support for the escape
C                   function specified.  Range 0,1,2.

C NARRATIVE        -An integer value indicating 2=hardware supported,
C                   1=software supported, 0=unsupported is returned in
C                   SUPPORT for the escape function ESCPCD.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      INTEGER ESCPCD,SUPPORT
      IF (ESCPCD.EQ.2100) THEN
         SUPPORT=2
      ELSEIF ((ESCPCD.GE.2101).AND.(ESCPCD.LE.2110)) THEN
         SUPPORT=2
C ELSE THERE IS NO SUPPORT OF ANY OTHER ESCAPE CODES
      ELSE
         SUPPORT=0
      END IF
      RETURN
      END
      SUBROUTINE PSTBUF(NCHRS,OUT)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C PSTBUF           -Output PostScript data

C C. D. Brown      -DEC 1986 (Adapted from QMSBUF)

C ENVIRONMENT      -COMPUTER/DEVICE DEPENDENT, FORTRAN 77

C ENTRY CONDITIONS -NCHRS = integer number of characters in OUT.
C                          = 0 means flush the buffer.
C                   OUT = character string of input data
C                   KOUTFL = integer number of the graphics output file.

C CALLS            -

C EXIT CONDITIONS  -

C NARRATIVE        -The data in OUT is buffered for output to KOUTFL.
C                   The buffer is output when it is "full" or a buffer
C                   flush is requested by specifying NCHRS<=0.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      SAVE
      INTEGER NCHRS
      CHARACTER OUT*(*)
      character*132 lstout
      integer lstchr

C COMPUTER DEPENDENT COMMON VARIABLES AND CONSTANTS.
      include 'vcpstc.blk'

C **NOTE: BUFFER SIZE (IN BITS) MUST BE AN EXACT MULTIPLE OF 8 (8-BIT DATA
C  MUST END EXACTLY AT WORD BOUNDARY)
      INTEGER CHARLN,ICNT,REMAIN
      CHARACTER CBUF*130
C                                        CHARLN=BUFFER SIZE IN CHARS
      DATA ICNT/1/,CHARLN/130/,LSTCHR/-1/,LSTOUT/' '/

C ...Check that last output string does not match current output GDS
      if (lstchr .eq. nchrs) then
        if (lstout(:lstchr) .eq. out(:nchrs)) return
      end if
      lstchr = nchrs
      lstout(:nchrs) = out(:nchrs)

C COMPUTE REMAINING AVAILABLE CHARACTERS IN BUFFER
      REMAIN=CHARLN-ICNT+1

C CHECK FOR BUFFER FLUSH REQUEST OR NOT ENOUGH ROOM IN BUFFER.
      IF((NCHRS.LE.0).OR.(NCHRS.GT.REMAIN)) THEN
C                                        TEST IF THERE'S ANYTHING TO FLUSH.
         IF (ICNT.GT.1) THEN
C                                   PAD TO END OF RECORD AND OUTPUT THE BUFFER.
            IF (ICNT .LE. CHARLN) THEN
               CBUF(ICNT:CHARLN)=' '
               WRITE(KOUTFL,'(A)') CBUF(1:ICNT)
            ELSE
               WRITE(KOUTFL,'(A)') CBUF
            END IF
            ICNT=1
         ENDIF
      ENDIF

C ADD TO BUFFER
      IF (NCHRS.GT.0) THEN
         CBUF(ICNT:ICNT+NCHRS-1)=OUT(1:NCHRS)
         ICNT=ICNT+NCHRS
      ENDIF

      RETURN
      END
      SUBROUTINE PSTA2C(ASCI,CHARAC)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C PSTA2C           - CONVERT FROM ASCII TO CHARACTER

C P. Watterberg    - 19 Jan 1982

C ENVIRONMENT      - computer dependent, system dependent, fortran 77

C ENTRY CONDITIONS - ASCI is an integer representing an ascii character

C CALLS            -

C EXIT CONDITIONS  - CHARAC is the character represented by ASCI

C NARRATIVE        -

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      CHARACTER CHARAC*(*)
      INTEGER ASCI

      CHARAC = CHAR(ASCI)

      return
      end
      SUBROUTINE PSTI2C(INT,NDIGIT,ISTR)
C     C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C
C
C     PSTI2C           - convert positive integer to decimal character
C     string equivalent
C
C     ENVIRONMENT      - COMPUTER-INdependent
C
C     ENTRY CONDITIONS - int = positive integer to be converted
C     ndigit = number of digits to be produced in string
C     form (pad left with zeros)
C     istr = character string of at least ndigit characters
C
C     CALLS            -
C
C     EXIT CONDITIONS  - istr contains decimal-string equivalent of int
C     (ndigits left-justified in istr)
C
C     NARRATIVE        - This routine modified 10/89  S.L.Thompson
C
C     C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C
      INTEGER INT,NDIGIT
      CHARACTER ISTR*(*)
      CHARACTER*1 KA(10)
      DATA KA /'0','1','2','3','4','5','6','7','8','9'/
C
C     check input parameters
      INT1=MAX(INT,0)
      LENGTH=LEN(ISTR)
      NDIG1=MAX(1,MIN(LENGTH,NDIGIT))
      ISTR='00000000000000000000000000000000000000000'
      ND=LENGTH
      DO I=1,NDIG1
         J=INT1/10
         K=INT1-10*J
         ISTR(ND:ND)=KA(K+1)
         ND=ND-1
         INT1=J
      end do
      RETURN
      END
      SUBROUTINE PSTBBG

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C     Color background black for white paper device.
C     Should only be called from vdnwpg and viinit.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C     pstmlw controls minimum line width
C     kpstbg controls background coloring
C            = 0,  not colored (white ground from paper)
C            = 1,  colored black
C     kpstci controls black-white interchange (colors 0 & 7 only)
C            = 0,  no interchange
C            = 1,  colors interchanged
      COMMON /VCPSTA/ PSTMLW, KPSTBG, KPSTCI

C     mopoly controls polygon fill =0, on ;  =1, off
C     mocolr controls color =0, on ;  =1, off
      COMMON /VCPSTB/ MOPOLY, MOCOLR

C CURRENT ATTRIBUTE VALUES. (DEVICE-INDEPENDENT)
      REAL VECTOR(7)
      COMMON /VCATTR/ VECTOR
C     VECTOR(1)=FOREGROUND COLOR
C           (2)=BACKGROUND COLOR
C           (3)=INTENSITY
C           (4)=LINE STYLE
C           (5)=LINE WIDTH
C           (6)=CHARACTER BOX Y
C           (7)=CHARACTER BOX X

      COMMON /VCVEC1/ IVECT
      INTEGER IVECT
      COMMON /VCESCP/ PGFORM,PATNO,BORDER
      INTEGER PGFORM,PATNO,BORDER
      DIMENSION X(4),Y(4)
      PARAMETER (ONEN=0.99999)
      PARAMETER (ASP=0.75)

      IF(MOPOLY.EQ.0) THEN
        IF(PGFORM.EQ.0) THEN
          X(1)=0.
          X(2)=0.
          X(3)=ONEN
          X(4)=ONEN
          Y(1)=0.
          Y(2)=ASP
          Y(3)=ASP
          Y(4)=0.
        ELSE
          X(1)=0.
          X(2)=0.
          X(3)=ASP
          X(4)=ASP
          Y(1)=0.
          Y(2)=ONEN
          Y(3)=ONEN
          Y(4)=0.
        END IF
        KOLSAV=NINT(VECTOR(1))
        CALL VDSTFC(nint(vector(2)))
        CALL VIPOLY(X,Y,4)
        CALL VDSTFC(KOLSAV)
      END IF
      RETURN
      END
      SUBROUTINE PSTJOB
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C PSTJOB           - GET JOB ID AND ROUTING INFORMATION

C ENVIRONMENT      - COMPUTER-DEPENDENT FOR CTSS

C ENTRY CONDITIONS -

C CALLS            -

C EXIT CONDITIONS  - KJTIME - TIME HOLLERITH STRING
C                    KJDATE - DATE HOLLERITH STRING
C                    KUSRID - USER IDENTIFICATION
C                    KJROUT - ROUTING INFORMATION

C NARRATIVE        - THIS ROUTINE INQUIRES THE SYSTEM TO FIND THE ABOVE
C                    INFORMATION.  THE INFO IS PACKED INTO THE ARRAYS AS
C                    HOLLERITH (INTERNAL DISPLAY CODE) STRINGS.  A TERMI
C                    CHARACTER "\" IS APPENDED TO EACH STRING SO THE CAL
C                    ROUTINE CAN FIND THE END IF FOR SOME REASON THE LEN
C                    VARIABLES ARE NOT SUFFICIENT.

C  None of functions are used in pst driver

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C JOB ID INFORMATION. (HC1, DIC)
      include 'vcjob.blk'

C  FOR SECURITY MARKINGS, CTSS CODES NEED TO MAP TO THESE SILLY
C  OLD SCOPE SECURITY CODES

C      SCOPE 3 CODE

C         0    UNCL
C         1    UNDEFINED
C         2    UNDEFINED
C         3    PARD
C         4    C
C         5    CNSI
C         6    CFRD
C         7    CRD
C         8    S
C         9    SNSI
C        10    SFRD
C        11    SRD

C     GET CLASSIFICATION LEVEL
      KSECUR = 0

C     GET USER ID
      KUSRSZ = 8
      KUSRID(1)=0
      KUSRID(2)=0
      KUSRID(3)=0
      KUSRID(4)=0

C     GET JOB ID AND USERS NAME
      KJOBID(1) = 0
      KJOBID(2) = 0
      KJOBID(3) = 0
      KJOBID(4) = 0
      KIDSIZ = 24

C     GET BOX NUMBER
      KSZROU = 777

      KJROUT(1) = 0
      KJROUT(2) = 0
      KJROUT(3) = 0
      KJROUT(4) = 0

C GET MACHINE ID
      MACHIN(1) = 0
      MACHIN(2) = 0
      MACHIN(3) = 0
      MACLEN=1

C GET THE TIME AND DATE
      KJTIME(1)=0
      KJTIME(2)=0
      KJTIME(3)=0
      KJDATE(1)=0
      KJDATE(2)=0
      KJDATE(3)=0

      END
      SUBROUTINE PSTSEL(KARG)

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C     Select type of desired output. Four options are

C                                                      device number
C   1. black & white, batch, no poly fill                  799.1
C   2. black & white, interactive, no poly                 799.2
C   3. black & white, batch, poly fill                     799.3
C   4. black & white, interactive, poly fill               799.4
C   5. color, batch                                        799.5
C   6. color, interactive                                  799.6
C   7. color, batch, black-white interchange               799.7
C   8. color, interactive, black-white interchange         799.8
C   9. color, batch, black background                      799.9
C   10.color, interactive, black background                799.11

C     A second function of this routine is to set the minimum line
C     width. For most systems the minimum width line is too narrow.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      SAVE
      CHARACTER*(*) KARG

C     pstmlw controls minimum line width
C     kpstbg controls background coloring
C            = 0,  not colored (white ground from paper)
C            = 1,  colored black
C     kpstci controls black-white interchange (colors 0 & 7 only)
C            = 0,  no interchange
C            = 1,  colors interchanged
      COMMON /VCPSTA/ PSTMLW, KPSTBG, KPSTCI

C     mopoly controls polygon fill =0, on ;  =1, off
C     mocolr controls color =0, on ;  =1, off
      COMMON /VCPSTB/ MOPOLY, MOCOLR

      COMMON /DEVCAP/ DEV(33)
      common /blotans/ BLTANS
      character*2 BLTANS

      CHARACTER*2 ANS,ARG
      DATA IONCE /0/
      ARG=KARG
      IF(IONCE.EQ.0) THEN
        KPSTBG=0
        KPSTCI=0
        IONCE=1
        IF(ARG.EQ.' ' .and. bltans .eq. ' ') THEN
          WRITE(*,10)
   10     FORMAT(/,' This VDI PostScript driver has seven options.',/,
     &     '     1. black & white, no polygon fill',/,
     &     '     3. black & white, polygon fill',/,
     &     '     5. color,',/,
     &     '     7. color, black-white interchange',/,
     $     '     8. gray-scale, black-white interchange',/,
     &     '     9. color, black background',/,
     $     '    10. gray-scale, black background',/,
     &     ' Enter option number')
          READ(5,'(A)',ERR=30) ANS
        ELSE
          if (arg .ne. ' ') then
            ANS=ARG
          else if (bltans .ne. ' ') then
            ans = bltans
          end if
        END IF
        IF(ANS.EQ.'6') THEN
            DEV(4)=256.
            DEV(23)=799.6
            DEV(26)=1.
            DEV(27)=256.
            DEV(32)=1.
            MOPOLY=0
            MOCOLR=0
        ELSEIF(ANS.EQ.'5') THEN
            DEV(4)=256.
            DEV(23)=799.5
            DEV(26)=1.
            DEV(27)=256.
            DEV(32)=1.
            MOPOLY=0
            MOCOLR=0
        ELSEIF(ANS.EQ.'2') THEN
            DEV(4)=1.
            DEV(23)=799.2
            DEV(26)=1.
            DEV(27)=1.
            DEV(32)=0.
            MOPOLY=1
            MOCOLR=1
        ELSEIF(ANS.EQ.'1') THEN
            DEV(4)=1.
            DEV(23)=799.1
            DEV(26)=1.
            DEV(27)=1.
            DEV(32)=0.
            MOPOLY=1
            MOCOLR=1
        ELSEIF(ANS.EQ.'4') THEN
            DEV(4)=1.
            DEV(23)=799.4
            DEV(26)=1.
            DEV(27)=1.
            DEV(32)=0.
            MOPOLY=0
            MOCOLR=1
        ELSEIF(ANS.EQ.'3') THEN
            DEV(4)=1.
            DEV(23)=799.3
            DEV(26)=1.
            DEV(27)=1.
            DEV(32)=0.
            MOPOLY=0
            MOCOLR=1
        ELSEIF(ANS.EQ.'7') THEN
            DEV(4)=256.
            DEV(23)=799.7
            DEV(26)=1.
            DEV(27)=256.
            DEV(32)=1.
            MOPOLY=0
            MOCOLR=0
            KPSTCI=1
        ELSEIF(ANS.EQ.'8') THEN
            DEV(4)=256.
            DEV(23)=799.8
            DEV(26)=1.
            DEV(27)=256.
            DEV(32)=1.
            MOPOLY=0
            MOCOLR=0
            KPSTCI=1
        ELSEIF(ANS.EQ.'9') THEN
            DEV(4)=256.
            DEV(23)=799.9
            DEV(26)=1.
            DEV(27)=256.
            DEV(32)=1.
            MOPOLY=0
            MOCOLR=0
            KPSTBG=1
        ELSEIF(ANS.EQ.'10') THEN
            DEV(4)=256.
            DEV(23)=799.10
            DEV(26)=1.
            DEV(27)=256.
            DEV(32)=1.
            MOPOLY=0
            MOCOLR=0
            KPSTBG=1
        ELSE
            GO TO 30
        END IF
        GO TO 50
   30      WRITE(6,40)
   40      FORMAT(' Bad input - defaulting to 7')
            DEV(4)=256.
            DEV(23)=799.7
            DEV(26)=1.
            DEV(27)=256.
            DEV(32)=1.
            MOPOLY=0
            MOCOLR=0
            KPSTCI=1
           ARG=' '
   50   CONTINUE
*- INCLUDE PSTMLW
C       set minimum line width (range 0 to 1)
        PSTMLW=0.025
*-
      END IF
      RETURN
      END
      SUBROUTINE PSTINI

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C     Generate any system dependent records which must be at the first
C     of PostScript output file. For example, a SUN laser printer
C     requires the first record of the file to be %! for the file
C     recognized as a PostScript file. This routine writes these
C     initial records.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      SAVE
      DATA KWAY /0/
      IF(KWAY.EQ.0) THEN
        KWAY=1

C       generate first records in output file

*- INCLUDE PSTHEAD
C       the following is for a SUN UNIX system
C       record is a comment except for sun lpr
        CALL PSTBUF(14,'%!PS-Adobe-2.0')
C       clear line buffer
        CALL PSTBUF(0,' ')
*-
      END IF

      RETURN
      END
