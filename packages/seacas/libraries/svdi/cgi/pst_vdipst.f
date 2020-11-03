C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C VDMOVA VDMONI VDGNAM VBIQDV VBIQPK VDLINA VDTEXT VDPNTA VDPOLY VDIQCP VDSTOS
C WPSTMV WPSTMO WPSTGN WPSTIV WPSTQP WPSTLN WPSTTX WPSTPT WPSTPY WPSTCP WPSTOS

C VDIQOS VDSTFC VDSTBC VDSTIN VDSTLS VDSTLW VDSTCS VDAABU VDALOC VDABGL VDAKGL
C WPSTIO WPSTFC WPSTBC WPSTIN WPSTLS WPSTLW WPSTCS WPSTBU WPSTLO WPSTBL WPSTKL

C VDSTLA VDINIT VDFRAM VDTERM VDIQDC VDNWPG VDBELL VDWAIT VDBUFL VDSTCO VDIQCO
C WPSTLA WPSTNT WPSTFR WPSTTR WPSTDC WPSTPG WPSTBE WPSTWT WPSTFL WPSTCO WPSTIC

C VDESCP VDIQES VDIQND VIMOVA VILINA VIPNTA VITEXT VIINIT VITERM VINWPG CDRCOM
C WPSTES WPSTIE WPSTID WPSTIM WPSTIL WPSTIP WPSTIX WPSTII WPSTIT WPSTIG CDRCOM

C VCJOB  VCONOD VBERRH VDLOGE CDRWFS CDRRFS CDROFS CDROF3 CDRCFS CDROFF CDROAB
C  VCJOB WPSTON WPSTER WPSTLE WPSTWF WPSTRF WPSTOF WPSTO3 WPSTCF WPSTFF WPSTAB

C BGPBUF QMSBUF QMSBU1 DDCBUF H75BUF BTKBUF NMTBUF ONHBUF VBIMBF VBPKG  VBDEV
C WPSTBF WPSTQM WPSTBF WPSTBF WPSTBF WPSTBF WPSTBF WPSTOH WPSTIB WPSTPK WPSTDV

C VDIQRS VDSTMP VDSTRS VDSTRV VDBRGB VDFRGB VDPIXL VDPIXI VDRPIX VDRPXI VDRSCL
C WPSTQR WPSTMP WPSTRS WPSTRV WPSTBG WPSTFG WPSTPX WPSTPI WPSTRP WPSTRI WPSTRL

C VDIQCI VBSTMP VCNDCM VCATTR VIFRAM VCCRPS VCESCP DEVCAP VCSCAL VCDDIM VIPOLY
C WPSTCI WPST01 WPST02 WPST03 WPST04 WPST05 WPST06 WPST07 WPST08 WPST09 WPST10

C VCVECT VBLINA VBVECT VBSTLS PSTBUF
C WPST11 WPST12 WPST13 WPST14 WPST15

      SUBROUTINE WPST01( IMAP )
      GOTO (1,2,3,4,5),IMAP

      CALL WPSTMP('UNKNOWN')
      RETURN

    1 CALL WPSTMP('1-TO-1')
      RETURN

    2 CALL WPSTMP('REPLICATE')
      RETURN

    3 CALL WPSTMP('VIRTUAL')
      RETURN

    4 CALL WPSTMP('NODISTORT')
      RETURN

    5 CALL WPSTMP('FREEFORM')
      RETURN

      END
      SUBROUTINE WPSTRS(I1,I2)
c*************************************************************************
c      This routine is to satisfy entry points used with raster vdi stuff
c      but not with regular vdi stuff.  This is done so raster vdi programs
c      can link with regular vdi.
c*************************************************************************
      CHARACTER C1*(*)
      REAL RA1(1),RA2(1),RA3(1)
      INTEGER IA1(1),IA2(1)
      ENTRY WPSTRV(R1,R2,R3,R4)
      ENTRY WPSTMP(C1)
      ENTRY WPSTPX(I1,I2,RA1,RA2,RA3,I3)
      ENTRY WPSTRP(I1,I2,I3,RA1,RA2,RA3,IA1)
      ENTRY WPSTPI(I1,I2,IA1,I3)
      ENTRY WPSTRI(I1,I2,I3,IA1,IA2)
      ENTRY WPSTQR(I1,RA1)
      ENTRY WPSTRL
      ENTRY WPSTFG(R1,R2,R3)
      ENTRY WPSTBG(R1,R2,R3)
      ENTRY WPSTCI(R1,R2,R3,I1)
      RETURN
      END
      SUBROUTINE WPSTER(ERRNUM,ERRSEV)
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
      CALL WPSTLE(ERRNUM,ERRSEV)

C CHECK FOR FATAL ERROR.
      IF(ERRSEV.GT.12) STOP

      RETURN
      END
      SUBROUTINE WPSTGN(NAME)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDGNAM           -Name the graphics output file

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

      CHARACTER*(*) NAME
      INTEGER LENGTH,ISTART,IEND,I
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
      IF(ISTART .GT. 0)THEN
        LENGTH = IEND-ISTART+1
        CALL CDRGNM(NAME(ISTART:IEND),LENGTH)
      ENDIF
      RETURN
      END
      SUBROUTINE WPSTNT(ASPECT,JUSTIF)
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

C CALLS            -CDRJOB, VBERRH, VIINIT.

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

      INTEGER*4 MACHIN(3),MACLEN
      INTEGER*4 KIDSIZ,KJOBID(4),KUSRSZ,KUSRID(4),KSZROU
      INTEGER*4 KJROUT(4),KSECUR,KJTIME(4),KJDATE(4)
      COMMON / VCJOB/ KIDSIZ,KJOBID,KUSRSZ,KUSRID,KSZROU,
     1               KJROUT,KSECUR,KJTIME,KJDATE,MACHIN,MACLEN

C THIS IS JUST A DUMMY ROUTINE WHICH DOES NOTHING BUT CALL VIINIT.
C THIS ORGANIZATION FACILITATES ADDING SECURITY MARKINGS TO SVDI.
      CALL WPSTII(ASPECT,JUSTIF)

      RETURN
      END
      SUBROUTINE WPSTID(XNDC,YNDC)
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

      REAL XNDCMX,YNDCMX
      COMMON /WPST02/ XNDCMX,YNDCMX

C RETURN THE MAXIMUM VALID NDC VALUES.
      XNDC=XNDCMX
      YNDC=YNDCMX

      RETURN
      END
      SUBROUTINE WPSTIO(ATTARR)
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

      REAL VECTOR(7)
      COMMON /WPST03/ VECTOR

      INTEGER I

      DO 100 I=1,7
         ATTARR(I)=VECTOR(I)
  100 CONTINUE

      RETURN
      END
      SUBROUTINE WPSTLN(X,Y)
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
      CALL WPSTIL(X,Y)

      RETURN
      END
      SUBROUTINE WPSTLE(ERRNUM,ERRSEV)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDLOGE           -Log Error.

C R.W.Simons       -08APR81
C K.M.Erickson     -8OCT84 - add buffer flush

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   Device-independent.

C ENTRY CONDITIONS -ERRNUM = integer error number.
C                   ERRSEV = integer error severity.

C CALLS            -CDRTBK, VDBUFL

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

      INTEGER*4 KWRTFL,KRDFL,KOUTFL,KINFL,KWRDSZ,KBYTEL,KCPW,KBAUD,
     1KCOMTP
      COMMON /CDRCOM/ KWRTFL,KRDFL,KOUTFL,KINFL,KWRDSZ,KBYTEL,KCPW,
     1 KBAUD,KCOMTP

C flush buffer before we do a write
      CALL WPSTFL

C WRITE THE ERROR TO THE LISTING.
      WRITE(KWRTFL,10)ERRNUM,ERRSEV
   10 FORMAT(' SVDI ERROR NUMBER ',I5,'   SEVERITY CODE ',I5)

C TRACEBACK.
      CALL CDRTBK

      RETURN
      END
      SUBROUTINE WPSTMO(ISTATE)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDMONI           -Logs Usage Information..

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   Device-independent.

C ENTRY CONDITIONS -ISTATE = 0 - initialization
C                            1 - new page
C                            2 - terminate

C CALLS            -CDRMON

C EXIT CONDITIONS  -

c NARRATIVE        -For ISTATE=0, job information is initialized, and
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

      CHARACTER *8 C1,C2,MDEV,MPKG
      INTEGER ISTATE,MPAGES
      SAVE MDEV,MPKG,MPAGES
      DATA MPKG /'        '/
      DATA MDEV /'        '/
      DATA MPAGES /0/

      IF(ISTATE.EQ.0) THEN
          CALL CDRELA(0)
      ELSEIF (ISTATE.EQ.1) THEN
          MPAGES=MPAGES+1
      ELSE
          CALL CDRELA(1)
          CALL CDRMON(MDEV,MPKG,MPAGES)
      ENDIF
      RETURN
C Usage Monitoring Information

      ENTRY WPSTPK (C1)
      MPKG = C1
      RETURN
      ENTRY WPSTDV (C2)
      MDEV = C2
      RETURN
      ENTRY WPSTQP(C1)
      C1 = MPKG
      RETURN
      ENTRY WPSTIV(C2)
      C2 = MDEV
      RETURN
      END
      SUBROUTINE WPSTMV(X,Y)
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
      CALL WPSTIM(X,Y)

      RETURN
      END
      SUBROUTINE WPSTPG
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
      CALL WPSTIG

      RETURN
      END
      SUBROUTINE WPSTPT(X,Y)
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
      CALL WPSTIP(X,Y)

      RETURN
      END
      SUBROUTINE WPSTPY(XARRAY,YARRAY,NPTS)
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
      CALL WPST10(XARRAY,YARRAY,NPTS)

      RETURN
      END
      SUBROUTINE WPSTOS(ATTARR)
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

      REAL VECTOR(7)
      COMMON /WPST03/ VECTOR

C CALL EACH OF THE INDIVIDUAL ATTRIBUTE SETTING ROUTINES.
C CHECK FOR VALIDITY OF INPUT VALUES WILL BE DONE IN EACH INDIVIDUAL
C ROUTINE.
      CALL WPSTFC(INT(ATTARR(1)))
      CALL WPSTBC(INT(ATTARR(2)))
      CALL WPSTIN(ATTARR(3))
      CALL WPSTLS(INT(ATTARR(4)))
      CALL WPSTLW(ATTARR(5))
      CALL WPSTCS(ATTARR(6))

      RETURN
      END
      SUBROUTINE WPSTTR
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
      CALL WPSTIT

      RETURN
      END
      SUBROUTINE WPSTTX(LENGTH,CHARS)
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
C                   output; mapping lower; case to upper case letters is
C                   considered reasonable.  Attributes foreground color,
C                   background color, intensity, and  character size
C                   apply.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      INTEGER LENGTH,CHARS(136)

C THIS IS JUST A DUMMY ROUTINE WHICH DOES NOTHING BUT CALL VITEXT.
C THIS ORGANIZATION FACILITATES ADDING SECURITY NARKINGS TO SVDI.
      CALL WPSTIX(LENGTH,CHARS)

      RETURN
      END
      SUBROUTINE WPSTFR(ITYPE)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDFRAM           - Draw header or trailer frame

C P. Watterberg    - 27 Aug 81

C ENVIRONMENT      -Computer-independent, system-independent, FORTRAN 77

C ENTRY CONDITIONS - ITYPE = 0   for header frame
C                          = 1   for trailer frame

C CALLS            - VIFRAM

C EXIT CONDITIONS  -

C NARRATIVE        - Calls vifram to get time and date from the
c                    system via the computer-dependent routine CDRTOD(entry
c                    point in CDRJOB) and writes it on an identification frame.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      INTEGER ITYPE

      CALL WPST04(ITYPE)
      RETURN
      END
      SUBROUTINE WPST04(ITYPE)
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
      SUBROUTINE WPSTBU(BTNNUM)
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

      INTEGER*4 KWRTFL,KRDFL,KOUTFL,KINFL,KWRDSZ,KBYTEL,KCPW,KBAUD,
     1KCOMTP
      COMMON /CDRCOM/ KWRTFL,KRDFL,KOUTFL,KINFL,KWRDSZ,KBYTEL,KCPW,
     1 KBAUD,KCOMTP

C BATCH DEVICES DON'T NEED TO FLUSH BUFFERS.

C READ A RECORD FROM COMPUTER DEPENDENT FILE KINFL IN I5 FORMAT.
      READ(KINFL,10) BTNNUM
   10 FORMAT(I5)

C CHECK FOR VALID BTNNUM.
C RANGE FOR BATCH DEVICES IS 1-99999.  IF OUT OF RANGE, MAP IT BACK IN:
C MAPPING (-1)-(-9999) TO 1-9999 AND MAPPING 0 TO 10000.
      IF(BTNNUM.LT.0) BTNNUM=-BTNNUM
      IF(BTNNUM.EQ.0) BTNNUM=10000

      RETURN
      END
      SUBROUTINE WPSTBL(BTNNUM,X,Y)
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

      INTEGER*4 KWRTFL,KRDFL,KOUTFL,KINFL,KWRDSZ,KBYTEL,KCPW,KBAUD,
     1KCOMTP
      COMMON /CDRCOM/ KWRTFL,KRDFL,KOUTFL,KINFL,KWRDSZ,KBYTEL,KCPW,
     1 KBAUD,KCOMTP

C BATCH DEVICES DON'T NEED TO FLUSH BUFFERS.

C READ A RECORD FROM COMPUTER DEPENDENT FILE KINFL IN I5,2F10.7 FORMAT.
      READ(KINFL,10) BTNNUM,X,Y
   10 FORMAT(I5,2F10.7)

C CHECK FOR VALID BTNNUM.
C RANGE FOR BATCH DEVICES IS 1-99999.  IF OUT OF RANGE, MAP IT BACK IN:
C MAPPING (-1)-(-9999) TO 1-9999 AND MAPPING 0 TO 10000.
      IF(BTNNUM.LT.0) BTNNUM=-BTNNUM
      IF(BTNNUM.EQ.0) BTNNUM=10000

      RETURN
      END
      SUBROUTINE WPSTKL(CHAR,X,Y)
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

      INTEGER IN,CHR

      INTEGER*4 KWRTFL,KRDFL,KOUTFL,KINFL,KWRDSZ,KBYTEL,KCPW,KBAUD,
     1KCOMTP
      COMMON /CDRCOM/ KWRTFL,KRDFL,KOUTFL,KINFL,KWRDSZ,KBYTEL,KCPW,
     1 KBAUD,KCOMTP

C BATCH DEVICES DON'T NEED TO FLUSH BUFFERS.

C READ A CHARACTER AND AN X,Y POSITION FROM COMPUTER DEPENDENT FILE
C KINFL WITH FORMAT A1,2F10.7.
      READ(KINFL,10) CHR,X,Y
   10 FORMAT(A1,2F10.7)

C CONVERT CHARACTER TO INTEGER ASCII AND CHECK FOR VALID RANGE.
      CALL CDR1CH(1,CHR,IN)
      CALL CDRCVT(IN,CHAR)
      IF(CHAR.LT.32.OR.CHAR.GT.126) CHAR=32

      RETURN
      END
      SUBROUTINE WPSTLO(X,Y)
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

      INTEGER*4 KWRTFL,KRDFL,KOUTFL,KINFL,KWRDSZ,KBYTEL,KCPW,KBAUD,
     1KCOMTP
      COMMON /CDRCOM/ KWRTFL,KRDFL,KOUTFL,KINFL,KWRDSZ,KBYTEL,KCPW,
     1 KBAUD,KCOMTP

C BATCH DEVICES DON'T NEED TO FLUSH BUFFERS.

C READ AN X,Y POSITION FROM COMPUTER DEPENDENT FILE
C KINFL WITH FORMAT 2F10.7.
      READ(KINFL,10) X,Y
   10 FORMAT(2F10.7)

      RETURN
      END
      SUBROUTINE WPSTBE
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
      SUBROUTINE WPSTFL
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
      SUBROUTINE WPSTLA(LOCX,LOCY)
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
      SUBROUTINE WPSTWT
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
      SUBROUTINE WPSTIC(NUM,INDEX,CLRARY,CLRMOD)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDIQCO           -Inquire Color Table.

C R.W.Simons       -08APR81
C H. S. LAUSON      29MAY86 - changed for current HLS interpretation

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   All Black and White Devices. (LXY, HC1, ALP)

C ENTRY CONDITIONS -NUM = integer number of color indexes to inquire.
C                   Range 1-256.
C                   INDEX = integer array of indexes to inquire.  Range
C                   0-255.
C                   CLRMOD = integer color model to be used.  Range 0,1.

C CALLS            -VBERRH

C EXIT CONDITIONS  -CLRARY = real array of 3 by NUM elements returning
C                   the values of the components of the indexes inquired.
C                   Range for RGB: red 0.0-1.0
C                                  green 0.0-1.0
C                                  blue 0.0-1.0
C                   Range for HLS: hue 0.0-360.0
C                                  lightness 0.0-1.0
C                                  saturation 0.0-1.0

C NARRATIVE        -Inquire one or more color table entries.  NUM and
C                   INDEX specify how many and which indexes are being
C                   inquired.  CLRMOD specifies which color model
C                   (0=RGB, 1=HLS) should be used in constructing values
C                   to return in CLRARY.  A device which does not
C                   support a color table index specified will
C                   return -1.0 in the first element of the CLRARY value
C                   for that index.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      INTEGER NUM,INDEX(NUM),CLRMOD
      REAL CLRARY(3,NUM)

C CHECK FOR VALID NUM.
      IF(NUM.LT.1.OR.NUM.GT.256) THEN
         CALL WPSTER(723,5)
         GOTO 999
      END IF

C CHECK FOR VALID CLRMOD.
      IF(CLRMOD.NE.0.AND.CLRMOD.NE.1) THEN
         CALL WPSTER(725,5)
         GOTO 999
      END IF

C CHECK FOR VALID INDEXES.
      DO 100 I=1,NUM
         INDEXN=INDEX(I)
         IF(INDEXN.LT.0.OR.INDEXN.GT.255) THEN
            CALL WPSTER(724,5)
            GOTO 100
         END IF
         IF(INDEXN.EQ.0) THEN
C INDEX 0 IS ALWAYS BLACK
            IF (CLRMOD.EQ.0) THEN
                CLRARY(1,I)=0.0
              ELSE
                CLRARY(1,I)=120.
            ENDIF
            CLRARY(2,I)=0.0
            CLRARY(3,I)=0.0
         ELSE IF(INDEXN.EQ.7) THEN
C INDEX 7 IS ALWAYS WHITE
            IF (CLRMOD.EQ.0) THEN
                CLRARY(1,I)=1.0
                CLRARY(3,I)=1.0
              ELSE
                CLRARY(1,I)=120.
                CLRARY(3,I)=0.0
            ENDIF
            CLRARY(2,I)=1.0
         ELSE
C NO OTHER INDEXES ARE SUPPORTED.
            CLRARY(1,I)=-1.0
         END IF
  100 CONTINUE

  999 RETURN
      END
      SUBROUTINE WPSTCP(X,Y)
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

      REAL XCP,YCP
      COMMON /WPST05/ XCP,YCP

C ASSIGN THE CP TO X,Y.
      X=XCP
      Y=YCP

      RETURN
      END
      SUBROUTINE WPSTBC(COLOR)
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

      REAL VECTOR(7)
      COMMON /WPST03/ VECTOR

C CHECK FOR VALID COLOR.
      IF(COLOR.LT.0.OR.COLOR.GT.255) THEN
         CALL WPSTER(724,5)
         GOTO 999
      END IF

C ONLY THE SINGLE BACKGROUND COLOR 7 (WHITE) IS SUPPORTED,
C SO NO ACTION IS NECESSARY.

  999 RETURN
      END
      SUBROUTINE WPSTCO(NUM,INDEX,CLRARY,CLRMOD)
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

C CHECK FOR VALID NUM.
      IF(NUM.LT.1.OR.NUM.GT.256) THEN
         CALL WPSTER(723,5)
         GOTO 999
      END IF

C CHECK FOR VALID CLRMOD.
      IF(CLRMOD.NE.0.AND.CLRMOD.NE.1) THEN
         CALL WPSTER(725,5)
         GOTO 999
      END IF

C CHECK FOR VALID INDEXES.
      DO 100 I=1,NUM
         INDEXN=INDEX(I)
         IF(INDEXN.LT.0.OR.INDEXN.GT.255) THEN
            CALL WPSTER(724,5)
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
               CALL WPSTER(727,5)
               GOTO 100
            END IF

C ONLY TWO INDEXES ARE SUPPORTED:
C    0 WHICH IS THE DEFAULT FOREGROUND COLOR AND MUST ALWAYS REMAIN
C      BLACK (0.0, 0.0, 0.0)
C    7 WHICH IS THE DEFAULT BACKGROUND COLOR AND MUST ALWAYS REMAIN
C      WHITE (1.0, 1.0, 1.0)
C THEREFORE, NO ACTION IS NECESSARY.
         ELSE
            IF(CLRAR1.LT.0..OR.CLRAR1.GT.360.
     X      .OR.CLRAR2.LT.0..OR.CLRAR2.GT.1.
     X      .OR.CLRAR3.LT.0..OR.CLRAR3.GT.1.) THEN
               CALL WPSTER(727,5)
               GOTO 100
            END IF

C ONLY TWO INDEXES ARE SUPPORTED:
C    0 WHICH IS THE DEFAULT FOREGROUND COLOR AND MUST ALWAYS REMAIN
C      BLACK (0.0, 0.0, 0.0)
C    7 WHICH IS THE DEFAULT BACKGROUND COLOR AND MUST ALWAYS REMAIN
C      WHITE (1.0, 1.0, 1.0)
C THEREFORE, NO ACTION IS NECESSARY.
         END IF
  100 CONTINUE

  999 RETURN
      END
      SUBROUTINE WPSTFC(COLOR)
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

      REAL VECTOR(7)
      COMMON /WPST03/ VECTOR

C CHECK FOR VALID COLOR.
      IF(COLOR.LT.0.OR.COLOR.GT.255) THEN
         CALL WPSTER(724,5)
         GOTO 999
      END IF

C ONLY THE SINGLE FOREGROUND COLOR 0 (BLACK) IS SUPPORTED,
C SO NO ACTION IS NECESSARY.

  999 RETURN
      END
      SUBROUTINE WPSTIN(INTEN)
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

      REAL VECTOR(7)
      COMMON /WPST03/ VECTOR

C CHECK FOR VALID INTEN.
      IF(INTEN.LT.0.0.OR.INTEN.GT.1.0) THEN
         CALL WPSTER(401,5)
         GOTO 999
      END IF

C ONLY THE SINGLE INTENSITY 1.0 (MAXIMUM) IS SUPPORTED,
C SO NO ACTION IS NECESSARY.

  999 RETURN
      END
      SUBROUTINE WPSTDC(INDEX,VALUE)
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
C           input at any time; this input; is then saved by the system
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
c       12.  HPP - Hewlett Packard Printer/Plotter 2671G
C       12.1 H75 - HP 7580
C       12.2 H72 - HP 7221C OR T
C       12.3 H74 - HP 7475A
C       14.  RET - Retrographics
C       15.  AP5 - Aps 5 Phototypesetter
c       16.  JP7 - JUPITER 7
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

      INTEGER INDEX
      REAL VALUE
C ESCAPE FLAGS
C PATNO AND BORDER USED BY VIPOLY FOR FILL PATTERN AND BORDER ON/OFF;
C DEFAULT; COMPLETE FILL WITH BORDER. PLC.
      COMMON /WPST06/ PGFORM,PATNO,BORDER
      INTEGER PGFORM,PATNO,BORDER

C INITIALIZE THE DEVICE CAPABILITIES VECTOR.
      COMMON /WPST07/DEV
      REAL DEV(33)
      DEV(1)  = 0.
      DEV(2)  = 1.
      DEV(3)  = 1.
      DEV(4)  = 1.
      DEV(5)  = 15.
      DEV(6)  = 2.
      DEV(7)  = 0.
      DEV(8)  = 0.
      DEV(9)  = 0.
      DEV(10) = 0.
      DEV(11) = 0.
      DEV(12) = 0.
      DEV(13) = 0.
      DEV(14) = 0.
      DEV(15) = 7230.
      DEV(16) = 5040.
      DEV(17) = 254.
      DEV(18) = 178.
      DEV(19) = 4.
      DEV(20) = 10.
      DEV(21) = 84.
      DEV(22) = 0.
      DEV(23) = 99.
      DEV(24) = 3.
      DEV(25) = 99999.
      DEV(26) = 0.
      DEV(27) = 1.
      DEV(28) = 0.
      DEV(29) = 0.
      DEV(30) = 5000.
      DEV(31) = 750.
      DEV(32) = 0.
      DEV(33) = 1.

C CHECK FOR VALID INDEX.
      IF(INDEX.LT.1.OR.INDEX.GT.33) THEN
         CALL WPSTER(726,5)
         GOTO 999
      END IF

C RETURN INDEXED VALUE.
      VALUE=DEV(INDEX)

  999 RETURN
      END
      SUBROUTINE WPSTII(ASPECT,JUSTIF)
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

      REAL ASPECT
      INTEGER JUSTIF

      REAL XNDCMX,YNDCMX
      COMMON /WPST02/ XNDCMX,YNDCMX
      REAL VECTOR(7)
      COMMON /WPST03/ VECTOR
      REAL XCP,YCP
      COMMON /WPST05/ XCP,YCP
      INTEGER*4 KWRTFL,KRDFL,KOUTFL,KINFL,KWRDSZ,KBYTEL,KCPW,KBAUD,
     1KCOMTP
      COMMON /CDRCOM/ KWRTFL,KRDFL,KOUTFL,KINFL,KWRDSZ,KBYTEL,KCPW,
     1 KBAUD,KCOMTP
      REAL XSCALE,YSCALE
      COMMON /WPST08/ XSCALE,YSCALE
      INTEGER XPAD,YPAD,XDEVIC,YDEVIC
      COMMON /WPST09/ XPAD,YPAD,XDEVIC,YDEVIC
      COMMON /WPST11/ IVECT,COORD
      CHARACTER COORD*20
      INTEGER IVECT
      INTEGER*4 MACHIN(3),MACLEN
      INTEGER*4 KIDSIZ,KJOBID(4),KUSRSZ,KUSRID(4),KSZROU
      INTEGER*4 KJROUT(4),KSECUR,KJTIME(4),KJDATE(4)
      COMMON / VCJOB/ KIDSIZ,KJOBID,KUSRSZ,KUSRID,KSZROU,
     1               KJROUT,KSECUR,KJTIME,KJDATE,MACHIN,MACLEN

      COMMON /WPST07/ DEV
      REAL DEV(33)
C ESCAPE FLAGS
C PATNO AND BORDER USED BY VIPOLY FOR FILL PATTERN AND BORDER ON/OFF;
C DEFAULT; COMPLETE FILL WITH BORDER. PLC.
      COMMON /WPST06/ PGFORM,PATNO,BORDER
      INTEGER PGFORM,PATNO,BORDER
      CHARACTER XCOORD*5,YCOORD*5

C SET DEFAULT ATTRIBUTE VALUES.  ALL ARE DEVICE DEPENDENT EXCEPT
C VECTOR(4)=0.0.
C     VECTOR(1)=FOREGROUND COLOR - BLACK
C           (2)=BACKGROUND COLOR - WHITE
C           (3)=INTENSITY        - MAXIMUM
C           (4)=LINE STYLE       - SOLID
C           (5)=LINE WIDTH       - ABOUT 1/72 INCHES
C           (6)=CHARACTER BOX Y  - ABOUT 1/10 INCHES
C           (7)=CHARACTER BOX X  - 5/7 OF BOX-Y

      IVECT=0
      VECTOR(1) = 0.
      VECTOR(2) = 7.
      VECTOR(3) = 1.
      VECTOR(4) = 0.
      VECTOR(5) = 0.06255
      VECTOR(6) = 0.01
      VECTOR(7) = 0.
      PATNO  = 20
      BORDER = 1
      PGFORM = 0
      XCP    = 0.
      YCP    = 0.
C PROTECT INPUT PARAMETERS FROM BEING CHANGED.
      ASPEC1=ASPECT
      JUSTI1=JUSTIF

C CHECK FOR VALID ASPECT.  IF(ASPECT.LT.0.0) THEN CALL VBERRH(721,5),
C AND USE DEFAULT ASPECT.
      IF(ASPECT.LT.0.0) THEN
         CALL WPSTER(721,5)
         ASPEC1=0.0
      END IF

C CHECK FOR VALID JUSTIF.  IF(JUSTIF.LT.0 .OR. JUSTIF.GT.9) THEN
C CALL VBERRH(720,5), AND USE DEFAULT JUSTIF.
      IF(JUSTIF.LT.0.OR.JUSTIF.GT.9) THEN
         CALL WPSTER(720,5)
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
      XINCH=10.0
      YINCH=7.5
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
         XNDCMX=AMIN1(1.,ASPEC1)
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

C  SET UP MONITORING INFORMATION
      CALL WPSTDV('C PST   ')
      CALL WPSTMO(0)

C OPEN OUTPUT FILE
      CALL WPSTOF(KOUTFL)

C INITIALIZE THE QMS PS JET
      CALL WPST15(2,'%!')
      CALL WPST15(1,CHAR(13))
      CALL WPST15(1,CHAR(10))
      CALL WPST15(27,'/y {/Courier findfont} def ')
      CALL WPST15(27,'/x {scalefont setfont} def ')
      CALL WPST15(45,'initgraphics /m {moveto} def /l {lineto} def ')
      CALL WPST15
     *     (50,'/c {closepath} def /v {save} def /r {restore} def ')
      CALL WPST15
     *    (54,'/f {eofill} def /s {stroke} def /w {setlinewidth} def ')
      CALL WPST15(31,'/h {setdash} def /t {show} def ')
      CALL WPST15(33,'/d {gsave} def /u {grestore} def ')
      CALL WPST15(14,'1 setlinejoin ')
C                                       SET PAGE FORMAT (LANDSCAPE/PORTRAIT)
       IF (PGFORM.EQ.0) THEN
          CALL WPST15(4,'/o {')
          CALL WPST15(10,'90 rotate ')
          CALL CDRI2C(0,4,XCOORD)
          CALL CDRI2C(INT(YDEVIC+YDEVIC*3./32.),4,YCOORD)
          COORD = ' '//XCOORD(1:3)//'.'//XCOORD(4:4)//' -'//
     1    YCOORD(1:3)//'.'//YCOORD(4:4)
          CALL WPST15( 13,COORD)
          CALL WPST15(11,' translate ')
          CALL WPST15(6,'} def ')
          YPAD = -YPAD
       ELSE
          CALL WPST15(17,'/o {newpath} def ')
       ENDIF
c        define the postscript current position
      CALL WPST15(35,'/p {showpage} def 1 setlinecap v o ')
      CALL WPST13(0,XCP,YCP)

C INIT LINE WIDTH,CHARACTER SIZE
      CALL WPSTLW(VECTOR(5))
      CALL WPSTCS(VECTOR(6))
      RETURN
      END
      SUBROUTINE WPSTIT
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

      INTEGER*4 KWRTFL,KRDFL,KOUTFL,KINFL,KWRDSZ,KBYTEL,KCPW,KBAUD,
     1KCOMTP
      COMMON /CDRCOM/ KWRTFL,KRDFL,KOUTFL,KINFL,KWRDSZ,KBYTEL,KCPW,
     1 KBAUD,KCOMTP

c  put out the last page and restore postscript environment so
c                            nothing is left on the stack
      CALL WPSTIG
      CALL WPST15(2,'r ')
C FLUSH BUFFER
      CALL WPST15(0,' ')
C CLOSE OUTPUT FILE
      CALL WPSTCF(KOUTFL,1)
      CALL WPSTMO(2)

      RETURN
      END
      SUBROUTINE WPSTIM(X,Y)
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
      CALL WPST13(0,X,Y)

      RETURN
      END
      SUBROUTINE WPSTIP(X,Y)
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

      CALL WPSTIM(X,Y)
      CALL WPSTIL(X,Y)

      RETURN
      END
      SUBROUTINE WPST10(XARRAY,YARRAY,NPTS)
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

c ESCAPE FLAGS
C PATNO AND BORDER USED BY VIPOLY FOR PATTERN FILL AND BORDER ON/OFF. DEFAULT
C COMPLETE FILL AND BORDER ON
      REAL VECTOR(7)
      COMMON /WPST03/ VECTOR
      COMMON /WPST06/ PGFORM,PATNO,BORDER
      INTEGER PGFORM,PATNO,BORDER
      COMMON /WPST11/ IVECT,COORD
      CHARACTER COORD*20
      INTEGER IVECT

C CHECK FOR VALID N
      IF (NPTS.LT.1 .OR. NPTS.GT.1490) THEN
         CALL WPSTER(802,5)
         GO TO 999
      END IF

C IF A SET OF VECTORS WAS IN PROCESS, ISSUE STROKE COMMAND TO DRAW THEM
C Start a new path.

      CALL WPST15(2,'s ')
      IVECT=0

      CALL WPST15(2,'r ')
      CALL WPST15(0,' ')
      CALL WPST15(4,'v o ')
      CALL WPSTLW(VECTOR(5))
      CALL WPSTCS(VECTOR(6))

C DRAW POLYGON VECTORS

C MOVE TO FIRST POINT
      CALL WPSTMV(XARRAY(1),YARRAY(1))

C CALL VDLINA TO DRAW POINTS FROM 1ST POINT TO NTH POINT
      DO 100 I=2,NPTS
        CALL WPSTLN(XARRAY(I),YARRAY(I))
100   CONTINUE

C THEN DRAW A LINE TO THE FIRST POINT TO CLOSE THE POLYGON
      CALL WPSTLN(XARRAY(1),YARRAY(1))

C CLOSE THE POLYGON, GRAPHICS SAVE, FILL IT, GRAPHICS RESTORE, STROKE
C    TO PROVIDE THE SAME FILLED AREA AS IF IT WERE FILLED WITH VECTORS
C    THEN RESTORE AND SAVE POSTSCRIPT ENVIRONMENT TO AVOID INPUT BUFFER OVERFLOW
      CALL WPST15(13,' c d f u s r ')
      CALL WPST15(0,' ')
      CALL WPST15(4,'v o ')
      CALL WPSTLW(VECTOR(5))
      CALL WPSTCS(VECTOR(6))

C INIT THE CURRENT POSITION WITHIN POSTSCRIPT
      CALL WPSTMV(XARRAY(NPTS),YARRAY(NPTS))
      IVECT=0

999   RETURN
      END
      SUBROUTINE WPSTIG
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

      COMMON /WPST11/ IVECT,COORD
      CHARACTER COORD*20
      INTEGER IVECT
      REAL VECTOR(7)
      COMMON /WPST03/ VECTOR

c  stroke the path in case there are any vector and show text
         CALL WPST15(2,'s ')
         IVECT=0

c   showpage and restore postscript environment to avoid buffer overflow
c            flush buffer because save and restore won't work back-to-back

      CALL WPST15(4,'p r ')
      CALL WPST15(0,' ')
      CALL WPST15(4,'v o ')
      CALL WPSTLW(VECTOR(5))
      CALL WPSTCS(VECTOR(6))
      CALL WPSTMO(1)

      RETURN
      END
      SUBROUTINE WPSTES(ESCPCD,N,ARGS)
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
      INTEGER*4 KWRTFL,KRDFL,KOUTFL,KINFL,KWRDSZ,KBYTEL,KCPW,KBAUD,
     1KCOMTP
      COMMON /CDRCOM/ KWRTFL,KRDFL,KOUTFL,KINFL,KWRDSZ,KBYTEL,KCPW,
     1 KBAUD,KCOMTP
      REAL XSCALE,YSCALE
      COMMON /WPST08/ XSCALE,YSCALE
      INTEGER XPAD,YPAD,XDEVIC,YDEVIC
      COMMON /WPST09/ XPAD,YPAD,XDEVIC,YDEVIC

C USED BY VIPOLY FOR PATTERN FILL AND BORDER ON/OFF. DEFAULT COMPLETE FILL
C AND BORDER ON. PLC.
      COMMON/WPST06/PGFORM,PATNO,BORDER
      INTEGER PGFORM,PATNO,BORDER

C CHECK FOR VALID N.
      IF(N.LT.0) THEN
         CALL WPSTER(802,5)
         GOTO 999
      END IF

C 2100 - PAGE FORMAT (0=LANDSCAPE,1=PORTRAIT)
      IF (ESCPCD.EQ.2100) THEN
         IF (ARGS(1).EQ.0) THEN
            PGFORM=0
         ELSE
            PGFORM=1
         ENDIF
      ENDIF

 999  RETURN
      END
      SUBROUTINE WPSTIL (X,Y)
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

C draw
      ENTRY WPST12(X,Y)
      CALL WPST13(1,X,Y)

      RETURN
      END

      SUBROUTINE WPST13(IPEN,X,Y)
c****************************************************
c vbvect - do move or draw to x,y (depending on ipen)

c     ipen = 0 for move, 1 for draw
c     x,y = NDC coordinates to be moved/drawn to

c******************************************************

      REAL X,Y,XOFF,YOFF
      CHARACTER CTEMP*20,XCOORD*5,YCOORD*5
      REAL XSCALE,YSCALE
      COMMON /WPST08/ XSCALE,YSCALE
      INTEGER XPAD,YPAD,XDEVIC,YDEVIC
      COMMON /WPST09/ XPAD,YPAD,XDEVIC,YDEVIC
      REAL XCP,YCP
      COMMON /WPST05/ XCP,YCP
      REAL VECTOR(7)
      COMMON /WPST03/ VECTOR
      COMMON /WPST11/ IVECT,COORD
      CHARACTER COORD*20
      INTEGER IVECT

c compute new point in dev. coord.
c     convert to floating offsets
      XOFF=XPAD
      YOFF=YPAD

      IXDC=X*XSCALE+XOFF
      IYDC=Y*YSCALE+YOFF

c        write(xcoord,'(i5)')ixdc
c        write(ycoord,'(i5)')iydc
c                                ...include both x,y
      CALL CDRI2C(IXDC,4,XCOORD)
      CALL CDRI2C(IYDC,4,YCOORD)
      COORD = XCOORD(1:3)//'.'//XCOORD(4:4)//' '//
     1 YCOORD(1:3)//'.'//YCOORD(4:4)

c pack up move/draw command, send it down
      IF (IPEN.EQ.0) THEN
         CTEMP= COORD(1:11) // ' m '
      ELSE
         CTEMP= COORD(1:11) // ' l '
      ENDIF
      CALL WPST15(14,CTEMP)
c                          ...count the coordinate pair
      IVECT=IVECT+1

c  stroke the path if we are approaching the 1500-coord pair limit
c                also restore and save postscript environment to avoid
c                input buffer overflow (must have a c/r between restore
c                and save)
      IF(IVECT.GT.1400) THEN
         CALL WPST15(4,'s r ')
         CALL WPST15(0,' ')
         CALL WPST15(4,'v o ')
         CALL WPSTLW(VECTOR(5))
c        ...reset the vector count - vdstls (called by vdstlw)
c                 reinitted the current posn
         IVECT=1
      ENDIF

C UPDATE CURRENT POSITION
      XCP=X
      YCP=Y

      RETURN
      END
      SUBROUTINE WPSTIX(LENGT1,CHARS)
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
C                    output; mapping lower; case to upper case letters is
C                    considered reasonable.  Attributes foreground color,
C                    background color, intensity, and character size
C                    apply.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      INTEGER LENGT1, CHARS(136), LENGTH

      INTEGER XPAD,YPAD,XDEVIC,YDEVIC
      COMMON /WPST09/ XPAD,YPAD,XDEVIC,YDEVIC
      CHARACTER CTEMP*150,STR*3
      REAL XCP,YCP
      COMMON /WPST05/ XCP,YCP
      REAL VECTOR(7)
      COMMON /WPST03/ VECTOR
      REAL XSCALE,YSCALE
      COMMON /WPST08/ XSCALE,YSCALE

c          check for valid length.

      LENGTH = LENGT1
      IF(LENGTH.LT.1) THEN
         CALL WPSTER(212,5)
         GO TO 999
      END IF

c          if(length.gt.136) then call vberrh(213,5), and use the
c          maximum length of 136.

      IF(LENGTH.GT.136) THEN
         CALL WPSTER(213,5)
         LENGTH = 136
      ENDIF

      CTEMP='('
      LENOUT=1

c          loop through length characters.

      DO 100 I=1,LENGTH

c          check for valid chars.

c          ignore control characters, except for:
c          8 is backspace
c          10 is linefeed
c          13 is carriage return

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
               CALL WPSTER(208,5)
               GOTO 100
            ENDIF

c           finish the string, emulate the control char, and start a new one

c           send the buffered chars to the printer if there are any
            IF(LENOUT.NE.1) THEN
               CTEMP=CTEMP(1:LENOUT)//') t '
               LENOUT=LENOUT+4
               CALL WPST15(LENOUT,CTEMP)
C              reset the cp from the characters
               XCP=XCP+(LENOUT-5)*VECTOR(7)
            ENDIF

c           calculate the new current position after the control char
            XCP=XCP+DX
            YCP=YCP+DY
            CALL WPST13(0,XCP,YCP)

c           start a new string
            CTEMP='('
            LENOUT=1

         ELSE

c           Char value is 32-126 inclusive.  Put \ before these:
c              92 is \
c              40 is (
c              41 is )

            IF(CHARS(I).EQ.40.OR.CHARS(I).EQ.41.OR.CHARS(I).EQ.92) THEN
               CTEMP=CTEMP(1:LENOUT)//char(92)
               LENOUT=LENOUT+1
            ENDIF

c           now pack the chars into the buffer

            CALL CDRA2C(CHARS(I),STR)
            CTEMP=CTEMP(1:LENOUT)//STR(1:1)
            LENOUT=LENOUT+1
         ENDIF

  100 CONTINUE

c          send the chars to the printer

      CTEMP=CTEMP(1:LENOUT)//') t '
      LENOUT=LENOUT+4
      CALL WPST15(LENOUT,CTEMP)

C          reset the cp from the characters

      XCP=XCP+(LENOUT-5)*VECTOR(7)

  999 RETURN
      END
      SUBROUTINE WPSTLS(LINSTY)
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
      COMMON /WPST11/ IVECT,COORD
      CHARACTER COORD*20
      INTEGER IVECT
      CHARACTER CTEMP*30,STRL*3,STRS*3,STRG*3

      REAL VECTOR(7)
      COMMON /WPST03/ VECTOR
      REAL XSCALE,YSCALE
      COMMON /WPST08/ XSCALE,YSCALE
      REAL XCP,YCP
      COMMON /WPST05/ XCP,YCP
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

      ENTRY WPST14(LINSTY)

C CHECK FOR VALID LINSTY.
      IF(LINSTY.LT.0.OR.LINSTY.GT.5) THEN
         CALL WPSTER(401,5)
         VECTOR(4) = 0
         GOTO 999
      END IF

      CALL WPST15(2,'s ')
C GENERATE THE LINESTYLE COMMANDS
      IF(LINSTY.EQ.0) THEN
        CALL WPST15(7,'[] 0 h ')
      ENDIF

c calculate the linewidth -- it's needed below in every case

c        actual xscale is xscale*.1; linewidth=1 => .01 in NDC
         LW=VECTOR(5)
         LW=XSCALE*VECTOR(5)*.001
c        a linewidth of zero isn't good with postscript
         IF(LW.LT.1.) LW=1.

c     from here on, set; up patterns that depend on the linewidth and
c          the extra length added to the line segment
c          by the hemispherical end cap

      IF(LINSTY.EQ.1) THEN
         ILL=NINT(0.5*LW)
         IGAP=NINT(3.*LW)
         CALL CDRI2C(ILL,3,STRL)
         CALL CDRI2C(IGAP,3,STRG)
         CTEMP='['//STRL(1:3)//' '//STRG(1:3)//'] 0 h '
         CALL WPST15(14,CTEMP)

      ELSE IF(LINSTY.EQ.2) THEN
         ILL=NINT(18.*LW)
         JLL=NINT(1.5*LW)
         IGAP=NINT(3.*LW)
         CALL CDRI2C(ILL,3,STRL)
         CALL CDRI2C(JLL,3,STRS)
         CALL CDRI2C(IGAP,3,STRG)
         CTEMP='['//STRS(1:3)//' '//STRG(1:3)//' '//STRL(1:3)
     *            //' '//STRG(1:3)//'] 0 h '
         CALL WPST15(22,CTEMP)
c         call pstbuf(14,'[2 2 6 2] 0 h ')

      ELSE IF(LINSTY.EQ.3) THEN
         ILL=NINT(6.*LW)
         IGAP=NINT(7.*LW)
         CALL CDRI2C(ILL,3,STRL)
         CALL CDRI2C(IGAP,3,STRG)
         CTEMP='['//STRL(1:3)//' '//STRG(1:3)//'] 0 h '
         CALL WPST15(14,CTEMP)
c         call pstbuf(8,'[4] 0 h ')

      ELSE IF(LINSTY.EQ.4) THEN
         ILL=NINT(24.*LW)
         IGAP=NINT(18.*LW)
         CALL CDRI2C(ILL,3,STRL)
         CALL CDRI2C(IGAP,3,STRG)
         CTEMP='['//STRL(1:3)//' '//STRG(1:3)//'] 0 h '
         CALL WPST15(14,CTEMP)
c         call pstbuf(8,'[8] 0 h ')

      ELSE IF(LINSTY.EQ.5) THEN
         ILL=NINT(12.*LW)
         IGAP=NINT(10.*LW)
         CALL CDRI2C(ILL,3,STRL)
         CALL CDRI2C(IGAP,3,STRG)
         CTEMP='['//STRL(1:3)//' '//STRG(1:3)//'] 0 h '
         CALL WPST15(14,CTEMP)

      ENDIF

c     redefine the postscript current position

c     the code below is equivalent to
c      call vbvect(0,xcp,ycp)
c     but can't do it because vbvect calls vdstlw which calls this routine

      CTEMP=COORD(1:11)//' m '
      CALL WPST15(14,CTEMP)

      VECTOR(4)=LINSTY

  999 RETURN
      END
      SUBROUTINE WPSTCS(YSIZE)
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
      CHARACTER STR*6,CTEMP*10

      REAL VECTOR(7)
      COMMON /WPST03/ VECTOR
      REAL XSCALE,YSCALE
      COMMON /WPST08/ XSCALE,YSCALE

C CHECK FOR VALID YSIZE.
      IF(YSIZE.LT.0.0.OR.YSIZE.GT.1.0) THEN
         CALL WPSTER(401,5)
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

c output the postscript command

      CALL CDRI2C(IYSIZE,4,STR)
c     iysize is in tenths of device units
      CTEMP='y '//STR(1:3)//' x '
      CALL WPST15(8,CTEMP)

  999 RETURN
      END
      SUBROUTINE WPSTLW(LINWTH)
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

      INTEGER*4 KWRTFL,KRDFL,KOUTFL,KINFL,KWRDSZ,KBYTEL,KCPW,KBAUD,
     1KCOMTP
      COMMON /CDRCOM/ KWRTFL,KRDFL,KOUTFL,KINFL,KWRDSZ,KBYTEL,KCPW,
     1 KBAUD,KCOMTP
      REAL VECTOR(7)
      COMMON /WPST03/ VECTOR
      REAL XSCALE,YSCALE
      COMMON /WPST08/ XSCALE,YSCALE

C CHECK FOR VALID LINWTH.
      IF(LINWTH.LT.0.0.OR.LINWTH.GT.1.) THEN
         CALL WPSTER(401,5)
         GOTO 999
      END IF

C CONVERT LINE-WIDTH TO NDC
      LW=LINWTH*.01

C CONVERT WIDTH TO DEVICE COORDINATES AND ADD A DIGIT; NEED IT; TO HUNDREDTHS
      ILW=NINT(XSCALE*LW*10.)
C     A LINEWIDTH OF ZERO WORKS ONLY PART OF THE TIME
      IF(ILW.LT.10) ILW=10

C SET LINE WIDTH
      CALL CDRI2C(ILW,5,STR)
      CTEMP='s '//STR(1:3)//'.'//STR(4:5)//' w '
      CALL WPST15(11,CTEMP)

      VECTOR(5)=LINWTH

c     since linestyle uses the linewidth in setting the pattern, call it

      LINSTY=VECTOR(4)
      CALL WPST14(LINSTY)
  999 RETURN
      END
      SUBROUTINE WPSTIE(ESCPCD,SUPPRT)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDIQES           -Inquire Escape.

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77

C ENTRY CONDITIONS -ESCPCD = integer escape function code.

C CALLS            -

C EXIT CONDITIONS  -SUPPRT = integer level of support for the escape
C                   function specified.  Range 0,1,2.

C NARRATIVE        -An integer value indicating 2=hardware supported,
C                   1=software supported, 0=unsupported is returned in
C                   SUPPRT for the escape function ESCPCD.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      INTEGER ESCPCD,SUPPRT
      IF (ESCPCD.EQ.2100) THEN
         SUPPRT=2
C ELSE THERE IS NO SUPPORT OF ANY OTHER ESCAPE CODES
      ELSE
         SUPPRT=0
      END IF
      RETURN
      END
