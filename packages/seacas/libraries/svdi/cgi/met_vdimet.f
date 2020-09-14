C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C VDMOVA VDMONI VDGNAM VBIQDV VBIQPK VDLINA VDTEXT VDPNTA VDPOLY VDIQCP VDSTOS
C WMETMV WMETMO WMETGN WMETIV WMETQP WMETLN WMETTX WMETPT WMETPY WMETCP WMETOS

C VDIQOS VDSTFC VDSTBC VDSTIN VDSTLS VDSTLW VDSTCS VDAABU VDALOC VDABGL VDAKGL
C WMETIO WMETFC WMETBC WMETIN WMETLS WMETLW WMETCS WMETBU WMETLO WMETBL WMETKL

C VDSTLA VDINIT VDFRAM VDTERM VDIQDC VDNWPG VDBELL VDWAIT VDBUFL VDSTCO VDIQCO
C WMETLA WMETNT WMETFR WMETTR WMETDC WMETPG WMETBE WMETWT WMETFL WMETCO WMETIC

C VDESCP VDIQES VDIQND VIMOVA VILINA VIPNTA VITEXT VIINIT VITERM VINWPG CDRCOM
C WMETES WMETIE WMETID WMETIM WMETIL WMETIP WMETIX WMETII WMETIT WMETIG CDRCOM

C VCJOB  VCONOD VBERRH VDLOGE CDRWFS CDRRFS CDROFS CDROF3 CDRCFS CDROFF CDROAB
C  VCJOB VCONOD WMETER WMETLE WMETWF WMETRF WMETOF WMETO3 WMETCF WMETFF WMETAB

C BGPBUF QMSBUF QMSBU1 DDCBUF H75BUF BTKBUF NMTBUF VBIMBF VBPKG  VBDEV  VDIQRS
C WMETBF WMETQM WMETBF WMETBF WMETBF WMETBF WMETBF WMETIB WMETPK WMETDV WMETQR

C VDSTMP VDSTRS VDSTRV VDBRGB VDFRGB VDPIXL VDPIXI VDRPIX VDRPXI VDRSCL VDIQCI
C WMETMP WMETRS WMETRV WMETBG WMETFG WMETPX WMETPI WMETRP WMETRI WMETRL WMETCI

C VBSTMP VIFRAM VCNDCM VCATTR VBINI1 VB2HLS VB2RGB VCCOLT VCCRPS VCSCAL VCDDIM
C WMET01 WMET02 WMET03 WMET04 WMET05 WMET06 WMET07 WMET08 WMET09 WMET10 WMET11

C VIPOLY VBOUT
C WMET12 WMET13

      SUBROUTINE WMET01( IMAP )
      integer*4 imap

      GOTO (1,2,3,4,5),IMAP

      CALL WMETMP('UNKNOWN')
      RETURN

    1 CALL WMETMP('1-TO-1')
      RETURN

    2 CALL WMETMP('REPLICATE')
      RETURN

    3 CALL WMETMP('VIRTUAL')
      RETURN

    4 CALL WMETMP('NODISTORT')
      RETURN

    5 CALL WMETMP('FREEFORM')
      RETURN

      END
      SUBROUTINE WMETRS(I1,I2)
c*************************************************************************
c      This routine is to satisfy entry points used with raster vdi stuff
c      but not with regular vdi stuff.  This is done so raster vdi programs
c      can link with regular vdi.
c*************************************************************************
      INTEGER*4 i1, i2, i3
      CHARACTER C1*(*)
      REAL*4 RA1(1),RA2(1),RA3(1)
      INTEGER*4 IA1(1),IA2(1)
      real*4 r1, r2, r3, r4

      ENTRY WMETRV(R1,R2,R3,R4)
      ENTRY WMETMP(C1)
      ENTRY WMETPX(I1,I2,RA1,RA2,RA3,I3)
      ENTRY WMETRP(I1,I2,I3,RA1,RA2,RA3,IA1)
      ENTRY WMETPI(I1,I2,IA1,I3)
      ENTRY WMETRI(I1,I2,I3,IA1,IA2)
      ENTRY WMETQR(I1,RA1)
      ENTRY WMETRL
      ENTRY WMETFG(R1,R2,R3)
      ENTRY WMETBG(R1,R2,R3)
      ENTRY WMETCI(R1,R2,R3,I1)
      RETURN
      END
      SUBROUTINE WMETFR(ITYPE)
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

      INTEGER*4 ITYPE

      CALL WMET02(ITYPE)
      RETURN
      END
      SUBROUTINE WMET02(ITYPE)
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

      INTEGER*4 ITYPE

      RETURN
      END
      SUBROUTINE WMETGN(NAME)
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
      INTEGER*4 LENGTH,ISTART,IEND,I
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
      SUBROUTINE WMETNT(ASPECT,JUSTIF)
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

      REAL*4 ASPECT
      INTEGER*4 JUSTIF

      CALL WMETII(ASPECT,JUSTIF)

      RETURN
      END
      SUBROUTINE WMETLE(ERRNUM,ERRSEV)
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

      INTEGER ERRNUM
      INTEGER ERRSEV

      INTEGER*4 KWRTFL,KRDFL,KOUTFL,KINFL,KWRDSZ,KBYTEL,KCPW,KBAUD,
     1KCOMTP
      COMMON /CDRCOM/ KWRTFL,KRDFL,KOUTFL,KINFL,KWRDSZ,KBYTEL,KCPW,
     1 KBAUD,KCOMTP

C flush buffer before we do a write
      CALL WMETFL

C WRITE THE ERROR TO THE LISTING.
      WRITE(KWRTFL,10)ERRNUM,ERRSEV
   10 FORMAT(' SVDI ERROR NUMBER ',I5,'   SEVERITY CODE ',I5)

C TRACEBACK.
      CALL CDRTBK

      RETURN
      END
      SUBROUTINE WMETID(XNDC,YNDC)
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

      REAL*4 XNDC,YNDC

      REAL*4 XNDCMX,YNDCMX
      COMMON /WMET03/ XNDCMX,YNDCMX

C RETURN THE MAXIMUM VALID NDC VALUES.
      XNDC=XNDCMX
      YNDC=YNDCMX

      RETURN
      END
      SUBROUTINE WMETIO(ATTARR)
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

      REAL*4 ATTARR(7)

      REAL*4 VECTOR(7)
      COMMON /WMET04/ VECTOR

      INTEGER*4 I

      DO 100 I=1,7
         ATTARR(I)=VECTOR(I)
  100 CONTINUE

      RETURN
      END
      SUBROUTINE WMETLN(X,Y)
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

      REAL*4 X,Y

C THIS IS JUST A DUMMY ROUTINE WHICH DOES NOTHING BUT CALL VILINA.
C THIS ORGANIZATION FACILITATES ADDING SECURITY MARKINGS TO SVDI.
      CALL WMETIL(X,Y)

      RETURN
      END
      SUBROUTINE WMETMO(ISTATE)
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
      INTEGER ISTATE
      INTEGER*4 MPAGES
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

      ENTRY WMETPK (C1)
      MPKG = C1
      RETURN
      ENTRY WMETDV (C2)
      MDEV = C2
      RETURN
      ENTRY WMETQP(C1)
      C1 = MPKG
      RETURN
      ENTRY WMETIV(C2)
      C2 = MDEV
      RETURN
      END
      SUBROUTINE WMETMV(X,Y)
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

      REAL*4 X,Y

C THIS IS JUST A DUMMY ROUTINE WHICH DOES NOTHING BUT CALL VIMOVA.
C THIS ORGANIZATION FACILITATES ADDING SECURITY MARKINGS TO SVDI.
      CALL WMETIM(X,Y)

      RETURN
      END
      SUBROUTINE WMETPG
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
      CALL WMETIG

      RETURN
      END
      SUBROUTINE WMETPT(X,Y)
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

      REAL*4 X,Y

C THIS IS JUST A DUMMY ROUTINE WHICH DOES NOTHING BUT CALL VIPNTA.
C THIS ORGANIZATION FACILITATES ADDING SECURITY MARKINGS TO SVDI.
      CALL WMETIP(X,Y)

      RETURN
      END
      SUBROUTINE WMETPY(XARRAY,YARRAY,NPTS)
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

      INTEGER*4 NPTS
      REAL*4 XARRAY(NPTS),YARRAY(NPTS)

C THIS IS JUST A DUMMY ROUTINE WHICH DOES NOTHING BUT CALL VIPOLY.
C THIS ORGANIZATION FACILITATES ADDING SECURITY MARKINGS TO SVDI.
      CALL WMET12(XARRAY,YARRAY,NPTS)

      RETURN
      END
      SUBROUTINE WMETOS(ATTARR)
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

      REAL*4 ATTARR(6)

      REAL*4 VECTOR(7)
      COMMON /WMET04/ VECTOR
      integer*4 i4

C CALL EACH OF THE INDIVIDUAL ATTRIBUTE SETTING ROUTINES.
C CHECK FOR VALIDITY OF INPUT VALUES WILL BE DONE IN EACH INDIVIDUAL
C ROUTINE.
      i4 = INT(ATTARR(1))
      CALL WMETFC(i4)
      i4 = INT(ATTARR(2))
      CALL WMETBC(i4)
      CALL WMETIN(ATTARR(3))
      i4 = INT(ATTARR(4))
      CALL WMETLS(i4)
      CALL WMETLW(ATTARR(5))
      CALL WMETCS(ATTARR(6))

      RETURN
      END
      SUBROUTINE WMETTR
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
      CALL WMETIT

      RETURN
      END
      SUBROUTINE WMETTX(LENGTH,CHARS)
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

      INTEGER*4 LENGTH,CHARS(136)

C THIS IS JUST A DUMMY ROUTINE WHICH DOES NOTHING BUT CALL VITEXT.
C THIS ORGANIZATION FACILITATES ADDING SECURITY NARKINGS TO SVDI.
      CALL WMETIX(LENGTH,CHARS)

      RETURN
      END
      SUBROUTINE WMETER(ERRNUM,ERRSEV)
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

      INTEGER ERRNUM
      INTEGER ERRSEV

C REPORT THE ERROR USING VDLOGE.
      CALL WMETLE(ERRNUM,ERRSEV)

C CHECK FOR FATAL ERROR.
      IF(ERRSEV.GT.12) STOP

      RETURN
      END
      SUBROUTINE WMETBU(BTNNUM)
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

      INTEGER*4 BTNNUM

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
      SUBROUTINE WMETBL(BTNNUM,X,Y)
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

      REAL*4 X,Y
      INTEGER*4 BTNNUM

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
      SUBROUTINE WMETKL(CHAR,X,Y)
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

      REAL*4 X,Y
      INTEGER*4 CHAR

      INTEGER*4 IN,CHR

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
      SUBROUTINE WMETLO(X,Y)
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

      REAL*4 X,Y

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
      SUBROUTINE WMETBE
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
      SUBROUTINE WMETLA(LOCX,LOCY)
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

      REAL*4 LOCX,LOCY

C BATCH DEVICES IGNORE THIS FUNCTION.

      RETURN
      END
      SUBROUTINE WMETWT
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
      SUBROUTINE WMET05(N,NSTR)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VBINI1           -Virtual Device Initialization String Output.

C R.W.Simons       -18MAY81

C ENVIRONMENT      -Computer-independent, System-independent, FORTRAN 77
C                   Hard Copy Format 1.

C ENTRY CONDITIONS -N = integer number of words in NSTR. (max=4)
C                   NSTR = integer array containing the string to be
C                   converted and output.  The last character must
C                   be the string terminator.

C CALLS            -CDRCVT,CDR1CH,VBOUT.

C EXIT CONDITIONS  -

C NARRATIVE        -This routine converts a string from
C                   internal computer-dependent format to
C                   ASCII and sends it to the device.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      INTEGER   N
      INTEGER*4 NSTR(4)
      integer i, j
      integer*4 itemp
      integer*4 itemp1, itemp2
      integer   itemp8

      INTEGER*4 KWRTFL,KRDFL,KOUTFL,KINFL,KWRDSZ,KBYTEL,KCPW,KBAUD,
     1KCOMTP
      COMMON /CDRCOM/ KWRTFL,KRDFL,KOUTFL,KINFL,KWRDSZ,KBYTEL,KCPW,
     1 KBAUD,KCOMTP

C KTERM = STRING TERMINATOR CHARACTER (\).
      INTEGER*4 KTERM
      DATA KTERM /92/

C LOOP THROUGH EACH CHARACTER IN EACH WORD OF NSTR.
      DO I=1,N
         DO J=1,KCPW
            CALL CDR1CH(J,NSTR(I),ITEMP)

C     CONVERT CHARACTER.
            CALL CDRCVT(ITEMP,ITEMP1)
C
C     CHECK FOR END-OF-STRING CHARACTER.
            IF (ITEMP1.EQ.KTERM) GO TO 20
C
C     SEND PAIRS OF CHARACTERS TO THE OUTPUT FILE.
            IF(MOD(J,2).EQ.1) THEN
               ITEMP2=ITEMP1
            ELSE
               ITEMP8=256*ITEMP2+ITEMP1
               CALL WMET13S(ITEMP8)
            ENDIF
         END DO
      END DO
C
C PAD WITH A BLANK IF NECESSARY TO MAKE NUMBER OF CHARS EVEN.
 20   CONTINUE
      IF(MOD(J,2).EQ.0) THEN
         ITEMP8=256*ITEMP2+32
         CALL WMET13S(ITEMP8)
      ENDIF

      RETURN
      END
      SUBROUTINE WMET06(RGB,MAXVAL,HLS)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VB2HLS            - Transform RGB to HLS

C P. Watterberg    - 2 APR 81

C ENVIRONMENT      - COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77

C ENTRY CONDITIONS - RGB = integer array with 3 elements specifying
C                          Red,   RGB(1),    range 0 - MAXVAL
C                          Green,  RGB(2),    range 0 - MAXVAL
C                          Blue, RGB(3),    range 0 - MAXVAL

C                    MAXVAL = integer, largest value that each of R, G or B
C                             can assume

C CALLS            - none

C EXIT CONDITIONS  - HLS = Real array with 3 elements specifying
C                          Hue,        HLS(1),   range 0. - 360.
C                          Lightness,  HLS(2),   range 0. - 1.
C                          Saturation, HLS(3),   range 0. - 1.

C NARRATIVE        - This routine converts RGB to HLS.  The interpretation
C                    of HLS is the one adopted by GSPC 79.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      integer   maxval
      REAL*4 HLS(3), temp
      INTEGER*4 RGB(3)
      integer*4 ired, igre, iblu, maxc, minc, isum, idif
      integer*4 maxlit

C          copy the inputs to locals

      IRED = RGB(1)
      IGRE = RGB(2)
      IBLU = RGB(3)

c          compute some useful quantities

      MAXC = MAX(IRED,IGRE,IBLU)
      MINC = MIN(IRED,IGRE,IBLU)
      ISUM = MAXC + MINC
      IDIF = MAXC - MINC
      MAXLIT = 2*MAXVAL

c          getting lightness is easy

      HLS(2) = DBLE(ISUM)/DBLE(MAXLIT)

c          getting saturation is a little more difficult

      IF(IDIF.EQ.0) THEN
          HLS(3) = 0.
        ELSE
          IF(ISUM.LE.MAXVAL) THEN
              HLS(3) = DBLE(IDIF)/ISUM
            ELSE
              HLS(3) = DBLE(IDIF)/(MAXLIT-ISUM)
            ENDIF
        ENDIF

c          getting hue is a little harder yet

      IF(IDIF.EQ.0) THEN
          TEMP = 0.
        ELSE
          TEMP = 60./IDIF
        ENDIF

c                     is the maximum color red?

      IF(MAXC.EQ.IRED) THEN
          HLS(1) = 120. + ((IGRE-MINC) - (IBLU-MINC))*TEMP

c                     is it green?

        ELSE IF(MAXC.EQ.IGRE) THEN
          HLS(1) = 240. + ((IBLU-MINC) - (IRED-MINC))*TEMP

c                     well then, it must be blue

        ELSE
          IF(IRED.GE.IGRE) THEN
              HLS(1) = (IRED-MINC)*TEMP
            ELSE
              HLS(1) = 360. - (IGRE-MINC)*TEMP
            ENDIF
        ENDIF

      RETURN
      END
      SUBROUTINE WMET07(HLS,RGB,MAXVAL)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VB2RGB           - Transform HLS to RGB

C P. Watterberg    - 30 MAR 81

C ENVIRONMENT      - COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77

C ENTRY CONDITIONS - HLS = Real array with 3 elements specifying
C                          Hue,        HLS(1),   range 0. - 360.
C                          Lightness,  HLS(2),   range 0. - 1.
C                          Saturation, HLS(3),   range 0. - 1.

C                    MAXVAL = integer, largest value that any of R, G or B
C                             can assume

C CALLS            - none

C EXIT CONDITIONS  - RGB = integer array with 3 elements specifying
C                          Red,   RGB(1),    range 0 - MAXVAL
C                          Green,  RGB(2),    range 0 - MAXVAL
C                          Blue , RGB(3),    range 0 - MAXVAL

C NARRATIVE        - This routine converts HLS to RGB.  The interpretation
C                    of HLS is the one adopted by GSPC 79.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      INTEGER MAXVAL
      REAL*4 HLS(3), LIT, HUE, SAT, f
      INTEGER*4 RGB(3), inten, irange, irang2, isplus
      integer*4 iv1, iv2, iv3, iv4, ib, ig, ir
      integer*4 ijump, zero
      integer*4 imaxval

      zero = 0
      imaxval = maxval

C          copy the inputs to locals

      HUE = HLS(1)
      LIT = HLS(2)
      SAT = HLS(3)

C          find out which major hue (0 - 5) we are interested in

      HUE = HUE/60.

C          avoid the maximum boundary conditions

      IF(HUE.GE.6.) HUE = 5.99
      IF(SAT.GE.1.) SAT = .99
      IF(LIT.GE.1.) LIT = .99

C          the conversions and convolutions that happen here are not
C          very easy to understand.  It's best to talk to Peter but
C          if you need to try to decipher it yourself, you might try
C          by first assuming a saturation of 1.  That way, irang2=irange,
C          isplus goes away and you are left with the outer shell of
C          rgb color cube to deal with.

C          ijump represents one of the six edges of the color cube
C          connecting the six major hues.

      IJUMP = HUE

c          f is the distance (0.-.999) along an edge between two major hues

      F = HUE - IJUMP
      INTEN = LIT*DBLE(2*IMAXVAL+1)

c          irange is the range a color may take on (i.e. maxval adjusted for
c                 intensity
c          irang2 is irange adjusted for saturation

      IRANGE = IMAXVAL - ABS(INTEN-IMAXVAL)
C ... This is done for the 8-byte systems so we can pass native int to mod intrinsic
      IRANGT = IRANGE
      IRANG2 = 2*(INT((IRANGT/2+1)*SAT)) + MOD(IRANGT,2)

c          isplus is an additive to account for saturation

      ISPLUS = (IRANGE-IRANG2)/2
      IV1 = MIN(INTEN,IMAXVAL) - ISPLUS
      IV2 = MAX(zero,INTEN-IMAXVAL) + ISPLUS
      IV3 = F*IRANG2 + .5 + IV2
      IV4 = (1.-F)*IRANG2 + .5 + IV2
      GOTO (610,620,630,640,650,660),IJUMP+1

  610 IB = IV1
      IG = IV2
      IR = IV3
      GOTO 670

  620 IR = IV1
      IG = IV2
      IB = IV4
      GOTO 670

  630 IR = IV1
      IB = IV2
      IG = IV3
      GOTO 670

  640 IG = IV1
      IB = IV2
      IR = IV4
      GOTO 670

  650 IG = IV1
      IR = IV2
      IB = IV3
      GOTO 670

  660 IB = IV1
      IR = IV2
      IG = IV4

  670 RGB(1) = IR
      RGB(2) = IG
      RGB(3) = IB

      RETURN
      END
      SUBROUTINE WMETFL
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDBUFL           -Buffer Flush.

C K.M.ERICKSON     -04 MAY 81

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77

C ENTRY CONDITIONS -

C CALLS            -VBOUT

C EXIT CONDITIONS  -

C NARRATIVE        -Assure that the picture is up-to-date by flushing
C                   buffers if necessary.  This is necessary to
C                   guarantee that the picture is in a certain state
C                   before interacting with it.
C                  -set terminal to alpha mode in order to facilitate
c                   fortran IO.
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

c  flush buffers
c modified 2-23-87 to be a dummy routine by JFM.
C**      CALL VBOUT(0,1)
      RETURN
      END
      SUBROUTINE WMETIC(NUM,INDEX,CLRARY,CLRMOD)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDIQCO           - Inquire Color Table.

C K.M.Erickson     - 04 May 81

C ENVIRONMENT      - COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77

C ENTRY CONDITIONS - NUM = integer number of color indexes to inquire.
C                    Range 1-256.
C                    INDEX = integer array of indexes to inquire.  Range
C                    0-255.
C                    CLRMOD = integer color model to be used.  Range 0,1.

C CALLS            - vb2hls

C EXIT CONDITIONS  - CLRARY = real array of 3 by NUM elements returning
C                    the values of the components of the indexes inquired.
C                    Range for rgb:  0. - 1.
C                    Range for hls:  hue   0. - 360.
C                                    l & s 0. - 1.

C NARRATIVE        - Inquire one or more color table entries.  NUM and
C                    INDEX specify how many and which indexes are being
C                    inquired.  CLRMOD specifies which color model
C                    (0=RGB, 1=HLS) should be used in constructing values
C                    to return in CLRARY.  A device which does not
C                    support the color table indexes specified will
C                    return -1. in the first element of the CLRARY value
C                    for that index.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      INTEGER*4  NUM,INDEX(NUM),CLRMOD, RGB(3)
      REAL*4 CLRARY(3,NUM)

      INTEGER*4 CLRTAB(256,3)
      COMMON /WMET08/ CLRTAB

C          check for valid num.

      integer*4 i, indexn

      IF(NUM.LT.1.OR.NUM.GT.256) THEN
         CALL WMETER(723,5)
         GOTO 999
      END IF

C check for valid clrmod.

      IF(CLRMOD.NE.0.AND.CLRMOD.NE.1) THEN
         CALL WMETER(725,5)
         GOTO 999
      END IF

C          check for valid indexes.

      DO I=1,NUM
         INDEXN=INDEX(I)
         IF(INDEXN.LT.0.OR.INDEXN.GT.255) THEN
            CALL WMETER(724,5)
            GOTO 100
         END IF
         RGB(1) =CLRTAB(INDEX(I)+1,1)
         RGB(2) =CLRTAB(INDEX(I)+1,2)
         RGB(3) =CLRTAB(INDEX(I)+1,3)
         IF(CLRMOD.EQ.0) THEN
             CLRARY(1,I) = RGB(1)/255.
             CLRARY(2,I) = RGB(2)/255.
             CLRARY(3,I) = RGB(3)/255.
           ELSE
             CALL WMET06(RGB,255,CLRARY(1,I))
           ENDIF
 100       continue
        end do
  999 RETURN
      END
      SUBROUTINE WMETCP(X,Y)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDIQCP           -Inquire Where Current Position Is.

C K.M. ERICKSON    -14 NOV 80

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77

C ENTRY CONDITIONS -

C CALLS            -

C EXIT CONDITIONS  -X,Y = real NDC position.

C NARRATIVE        -Return the value of current position.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C
      COMMON /WMET09/XCP,YCP
      REAL*4 XCP,YCP
      REAL*4 X,Y

C     ASSIGN THE CP TO X,Y

      X=XCP
      Y=YCP
      RETURN
      END
      SUBROUTINE WMETCO(NUM,INDEX,CLRARY,CLRMOD)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDSTCO           - Set Color Table.

C K.M.ERICKSON     -04 MAY 81

C ENVIRONMENT      - COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77

C ENTRY CONDITIONS - NUM = integer number of color indexes to be set.
C                    Range 1-256.
C                    INDEX = integer array of indexes to be set.  Range
C                    0-255.
C                    CLRARY = real array of 3 by NUM elements specifying
C                    the values of the components of the index to be
C                    set.
C                    Range for RGB: red   0.0 - 1.0
C                                   green 0.0 - 1.0
C                                   blue  0.0 - 1.0
C                    Range for HLS: hue        0.0 - 360.0
C                                   lightness  0.0 -   1.0
C                                   saturation 0.0 -   1.0
C                    Default:
C                    Index  Color  RGB Values
C                      0    black   0.,0.,0.
C                      1    red     1.,0.,0.
C                      2    green   0.,1.,0.
C                      3    yellow  1.,1.,0.
C                      4    blue    0.,0.,1.
C                      5    magenta 1.,0.,1.
C                      6    cyan    0.,1.,1.
C                      7    white   1.,1.,1.
C                    CLRMOD = integer color model being used.  Range 0,1.
C                    Default: 0 (RGB).

C CALLS            - VBOUT

C EXIT CONDITIONS  - The Dicomed color table has been set

C NARRATIVE        - Set one or more color table entries.  This is a
C                    dynamic setting, if the device will support it.
C                    "Dynamic" neans that primitives which have already
C                    been drawn will change their appearance when a
C                    dynamic setting is changed.  INDEX is the
C                    position (or array of positions) in the table
C                    (0-255).  CLRARY is a three-element vector (or 3 by
C                    NUM array) with the fractions (0.-1.) of RGB or
C                    values (0.0-360.0, 0.0-1.0, 0.0-1.0) of HLS.
C                    A translator for HLS to RGB will be used from
C                    GSPC79.  CLRMOD specifies which color model is being
C                    used (0=RGB, 1=HLS).
C                    All devices must support at least a single device
C                    dependent value for each of red, green, and blue and
C                    the corresponding values for hue, lightness, and
C                    saturation.  If unsupported values are specified,
C                    set to the closest supported values.
C                    All devices must support both RGB and HLS color
C                    models.
C                    All devices must support at least a single device
C                    dependent INDEX value in the range 0-7.  If an
C                    unsupported value is
C                    specified, it should be ignored.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      INTEGER*4 NUM, INDEX(NUM), CLRMOD, RGB(3), i
      REAL*4 CLRARY(3,NUM)

      INTEGER*4 CLRTAB(256,3)
      COMMON /WMET08/ CLRTAB

c  batch update mode--c700
      INTEGER IBATUP

c  send color table--c800
      INTEGER ISNDCO

      INTEGER*4 IBUF(6)
      DATA IBUF/36869,5*0/
c     ibuf  (1)     (2)      (3)    (4)     (5)      (6)
c          9005    index     z R    g Y     B M      C W

      DATA IBATUP/50944/
      DATA ISNDCO/51200/

c          check for valid NUM.

      IF(NUM.LT.1.OR.NUM.GT.256) THEN
         CALL WMETER(723,5)
         GOTO 999
      END IF

c          check for valid clrmod.

      IF(CLRMOD.NE.0.AND.CLRMOD.NE.1) THEN
         CALL WMETER(725,5)
         GOTO 999
      END IF

c   send batch update mode
      CALL WMET13S(IBATUP)
c          check for valid indexes.

      DO 100 I=1,NUM
         IF(INDEX(I).LT.0.OR.INDEX(I).GT.255) THEN
            CALL WMETER(724,5)
            GOTO 100
         END IF

         IBUF(2) = INDEX(I)
c          check for valid clrary.
c rgb

         IF(CLRMOD.EQ.0) THEN
            IF( CLRARY(1,I).LT.0..OR.CLRARY(1,I).GT.1.
     X      .OR.CLRARY(2,I).LT.0..OR.CLRARY(2,I).GT.1.
     X      .OR.CLRARY(3,I).LT.0..OR.CLRARY(3,I).GT.1.) THEN
               CALL WMETER(727,5)
               GOTO 100
            END IF
            IBUF(3)=INT(255.99*CLRARY(1,I))
c                             zero/red
            IBUF(4)=INT(255.99*CLRARY(2,I))*256
c                              green/yellow
            IBUF(5)=INT(255.99*CLRARY(3,I)) *256
c                              blue/magenta

         ELSE
c hls

            IF( CLRARY(1,I).LT.0..OR.CLRARY(1,I).GT.360.
     X      .OR.CLRARY(2,I).LT.0..OR.CLRARY(2,I).GT.1.
     X      .OR.CLRARY(3,I).LT.0..OR.CLRARY(3,I).GT.1.) THEN
               CALL WMETER(727,5)
               GOTO 100
            END IF
            CALL WMET07(CLRARY(1,I),RGB,15)
            IBUF(3)=INT(255.99*RGB(1))
c                            red
            IBUF(4)=INT(255.99*RGB(2))*256
c                            green
            IBUF(5)=INT(255.99*RGB(3))*256
c                            blue

         END IF
c  store color table values

         CLRTAB(IBUF(2)+1,1)=IBUF(3)
         CLRTAB(IBUF(2)+1,2)=IBUF(4)/256
         CLRTAB(IBUF(2)+1,3)=IBUF(5)/256

      CALL WMET13(6,IBUF)
  100 CONTINUE

  999 CALL WMET13S(ISNDCO)
c  send color table
      RETURN
      END
      SUBROUTINE WMETFC(COLOR)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDSTFC           -Set Foreground Color.

C K.M. ERICKSON    -12 NOV 80

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77

C ENTRY CONDITIONS -COLOR = integer color table index . Range 0-255.
C                   Default is device dependent, in range 0-7.

C CALLS            -

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

      INTEGER*4 COLOR

      INTEGER*4 COL(2)

      REAL*4 VECTOR(7)
      COMMON /WMET04/ VECTOR

      DATA COL/37121,0/

C CHECK FOR VALID COLOR.
      IF(COLOR.LT.0.OR.COLOR.GT.255) THEN
         CALL WMETER(724,5)
         GOTO 999
      END IF

C      VECTOR(1)=FOREGROUND COLOR
C            (2)=BACKGROUND COLOR
C            (3)=INTENSITY
C            (4)=LINE STYLE
C            (5)=LINE WIDTH
C            (6)=CHARACTER BOX Y
C            (7)=CHARACTER BOX X
      VECTOR(1)=COLOR
      COL(2)=COLOR
      CALL WMET13(2,COL)
999   RETURN
      END
      SUBROUTINE WMETLW(LINWTH)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDSTLW           -Set Line Width.

C K.M. ERICKSON    -14 NOV 80

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

      REAL*4 LINWTH

      INTEGER*4 LW(2)

      REAL*4 VECTOR(7)
      COMMON /WMET04/ VECTOR
C CURRENT ATTRIBUTE VALUES. (DEVICE-INDEPENDENT)
C      REAL*4 VECTOR(7)
C      COMMON /VCATTR/ VECTOR
C     VECTOR(1)=FOREGROUND COLOR
C           (2)=BACKGROUND COLOR
C           (3)=INTENSITY
C           (4)=LINE STYLE
C           (5)=LINE WIDTH
C           (6)=CHARACTER BOX Y
C           (7)=CHARACTER BOX X

      DATA LW/43265,0/

C CHECK FOR VALID LINWTH.
      IF(LINWTH.LT.0.0.OR.LINWTH.GT.1.0) THEN
         CALL WMETER(401,5)
         GOTO 999
      END IF

C  ASSIGN VECTOR(5)

      VECTOR(5)=LINWTH

C MAP 0.-1. INTO 0.-32767; FOR DICOMED
C 1. IS .01 IN NDC SPACE

C      LW(2)=LINWTH*32767*.01
      LW(2)=LINWTH*32767

C  SEND BGP COMMAND

      CALL WMET13(2,LW)

  999 RETURN
      END
      SUBROUTINE WMETIG
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VINWPG           -New Page.

C R.W.Simons       -15MAY81

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   Hard Copy Format 1.

C ENTRY CONDITIONS -

C CALLS            -VBOUT.

C EXIT CONDITIONS  -

C NARRATIVE        -Physically advance the medium or clear the screen,
C                   whichever is appropriate.  Also flood the screen
C                   with the background color on devices that support
C                   this.  The current position is not changed.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C OUTARY = NEW PAGE COMMAND
C        = 8505,FFFF,FFFF,FFFF,PAGE NUMBER,0 IN HEX
C        = 34053,65535,65535,65535,PAGE NUMBER,0
      INTEGER*4 OUTARY(6)
      DATA OUTARY /34053,65535,65535,65535,0,0/

C SEND A NEW PAGE COMMAND TO THE PLOT FILE.
C INCREMENT THE PAGE NUMBER.
      OUTARY(5)=OUTARY(5)+1
      CALL WMET13(6,OUTARY)
      CALL WMETMO(1)

      RETURN
      END
      SUBROUTINE WMETIM(X,Y)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VIMOVA           -Move Absolute.

C R.W.Simons       -15MAY81

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   Hard Copy Format 1.

C ENTRY CONDITIONS -X,Y = real NDC position.

C CALLS            -VBOUT.

C EXIT CONDITIONS  -XCP,YCP = real updated current position. (XNDC,YNDC)

C NARRATIVE        -Set current position to absolute NDC position X,Y.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      REAL*4 X,Y

      INTEGER*4 OUTARY(2)

C CURRENT POSITION. (LXY,HC1)
      REAL*4 XCP,YCP
      COMMON /WMET09/ XCP,YCP
C      Include '[VDIMAINT.COMMON]VCCRPS'

C SCALE FACTORS FOR NDC TO DC MAPPING. (LXY)
      REAL*4 XSCALE,YSCALE
      COMMON /WMET10/ XSCALE,YSCALE

C DC DIMENSIONS OF OFFSETS AND PICTURE. (LXY,HC1)
      REAL*4 XPAD,YPAD,XDEVIC,YDEVIC
      COMMON /WMET11/ XPAD,YPAD,XDEVIC,YDEVIC

C CONVERT TO SCREEN UNITS.
C SET BIT 15 OF X = 0 TO INDICATE A COORDINATE POSITIONING COMMAND.
C SET BIT 15 OF Y = 0 TO INDICATE A DRAW COMMAND.
C ASSUME X AND Y ARE < 32768, WHICH WILL BE TRUE IF XNDC AND YNDC
C ARE IN THE PROPER RANGE.
      OUTARY(1)=X*XSCALE+XPAD
      OUTARY(2)=Y*YSCALE+YPAD

C SEND THE COMMAND TO THE PLOT FILE.
      CALL WMET13(2,OUTARY)

C UPDATE CURRENT POSITION.
      XCP=X
      YCP=Y

      RETURN
      END
      SUBROUTINE WMETES(ESCPCD,N,ARGS)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDESCP           -Escape Code Routine.

C K.M. ERICKSON    -4 MAY 81

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77

C ENTRY CONDITIONS -ESCPCD = integer escape function code.
C                   N = integer number of arguments in ARG.  RANGE 0-.
C                   ARGS = real array of arguments for the escape
C                   function specified.

C CALLS            -vbout

C EXIT CONDITIONS  -

C NARRATIVE        -Invoke the nonstandard, device-dependent function
C                   ESCPCD.  N is the number of arguments used by this
C                   function and ARGS is a real array containing those
C                   arguments.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      INTEGER*4 ESCPCD,N, i
C N COULD BE EQUAL TO 0 SO:
      REAL*4 ARGS(N+1), ARG
      REAL*4 ONE

      INTEGER*4 IBUF(4)
      ONE = 1.0

C CHECK FOR VALID N.
      IF(N.LT.0) THEN
         CALL WMETER(802,5)
         GOTO 999
      END IF

C          meta file escapes 800 -

      IF(ESCPCD.EQ.800) THEN

C SEND ASPECT RATIO
         IBUF(1) = 33539
         IBUF(2) = 32767.*MIN(ARGS(1),ONE)
         IBUF(3) = 32767./MAX(ARGS(1),ONE)
         IBUF(4) = 0
         CALL WMET13(4,IBUF)
         ENDIF

C ALL OTHER ESCAPE CODES

C FIRST CHECK IF THIS IS AN ESCAPE WITH ALPHA ARGUMENTS

      IF(ESCPCD.EQ.250) THEN
      ELSE

C          ALL OTHER ESCAPES HAVE REAL NUMBER ARGUMENTS
C          SEND  82xx  01nn  ESCPCD  ARGS(1) ...  ARGS(N)
C          xx is 2*(N+1) and nn is 2*N+1
C          and each arg is sent as a fixed point real with sixteen bits
C          integer and sixteen bits fraction.

      IBUF(1) = 33280 + 2*(N+1)
      IBUF(2) = 256 + 2*N + 1
      IBUF(3) = ESCPCD
      CALL WMET13(3,IBUF)
      DO 10 I=1,N
         ARG = ARGS(I)
         IBUF(1) = ARG
         IBUF(2) = (ARG-IBUF(1))*32768.
         IF(ARGS(I).LT.0.) IBUF(1) = IBUF(1) + 32768
         CALL WMET13(2,IBUF)
   10 CONTINUE
      ENDIF

999   RETURN
      END
      SUBROUTINE WMET12(XARRAY,YARRAY,NPTS)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VIPOLY           -POLYGON FILL ROUTINE

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77

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

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      integer*4 i, nn, npts
      REAL*4 XARRAY(NPTS),YARRAY(NPTS)
      REAL*4 ATTARR(7)

C MAX NPTS IS 508.  Constraint imposed by postprocessor.
      INTEGER*4 OUTARY(2)
      INTEGER*4 zero, i4

C SCALE FACTORS FOR NDC TO DC MAPPING.  (LXY)
      REAL*4 XSCALE,YSCALE
      COMMON /WMET10/ XSCALE,YSCALE

C DC DIMENSIONS OF OFFSETS AND PICTURE.  (LXY,HC1)
      REAL*4 XPAD,YPAD,XDEVIC,YDEVIC
      COMMON /WMET11/ XPAD,YPAD,XDEVIC,YDEVIC

C CHECK FOR VALID N
      IF (NPTS.LT.1) THEN
         CALL WMETER(802,5)
         GO TO 999
      END IF

C SAVE CURRENT ATTRIBUTES
      CALL WMETIO(ATTARR)

C SET CURRENT LINESTYLE TO SOLID
      zero = 0
      CALL WMETLS(zero)

C BEGIN POLYGON COMMAND = AA00
C                       = 43520
C END POLYGON COMMAND = AB00
C                     = 43776

      CALL WMET13S(43520)

      NN=NPTS
C CHECK MAXIMUM POINTS LIMIT
      IF (NN.GT.508) NN=508
C CONVERT EACH X,Y TO SCREEN UNITS AND WRITE OUT
      DO 100 I = 1,NN
        OUTARY(1) = XARRAY(I) * XSCALE + XPAD
        OUTARY(2) = YARRAY(I) * YSCALE + YPAD
        CALL WMET13(2,OUTARY)
 100  CONTINUE

      CALL WMET13S(43776)

C MOVE SOMEWHERE TO UPDATE CURRENT POSITION
      CALL WMETIM(XARRAY(1),YARRAY(1))

C RESTORE LINESTYLE
      i4 = ATTARR(4)
      CALL WMETLS(i4)

999   RETURN
      END
      SUBROUTINE WMET13(NUMWDS,OUTARY)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VBOUT            -Output 16 Bits of Data.

CC ENVIRONMENT      -Computer-independent, System-independent, FORTRAN 77

C ENTRY CONDITIONS -NUMWDS = integer number of words in OUTARY.
C                          = 0 means flush the buffer.
C                   OUTARY = integer array of output data, 16 bits/word,
C                   right-justified.

C NARRATIVE        - This routine used to do all the work but due to
C                    complex computer, device and software (COMDQ)
C                    dependencies, the work has moved to the computer
C                    dependent, device dependent, COMQ dependent routine
C                    BGPBUF.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C DIMENSION OUTARY TO NUMWDS+1 TO AVOID PROBLEMS WHEN NUMWDS = 0.
      INTEGER   NUMWDS
      INTEGER*4 OUTARY(NUMWDS+1)

      CALL WMETBF(NUMWDS,OUTARY)
      RETURN
      END
      SUBROUTINE WMET13S(OUT)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VBOUT            -Output 16 Bits of Data.

CC ENVIRONMENT      -Computer-independent, System-independent, FORTRAN 77

C ENTRY CONDITIONS -NUMWDS = integer number of words in OUTARY.
C                          = 0 means flush the buffer.
C                   OUTARY = integer array of output data, 16 bits/word,
C                   right-justified.

C NARRATIVE        - This routine used to do all the work but due to
C                    complex computer, device and software (COMDQ)
C                    dependencies, the work has moved to the computer
C                    dependent, device dependent, COMQ dependent routine
C                    BGPBUF.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C DIMENSION OUTARY TO NUMWDS+1 TO AVOID PROBLEMS WHEN NUMWDS = 0.
      INTEGER*4 OUTARY(1)
      INTEGER OUT

      OUTARY(1) = OUT
      CALL WMETBF(1,OUTARY)
      RETURN
      END
      SUBROUTINE WMETDC(INDEX,VALUE)
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

      INTEGER*4 INDEX
      REAL*4 VALUE
      REAL*4 DEV(33)
      SAVE DEV

C     SET DEVICE CAPABILITIES
C     THE VALUES CONTAINED IN DEV ARE:

c ** Jan 16, 1991 -- Dino Pavlakos
c       changed polygon support level (entry# 24) from 2 to 3

       DATA DEV/ 0.,0.,256.,256.,4096.,31.,32767.,0.,0.,0.,
     *          0.,0.,0.,0.,32767.,32767.,0.0,0.0,8.,8.,
     *          0.,0.,11.,3.,508.,1.,16777216.,0.,0.,21298.,
     *          32767.,1.,0./

C CHECK FOR VALID INDEX.
      IF(INDEX.LT.1.OR.INDEX.GT.33) THEN
         CALL WMETER(726,5)
         GOTO 999
      END IF

C RETURN INDEXED VALUE.
      VALUE=DEV(INDEX)

  999 RETURN
      END
      SUBROUTINE WMETIE(ESCPCD,SUPPRT)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDIQES           -Inquire Escape.

C K.M.ERICKSON     -04 may 81

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77

C ENTRY CONDITIONS -ESCPCD = integer escape function code.

C CALLS            -

C EXIT CONDITIONS  -SUPPRT = integer level of support for the escape
C                   function specified.  Range 0,1,2.

C NARRATIVE        -An integer value indicating 2=hardware supported,
C                   1=software supported, 0=unsupported is returned in
C                   SUPPRT for the escape function ESCPCD.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      INTEGER*4 ESCPCD,SUPPRT
      IF(ESCPCD .GE. 200 .AND. ESCPCD .LE. 205)THEN
              SUPPRT=2
C ELSE THERE IS NO SUPPORT OF ANY OTHER ESCAPE CODES
        ELSE
              SUPPRT=0
      END IF
      RETURN
      END
      SUBROUTINE WMETBC(COLOR)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDSTBC           -Set Background Color.

C K.M. ERICKSON    -12 NOV 80

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77

C ENTRY CONDITIONS -COLOR = integer color table index. Range 0-255.
C                   Default: device dependent, in range 0-7.

C CALLS            -

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

      INTEGER*4 COLOR
      INTEGER*4 SETBC(2)

      REAL*4 VECTOR(7)
      COMMON /WMET04/ VECTOR

      DATA SETBC(1) / 37633 /

C CHECK FOR VALID COLOR.
      IF(COLOR.LT.0.OR.COLOR.GT.255) THEN
         CALL WMETER(724,5)
         GOTO 999
      END IF

C      VECTOR(1)=FOREGROUND COLOR
C            (2)=BACKGROUND COLOR
C            (3)=INTENSITY
C            (4)=LINE STYLE
C            (5)=LINE WIDTH
C            (6)=CHARACTER BOX Y
C            (7)=CHARACTER BOX X

      VECTOR(2)=COLOR
      SETBC(2) = COLOR
      CALL WMET13(2,SETBC)

999   RETURN
      END
      SUBROUTINE WMETCS(YSIZE)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDSTCS           -Set Character Size.

C K.M.ERICKSON     -04 May 81

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77

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

      REAL*4 YSIZE

      INTEGER*4 IBUF(4)

      REAL*4 VECTOR(7)
      COMMON /WMET04/ VECTOR

      COMMON/WMET10/XSCALE,YSCALE
      REAL*4 XSCALE,YSCALE

C CHECK FOR VALID YSIZE.
      IF(YSIZE.LT.0.0.OR.YSIZE.GT.1.0) THEN
         CALL WMETER(401,5)
         GOTO 999
      END IF

C      VECTOR(1)=FOREGROUND COLOR
C            (2)=BACKGROUND COLOR
C            (3)=INTENSITY
C            (4)=LINE STYLE
C            (5)=LINE WIDTH
C            (6)=CHARACTER BOX Y
C            (7)=CHARACTER BOX X

C  SET CHARACTER BOX = SPACING OF LETTERS TO A 5/7 BOX
      VECTOR(6)=YSIZE
      VECTOR(7)=VECTOR(6)*(5./7.)

C   SEND 15 BITS FOR THE HEIGHT AND WIDTH TO BGP. NOTE THAT CHARACTER SIZES
C   ARE MAPPED FROM THE SMALLEST TO THE LARGEST CHARACTER DEFINED TO BE A
C   CHARACTER FILLING THE SMALLEST DIMENSION OF THE SCREEN ASPECT RATIO.

C  SEND BGP COMMAND ,B202-HEIGHT-WIDTH
      IBUF(1)=45570
      IBUF(2)=.65*VECTOR(6)*YSCALE
      IBUF(3)=.65*VECTOR(7)*XSCALE
      CALL WMET13(3,IBUF)

999   RETURN
      END
      SUBROUTINE WMETIN(INTEN)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDSTIN           -Set Intensity.

C K.M. ERICKSON    -12 NOV 80

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77

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

      REAL*4 INTEN

      INTEGER*4 INTE(2)

      REAL*4 VECTOR(7)
      COMMON /WMET04/ VECTOR

      DATA INTE/37377,0/
C CHECK FOR VALID INTEN.
      IF(INTEN.LT.0.0.OR.INTEN.GT.1.0) THEN
         CALL WMETER(401,5)
         GOTO 999
      END IF

C     VECTOR(1)=FOREGROUND COLOR
C           (2)=BACKGROUND COLOR
C           (3)=INTENSITY
C           (4)=LINE STYLE
C           (5)=LINE WIDTH
C           (6)=CHARACTER BOX Y
C           (7)=CHARACTER BOX X

      VECTOR(3)=INTEN
C MAP INTEN VALUE OF 0.-1. INTO 0.-32767. (15 BITS OF INFO)
      INTE(2)=INTEN*32767.

      CALL WMET13(2,INTE)

999   RETURN
      END
      SUBROUTINE WMETLS(LINSTY)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VDSTLS           -Set Line Style.

C K.M. ERICKSON    -12 NOV 80

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

      INTEGER*4 LINSTY

      INTEGER*4 LS(2)

      INTEGER*4 LINS(6)
      REAL*4 VECTOR(7)
      COMMON /WMET04/ VECTOR
C CURRENT ATTRIBUTE VALUES. (DEVICE-INDEPENDENT)
C      REAL*4 VECTOR(7)
C     COMMON /VCATTR/ VECTOR

C     VECTOR(1)=FOREGROUND COLOR
C           (2)=BACKGROUND COLOR
C           (3)=INTENSITY
C           (4)=LINE STYLE
C           (5)=LINE WIDTH
C           (6)=CHARACTER BOX Y
C           (7)=CHARACTER BOX X

      DATA LS/43009,0/
      DATA LINS/32767,0,16382,5461,27305,16383/

C CHECK FOR VALID LINSTY.
      IF(LINSTY.LT.0.OR.LINSTY.GT.5)THEN
            CALL WMETER(401,5)
            GOTO 99999
            END IF
      VECTOR(4)=LINSTY

C  LINE STYLE (15 BITS) VARIES FROM 0. : DOTTED , 1-32755 :DASHED, 32767: SOLID

      LS(2)=LINS(LINSTY+1)

C  SEND BGP COMMAND

      CALL WMET13(2,LS)
99999 RETURN
      END
      SUBROUTINE WMETII(ASPECT,JUSTIF)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VIINIT           -Initialize SVDI.   Metafile

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   DICOMED

C ENTRY CONDITIONS -ASPECT = real ratio of X dimension to Y dimension.
C                   Range >0.0.  Default: 0. (device dependent).
C                   JUSTIF = integer justification of NDC space on the
C                   device.  Range 0-9.  Default: 0 (device dependent.)

C CALLS            -VBERRH,CDROFS,VBOUT,VIMOVA.

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

      REAL*4 ASPECT
      INTEGER*4 JUSTIF
      integer     iidsiz, iusrsz, iszrou
      integer*4 i, just, aspe
      integer istat

      INTEGER*4 CLRTAB(256,3)
      COMMON /WMET08/ CLRTAB

C MAXIMUM VALID NDC VALUES. (DEVICE-INDEPENDENT)
      REAL*4 XNDCMX,YNDCMX
      COMMON /WMET03/ XNDCMX,YNDCMX

C CURRENT ATTRIBUTE VALUES. (DEVICE-INDEPENDENT)
C     VECTOR(1)=FOREGROUND COLOR
C           (2)=BACKGROUND COLOR
C           (3)=INTENSITY
C           (4)=LINE STYLE
C           (5)=LINE WIDTH
C           (6)=CHARACTER BOX Y
C           (7)=CHARACTER BOX X
      REAL*4 VECTOR(7)
      COMMON /WMET04/ VECTOR

C CURRENT POSITION. (LXY,HC1)
      REAL*4 XCP,YCP
      COMMON /WMET09/ XCP,YCP
C      Include '[VDIMAINT.COMMON]VCCRPS'

C SCALE FACTORS FOR NDC TO DC MAPPING. (LXY,HC1)
      REAL*4 XSCALE,YSCALE
      COMMON /WMET10/ XSCALE,YSCALE
C      Include '[VDIMAINT.COMMON]VCSCAL'

      INTEGER*4 KWRTFL,KRDFL,KOUTFL,KINFL,KWRDSZ,KBYTEL,KCPW,KBAUD,
     1KCOMTP
      COMMON /CDRCOM/ KWRTFL,KRDFL,KOUTFL,KINFL,KWRDSZ,KBYTEL,KCPW,
     1 KBAUD,KCOMTP

C DC DIMENSIONS OF OFFSETS AND PICTURE. (LXY,HC1)
      REAL*4 XPAD,YPAD,XDEVIC,YDEVIC
      COMMON /WMET11/ XPAD,YPAD,XDEVIC,YDEVIC

      INTEGER*4 MACHIN(3),MACLEN
      INTEGER*4 KIDSIZ,KJOBID(4),KUSRSZ,KUSRID(4),KSZROU
      INTEGER*4 KJROUT(4),KSECUR,KJTIME(4),KJDATE(4)
      COMMON / VCJOB/ KIDSIZ,KJOBID,KUSRSZ,KUSRID,KSZROU,
     1               KJROUT,KSECUR,KJTIME,KJDATE,MACHIN,MACLEN

C DECLARE FILE INITIALIZATION COMMANDS.
      integer*4 idfile(2),isecur(3), zero(1)
      real*4 rzero

      DATA ISECUR/33282,1,0/

      RZERO= 0
      ZERO(1) = 0
      XPAD = 0
      YPAD = 0

C SET DEFAULT ATTRIBUTE VALUES.  ALL ARE DEVICE DEPENDENT EXCEPT
C VECTOR(4)=0.0.
      VECTOR(1)=7.
C           FOREGROUND COLOR - white
      VECTOR(2)=0.
C           BACKGROUND COLOR - black
      VECTOR(3)=1.0
C           INTENSITY        -
      VECTOR(4)=0.
C           LINE STYLE       - SOLID
      VECTOR(5)=.00024414
C           LINE WIDTH       - 1/4096
      VECTOR(6)=.01
C           CHARACTER BOX Y  - NORMAL PRINT SIZE (100 LINES/PAGE)
      VECTOR(7)=.00714286
C           CHARACTER BOX X  - NORMAL PRINT SIZE (VECTOR(6)*5/7)

C     ESTABLISH DEVICE UNITS (MAX ADDRESSABLE UNITS)

      XDEVIC=32767
      YDEVIC=32767

C     ASSIGN INPUT PARAMETERS TO ASPE AND JUST
      ASPE=ASPECT
      JUST=JUSTIF

C CHECK FOR VALID ASPECT.  IF(ASPECT.LT.0.0) THEN CALL VBERRH(721,5),
C AND USE DEFAULT ASPECT.

      IF(ASPE.LT.0.) THEN
        CALL WMETER(721,5)
        ASPE=XDEVIC/YDEVIC
      END IF

C     ESTABLISH ASPECT RATIO

C     IF=0 SET TO DEVICE DEPENDENT ASPECT RATIO(FOR dic  ASPECT=32148/21698

      IF(ASPE .EQ. 0.) ASPE = XDEVIC/YDEVIC

      IF (ASPE .GT. 1.) THEN
          XNDCMX = 1.
          YNDCMX = 1./ASPE
      ELSE
          XNDCMX=ASPE
          YNDCMX=1.
      END IF

C     DEFINE MAPPING FUNCTION FOR ANY DEVICE

C     ESTABLISH SCALE FACTOR FOR MAXIMUM SCREEN DIMENSIONS OF THE DEVICE

      XSCALE = XDEVIC /XNDCMX
      YSCALE= YDEVIC / YNDCMX
      XSCALE = MIN( XSCALE, YSCALE)
      YSCALE = XSCALE

C CHECK FOR VALID JUSTIF.  IF(JUSTIF.LT.0 .OR. JUSTIF.GT.9) THEN
C CALL VBERRH(720,5), AND USE DEFAULT JUSTIF.

      IF(JUST .LT. 0 .OR. JUST .GT. 9) CALL WMETER(720,5)

C  MAKE OUTPUT FILE FOR THE METAFILE BE UNIT 55
      KOUTFL=55

C  SET UP MONITORING INFORMATION
      CALL WMETDV('C MET   ')
      CALL WMETMO(0)

C INITIALIZE THE OUTPUT FILE.
      CALL WMETFF(KOUTFL,1440,1,ISTAT)

C COMPUTE LENGTH OF FILE ID INSTRUCTION.
      IIDSIZ=(KIDSIZ+1)/2
      IUSRSZ=(KUSRSZ+1)/2
      ISZROU=(KSZROU+1)/2
      IDFILE(1)=33792+12+IIDSIZ+IUSRSZ+ISZROU
      IDFILE(2)=KCOMTP
c SEND FILE ID.

      CALL WMET13(2,IDFILE)

C SEND DATE AND TIME.
      CALL WMET05(3,KJDATE)
      CALL WMET05(3,KJTIME)

C SEND LENGTH OF JOB ID AND JOB ID.
      CALL WMET13S(IIDSIZ)
      CALL WMET05(4,KJOBID)

C SEND LENGTH OF USER ID AND USER ID.
      CALL WMET13S(IUSRSZ)
      CALL WMET05(4,KUSRID)

C SEND LENGTH OF ROUTING INFO AND ROUTING INFO.
      CALL WMET13S(ISZROU)
      CALL WMET05(4,KJROUT)

C SEND SECURITY AND FLUSH BUFFER.
      ISECUR(3)=KSECUR
      CALL WMET13(3,ISECUR)
      zero(1) = 0
      CALL WMET13(0,zero)

C SEND ASPECT RATIO
      CALL WMET13S(33539)
      CALL WMET13S(INT(XNDCMX*XDEVIC))
      CALL WMET13S(INT(YNDCMX*YDEVIC))
      CALL WMET13S(0)

C          SET UP COLOR TABLE

      DO 10 I=2,256
         CLRTAB(I,1) = 255
         CLRTAB(I,2) = 255
         CLRTAB(I,3) = 255
   10 CONTINUE
      CLRTAB(1,1) = 0
      CLRTAB(1,2) = 0
      CLRTAB(1,3) = 0
      CLRTAB(2,2) = 0
      CLRTAB(2,3) = 0
      CLRTAB(3,1) = 0
      CLRTAB(3,3) = 0
      CLRTAB(4,3) = 0
      CLRTAB(5,1) = 0
      CLRTAB(5,2) = 0
      CLRTAB(6,2) = 0
      CLRTAB(7,1) = 0

C SET CURRENT POSITION TO (0.,0.)
      CALL WMETIM(rzero, rzero)

      RETURN
      END
      SUBROUTINE WMETIL(XNDC,YNDC)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VILINA           -Line Absolute.

C R.W.Simons       -08MAY81

C ENVIRONMENT      -Computer-independent, system-independent, FORTRAN 77
C                   Hard Copy Format 1.

C ENTRY CONDITIONS -XNDC,YNDC = real NDC position.

C CALLS            -

C EXIT CONDITIONS  -XCP,YCP = real updated current position. (XNDC,YNDC)

C NARRATIVE        -Draw a line from current position to absolute NDC
C                   position XNDC,YNDC and update current position.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      REAL*4 XNDC,YNDC

      INTEGER*4 OUTARY(2)

C CURRENT POSITION. (LXY,HC1)
      REAL*4 XCP,YCP
      COMMON /WMET09/ XCP,YCP
C      Include '[VDIMAINT.COMMON]VCCRPS'

C SCALE FACTORS FOR NDC TO DC MAPPING. (LXY)
      REAL*4 XSCALE,YSCALE
      COMMON /WMET10/ XSCALE,YSCALE

C DC DIMENSIONS OF OFFSETS AND PICTURE. (LXY,HC1)
      REAL*4 XPAD,YPAD,XDEVIC,YDEVIC
      COMMON /WMET11/ XPAD,YPAD,XDEVIC,YDEVIC

C CONVERT TO SCREEN UNITS.
C SET BIT 15 OF X = 0 TO INDICATE A COORDINATE POSITIONING COMMAND.
C SET BIT 15 OF Y = 1 TO INDICATE A DRAW COMMAND.
C ASSUME X AND Y ARE < 32768, WHICH WILL BE TRUE IF XNDC AND YNDC
C ARE IN THE PROPER RANGE.
      OUTARY(1)=XNDC*XSCALE+XPAD
      OUTARY(2)=YNDC*YSCALE+YPAD+32768

C SEND THE COMMAND TO THE PLOT FILE.
      CALL WMET13(2,OUTARY)

C UPDATE CURRENT POSITION.
      XCP=XNDC
      YCP=YNDC

      RETURN
      END
      SUBROUTINE WMETIP(X,Y)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VIPNTA           -Point Absolute.

C R.W.Simons       -15MAY81

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   Hard Copy Format 1.

C ENTRY CONDITIONS -X,Y = real NDC position.

C CALLS            -VIMOVA,VBOUT.

C EXIT CONDITIONS  -

C NARRATIVE        -Set current position to absolute NDC position X,Y
C                   and put a visible point there.  Attributes
C                   foreground color and intensity apply.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      REAL*4 X,Y

C MOVE TO THE POSITION SPECIFIED.
      CALL WMETIM(X,Y)
C plot marker at current position (A400 HEX = 41984) TO THE OUTPUT FILE.
      CALL WMET13S(41984)

      RETURN
      END
      SUBROUTINE WMETIT
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VITERM           -Terminate SVDI.

C R.W.Simons       -13MAY81

C ENVIRONMENT      -COMPUTER-INDEPENDENT, SYSTEM-INDEPENDENT, FORTRAN 77
C                   Hard Copy Format 1.

C ENTRY CONDITIONS -

C CALLS            -VBOUT,VBERRH.

C EXIT CONDITIONS  -

C NARRATIVE        -Terminate the SVDI.  Flush buffers, etc.

C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

      integer*4 i

      INTEGER*4 KWRTFL,KRDFL,KOUTFL,KINFL,KWRDSZ,KBYTEL,KCPW,KBAUD,
     1KCOMTP
      COMMON /CDRCOM/ KWRTFL,KRDFL,KOUTFL,KINFL,KWRDSZ,KBYTEL,KCPW,
     1 KBAUD,KCOMTP
      INTEGER FILLER

C FILLER IS 1440 BYTES WORTH OF DECIMAL 8100 PADD CHARACTERS
      DATA FILLER /33024/
      integer*4 zero(1)

C SEND AN END OF DATA COMMAND TO THE OUTPUT FILE.
C 8600 HEX = 34304 = END OF DATA COMMAND.
      CALL WMET13S(34304)

C UPON TERMINATION, WE WANT TO SEND AN EXTRA BUFFER FULL OF PADD
C CHARACTERS.  BY DOING THIS, WE CAN EASILY BUILD A STANDARD
C METAFILE WHICH CAN BE PASSED AROUND SYSTEMS.  IT MAY BE
C NECESSARY FOR SOME SYSTEMS TO THROW AWAY DATA AT THE END OF
C THE FILE, AND THIS WILL ENSURE THAT NOTHING WORTHWHILE GETS
C DISCARDED.

      zero(1) = 0
      DO 10 I=1,2048
         CALL WMET13S(FILLER)
  10  CONTINUE

C FLUSH OUTPUT BUFFERS.
      CALL WMET13(0,zero)

      CALL WMETCF(KOUTFL,1)
      CALL WMETMO(2)

      RETURN
      END
      SUBROUTINE WMETIX(LENGT1,CHARS)
C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C C

C VITEXT           - Text from Array.

C P. Watterberg    - 24 MAR 81

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

      INTEGER   JSPOT
      INTEGER*4 LENGT1, CHARS(136), LENGTH, i
      REAL*4 VECTOR(7)
      COMMON /WMET04/ VECTOR

      COMMON /WMET09/XCP,YCP
      REAL*4 XCP,YCP

      LOGICAL ODD
      INTEGER*4 TBUF(73)

c          check for valid length.

      LENGTH = LENGT1
      IF(LENGTH.LT.1) THEN
         CALL WMETER(212,5)
         GO TO 999
         END IF

c          if(length.gt.136) then call vberrh(213,5), and use the
c          maximum length of 136.

      IF(LENGTH.GT.136) THEN
         CALL WMETER(213,5)
         LENGTH = 136
         ENDIF

c          initialized the number of chars in tbuf,
c          the spot marker and the odd/even flag

      ODD = .TRUE.
      JSPOT = 1

c          loop through length characters.

      DO 100 I=1,LENGTH

c          check for valid chars.
c          ignore control characters, except backspace and linefeed.

         IF(CHARS(I).LT.32 .OR. CHARS(I).GT.126) THEN
                IF(CHARS(I).NE.8.AND.CHARS(I).NE.10) THEN
                   CALL WMETER(208,5)
                   GOTO 100
                ENDIF
            END IF

c          now pack the chars into the buffer

         IF(ODD)  THEN
             JSPOT = JSPOT + 1
             TBUF(JSPOT) = CHARS(I)*256
           ELSE
             TBUF(JSPOT) = TBUF(JSPOT) + CHARS(I)
           ENDIF
         ODD = .NOT. ODD

C          UPDATE THE CURRENT POSITION

         IF(CHARS(I).GE.32) THEN
            XCP = XCP + VECTOR(7)
         ELSE IF(CHARS(I).EQ.10) THEN
            YCP = YCP - VECTOR(6)
         ELSE
            XCP = XCP - VECTOR(7)
         ENDIF

  100 CONTINUE

c          send the chars to the bgp file

c  45056 :: b000      tbuf(1)=45056  + jspot -1
c  jspot is the number of words filled or partially filled up in tbuf
c    including the 1st word with the command

      TBUF(1)=JSPOT+45055
      CALL WMET13(JSPOT,TBUF)

  999 RETURN
      END
