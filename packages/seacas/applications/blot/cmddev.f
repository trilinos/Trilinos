C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CMDDEV (INLINE,
     &   VERB, IFLD, INTYP, CFIELD, IFIELD, RFIELD, *)
C=======================================================================

C   --*** CMDDEV *** (BLOT) Process device parameter command
C   --   Written by Amy Gilkey - revised 04/26/88
C   --
C   --CMDDEV processes a device parameter command. The commands are:
C   --   SOFTCHAR - sets the software vs. hardware characters flag
C   --   COLOR - sets the number of standard colors to use
C   --   SPECTRUM - sets the number of spectrum colors
C   --   FONT - sets font to use
C   --   SNAP - sets the number of frames to snap
C   --   AUTO - sets the device for automatic plotting
C   --
C   --Parameters:
C   --   INLINE - IN/OUT - the parsed input lines for the log file
C   --   VERB - IN/OUT - the command verb; set for SHOW
C   --   IFLD - IN/OUT - the field number
C   --   INTYP - IN - the input types from the free field reader
C   --   CFIELD - IN - the character fields
C   --   IFIELD - IN - the integer fields
C   --   RFIELD - IN - the real fields
C   --   * - return statement if command error; message is printed

      include 'params.blk'
      include 'colormap.blk'
      CHARACTER*(*) INLINE(*)
      CHARACTER*(*) VERB
      INTEGER     INTYP(*)
      CHARACTER*(*) CFIELD(*)
      INTEGER     IFIELD(*)
      REAL        RFIELD(*)
      include 'icrnbw.blk'

      LOGICAL MATSTR
      CHARACTER*80 ERRMSG
      CHARACTER*(MXSTLN) WORD, WORDM
      INTEGER ISON
      include 'shades.blk'

      IF (VERB .EQ. 'SOFTCHAR') THEN
         CALL FFADDC (VERB, INLINE(1))
         CALL FFONOF (IFLD, INTYP, CFIELD, ISON, *110)
         CALL FFADDO (ISON, INLINE(1))
         CALL FFINTG (IFLD, INTYP, IFIELD,
     &      'device number', 0, IDEV, *110)
         IF (IDEV .NE. 0) CALL FFADDI (IDEV, INLINE(1))
         CALL GRSPAR (VERB, IDEV, ISON, ERRMSG)
         IF (ERRMSG .NE. ' ') GOTO 100

      ELSE IF (VERB .EQ. 'FONT') THEN
         CALL FFADDC (VERB, INLINE(1))
         CALL FFCHAR (IFLD, INTYP, CFIELD, 'STICK', WORD)
         IF (MATSTR (WORD, 'STICK', 2)) THEN
            CALL FFADDC ('STICK', INLINE(1))
            IFONT = 1
         ELSE IF (MATSTR (WORD, 'SANSERIF', 2)) THEN
            CALL FFADDC ('SANSERIF', INLINE(1))
            IFONT = 2
         ELSE IF (MATSTR (WORD, 'ROMAN', 1)) THEN
            CALL FFADDC ('ROMAN', INLINE(1))
            IFONT = 3
         ELSE
            CALL PRTERR ('CMDERR',
     &         'Expected "STICK", "SANSERIF" or "ROMAN"')
            GOTO 110
         END IF
         CALL FFINTG (IFLD, INTYP, IFIELD,
     &      'device number', 0, IDEV, *110)
         IF (IDEV .NE. 0) CALL FFADDI (IDEV, INLINE(1))
         CALL GRSPAR (VERB, IDEV, IFONT, ERRMSG)
         IF (ERRMSG .NE. ' ') GOTO 100

      ELSE IF (VERB .EQ. 'COLOR') THEN
         CALL FFADDC (VERB, INLINE(1))
         CALL FFINTG (IFLD, INTYP, IFIELD,
     &      'number of colors', 0, NCOL, *110)
         CALL FFADDI (NCOL, INLINE(1))
         CALL GRSPAR (VERB, -1, NCOL, ERRMSG)
         IF (ERRMSG .NE. ' ') GOTO 100

      ELSE IF (VERB .EQ. 'SPECTRUM') THEN
         CALL FFADDC (VERB, INLINE(1))
         CALL FFINTG (IFLD, INTYP, IFIELD,
     &      'number of spectrum colors', 5, NCOL, *110)
         CALL FFADDI (NCOL, INLINE(1))
         ISINV = 0
 10      continue
         CALL FFCHAR (IFLD, INTYP, CFIELD, 'DEFAULT', WORD)
         CALL FFADDC (WORD, INLINE(1))
         IF (MATSTR (WORD, 'DEFAULT', 3)) THEN
           ISPEC = DEFAULT
         ELSE IF (MATSTR (WORD, 'RAINBOW', 3)) THEN
           ISPEC = RAINBW
         ELSE IF (MATSTR (WORD, 'VIRIDIS', 3)) THEN
           ISPEC = VIRDIS
         ELSE IF (MATSTR (WORD, 'GRAY', 2) .or.
     *            MATSTR (WORD, 'GREY', 2)) THEN
           ISPEC = GRAY
         ELSE IF (MATSTR (WORD, 'TERRAIN', 3) .or.
     *            MATSTR (WORD, 'TOPOGRA', 3)) THEN
           ISPEC = TERRAIN
         ELSE IF (MATSTR (WORD, 'HOT', 3) .or.
     *            MATSTR (WORD, 'HEATED', 3) .or.
     *            MATSTR (WORD, 'THERMAL', 3)) THEN
           ISPEC = IRON
         ELSE IF (MATSTR (WORD, 'ASTRO', 3)) THEN
           ISPEC = ASTRO
         ELSE IF (MATSTR (WORD, 'ZEBRA', 3)) THEN
           ISPEC = ZEBRA
         ELSE IF (MATSTR (WORD, 'COOL', 3)) THEN
           ISPEC = COOL
         ELSE IF (MATSTR (WORD, 'METAL', 3) .or.
     *       MATSTR (WORD, 'USER',  3)) THEN
           ISPEC = METAL
           IF (INTYP(IFLD) .EQ. 0) then
C ... User has specified a 'predefined' color.
             CALL FFCHAR (IFLD, INTYP, CFIELD, 'DEFAULT', WORD)
             CALL FFADDC (WORD, INLINE(1))
             CALL ABRSTR (WORDM, WORD, SHDLST)
             IF (WORDM .EQ. ' ') THEN
               WRITE (*, 10000) WORD
10000          FORMAT (1X, A, ' not a valid color name.')
               CALL SHOCMD ('Valid predefined colors', SHDLST)
               RMULT = 1.000
               GMULT = 1.000
               BMULT = 1.000
             ELSE
               IDCOL = LOCSTR(WORDM, NCLSHD, SHDLST)
               RMULT = shades(1,IDCOL)
               GMULT = shades(2,IDCOL)
               BMULT = shades(3,IDCOL)
             END IF
           else
C ... User has specified the RGB components of the color.
             CALL FFREAL (IFLD, INTYP, RFIELD,
     *         'Red Multiplier', 1.0, RMULT, *110)
             CALL FFADDR (RMULT, INLINE(1))
             CALL FFREAL (IFLD, INTYP, RFIELD,
     *         'Green Multiplier', 1.0, GMULT, *110)
             CALL FFADDR (GMULT, INLINE(1))
             CALL FFREAL (IFLD, INTYP, RFIELD,
     *         'Blue Multiplier', 1.0, BMULT, *110)
             CALL FFADDR (BMULT, INLINE(1))
           END IF
         ELSE IF (MATSTR (WORD, 'INVERSE', 3)) THEN
           ISINV = 1
           GO TO 10
         ELSE IF (MATSTR (WORD, 'HELP', 4)) THEN
           CALL PRTERR ('CMDSPEC',
     *       'Valid: RAINBOW, VIRIDIS, GRAY, TERRAIN, THERMAL, ASTRO,')
           CALL PRTERR ('CMDSPEC',
     *       '       COOL, METAL, a color, or enter rgb triplet')
           ERRMSG = 'Unknown Color Map'
           GO TO 100
         ELSE
           ERRMSG = 'Unknown Color Map'
           GO TO 100
         END IF
         CALL FFREAL (IFLD, INTYP, RFIELD,
     *     'SATURATION', 1.0, SATUR, *110)
         CALL FFADDR (SATUR, INLINE(1))
         CALL GRSPAR (VERB, -1, NCOL, ERRMSG)
         IF (ERRMSG .NE. ' ') GOTO 100

      ELSE IF (VERB .EQ. 'SNAP') THEN
         CALL FFADDC (VERB, INLINE(1))
         CALL FFINTG (IFLD, INTYP, IFIELD,
     &      'number of frames to snap', 1, N, *110)
         CALL FFADDI (N, INLINE(1))
         IDEV = 0
         CALL GRSPAR (VERB, IDEV, N, ERRMSG)
         IF (ERRMSG .NE. ' ') GOTO 100

      ELSE IF (VERB .EQ. 'AUTO') THEN
         CALL FFADDC (VERB, INLINE(1))
         CALL FFONOF (IFLD, INTYP, CFIELD, ISON, *110)
         CALL FFADDO (ISON, INLINE(1))
         IDEV = 0
         CALL GRSPAR (VERB, IDEV, ISON, ERRMSG)
         IF (ERRMSG .NE. ' ') GOTO 100
      END IF

      RETURN

  100 CONTINUE
      CALL PRTERR ('CMDERR', ERRMSG(:LENSTR(ERRMSG)))
  110 CONTINUE
      RETURN 1
      END
