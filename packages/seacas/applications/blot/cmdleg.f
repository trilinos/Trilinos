C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CMDLEG (INLINE, VERB, IFLD, INTYP, CFIELD, IFIELD,
     &           XYTYPE, MESHOK, DOQA, DOLEG, DOAXIS, DOBOX, CAPTN, *)
C=======================================================================

C   --*** CMDLEG *** (BLOT) Process legend labeling command
C   --   Written by Amy Gilkey - revised 04/21/88
C   --
C   --CMDLEG processes a legend labeling command.  The commands are:
C   --   QA      - controls whether the QA information is to be drawn
C   --   LEGEND  - controls whether non-QA legend information is to be drawn
C   --   AXIS    - controls whether the axis is to be numbered
C   --   CAPTION - sets the plot caption
C   --
C   --Parameters:
C   --   INLINE - I/O - the parsed input lines for the log file
C   --   VERB   - I/O - the command verb; set for SHOW
C   --   IFLD   - I/O - the field number
C   --   INTYP  - IN  - the input types from the free field reader
C   --   CFIELD - IN  - the character fields
C   --   IFIELD - IN  - the integer fields
C   --   XYTYPE - IN  - true iff curr program is an XY curve versus mesh plot
C   --   MESHOK - IN  - true iff mesh can be displayed
C   --   DOQA   - I/O - true iff QA information is to be drawn on legend
C   --                  (1) for mesh plots, (2) for curve plots
C   --   DOLEG  - I/O - true iff non-QA legend information is to be drawn
C   --                  (1) for mesh plots, (2) for curve plots
C   --   DOAXIS - I/O - true iff axis is to be labeled
C   --                  (1) for mesh plots, (2) for curve plots
C   --   CAPTN  - I/O - the three-line plot caption
C   --                  (1) for mesh plots, (2) for curve plots
C   --   *      - OUT - return statement if command error; message is printed

      include 'params.blk'
      CHARACTER*(*) INLINE(*)
      CHARACTER*(*) VERB
      INTEGER     INTYP(*)
      CHARACTER*(*) CFIELD(*)
      INTEGER     IFIELD(*)
      LOGICAL XYTYPE, MESHOK
      LOGICAL DOQA(2), DOLEG(2), DOAXIS(2), DOBOX
      CHARACTER*(*) CAPTN(3,2)

      INTEGER IDUM
      REAL RDUM
      CHARACTER*(MXSTLN) CDUM

      CHARACTER*(MXSTLN) WORD
      LOGICAL ISON

      DATA IDUM / 1 /

      IF (VERB .EQ. 'QA') THEN
         CALL FFADDC (VERB, INLINE(1))
         CALL FFONOF (IFLD, INTYP, CFIELD, ISON, *110)
         CALL FFADDO (ISON, INLINE(1))

      ELSE IF (VERB .EQ. 'LEGEND') THEN
         CALL FFADDC (VERB, INLINE(1))
         CALL FFONOF (IFLD, INTYP, CFIELD, ISON, *110)
         CALL FFADDO (ISON, INLINE(1))

      ELSE IF (VERB .EQ. 'AXIS') THEN
         CALL FFADDC (VERB, INLINE(1))
         CALL FFONOF (IFLD, INTYP, CFIELD, ISON, *110)
         CALL FFADDO (ISON, INLINE(1))

      ELSE IF (VERB .EQ. 'OUTLINE') THEN
         CALL FFADDC (VERB, INLINE(1))
         CALL FFONOF (IFLD, INTYP, CFIELD, ISON, *110)
         CALL FFADDO (ISON, INLINE(1))

      ELSE IF (VERB .EQ. 'CAPTION') THEN
         CALL FFADDC (VERB, INLINE(1))
         CALL FFINTG (IFLD, INTYP, IFIELD,
     &      'caption line number (1..3) or 0', 0, ISTART, *110)
         CALL FFADDI (ISTART, INLINE(1))
         IF ((ISTART .LT. 0) .OR. (ISTART .GT. 3)) THEN
            CALL PRTERR ('CMDERR',
     &         'Expected caption line number (1..3)')
            GOTO 110
         END IF
         IF (ISTART .EQ. 0) THEN
            ISTART = 1
            IEND = 3
         ELSE
            IEND = ISTART
         END IF

      ELSE
         GOTO 110
      END IF

      CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
      IF (WORD .EQ. ' ') THEN
         ITYPE = 1
         IF (XYTYPE) ITYPE = 2
      ELSE IF (WORD .EQ. 'MESH') THEN
         IF (.NOT. MESHOK) THEN
            CALL PRTERR ('CMDERR', 'Mesh is not defined')
            GOTO 110
         END IF
         ITYPE = 1
      ELSE IF (WORD .EQ. 'XY') THEN
         IF (.NOT. XYTYPE) THEN
            CALL PRTERR ('CMDERR',
     &         'XY option not allowed from mesh subprogram')
            GOTO 110
         END IF
         ITYPE = 2
      ELSE
         CALL PRTERR ('CMDERR', 'Expected "MESH" or "XY"')
         GOTO 110
      END IF
      IF (WORD .NE. ' ') CALL FFADDC (WORD, INLINE(1))

      IF (VERB .EQ. 'QA') THEN
         DOQA(ITYPE) = ISON

      ELSE IF (VERB .EQ. 'LEGEND') THEN
         DOLEG(ITYPE) = ISON

      ELSE IF (VERB .EQ. 'AXIS') THEN
         DOAXIS(ITYPE) = ISON

      ELSE IF (VERB .EQ. 'OUTLINE') THEN
         DOBOX = ISON

      ELSE IF (VERB .EQ. 'CAPTION') THEN
         DO 100 I = ISTART, IEND
            WRITE (WORD, '(''LINE '',(I1),''> '')') I
            LWORD = LENSTR (WORD)
            CALL GETINS ('line', IDUM, IDUM, IDUM, CDUM,
     &         IDUM, RDUM, CAPTN(I,ITYPE), IOSTAT,
     &         WORD, LWORD, *110)
            INLINE(I-ISTART+2) = CAPTN(I,ITYPE)
  100    CONTINUE
      END IF

      RETURN

  110 CONTINUE
      RETURN 1
      END
