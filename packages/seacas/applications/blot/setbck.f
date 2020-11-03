C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C============================================================================
      SUBROUTINE SETBCK (IFUNC, INLINE, IFLD, INTYP, CFIELD, *)
C============================================================================

C   --*** SETBCK ***  (BLOT) Process BACKGROUND command
C   --   Written by John Glick - 2/27/89
C   --
C   --Parameters:
C   --   IFUNC  - IN - = 1 if the call is to parse the BACKGROUND
C   --                     command.
C   --                 = 2 if the call is to set the background color
C   --                     with a call to a PLT routine.
C   --                 = 3 if the call is to reset the background color
C   --                     to the default color.
C   --   INLINE - IN/OUT - the parsed input lines for the log file
C   --   IFLD, INTYP, CFIELD, - IN/OUT - the free-field reader
C   --          index and character field.
C   --   * - return statement if command error; message is printed.

      INTEGER IFUNC
      CHARACTER*(*) INLINE(*)
      INTEGER IFLD, INTYP(*)
      CHARACTER*(*) CFIELD(*)

      include 'params.blk'
      include 'grpcom.blk'
      include 'plcolr.blk'
      include 'plcol2.blk'
      include 'grcol.blk'

      LOGICAL FFEXST, PLTICL, LDUM
      INTEGER IDCOL, LOCSTR

      CHARACTER*(MXSTLN) COLA, COLF
      include 'cmap-lst.blk'

      INTEGER LSTDV, DEFBCK, LSTBCK
      SAVE LSTDV, DEFBCK, LSTBCK

      DATA LSTDV, DEFBCK, LSTBCK / 0, 1, 0 /
      BCKGND = 'BLACK                           '

C *****************************************************************

      IF (IFUNC .EQ. 2) THEN

        IF (ICURDV .NE. LSTDV) THEN
C              Get mapping of colors
          DO 100 I = 1, NCOLOR-2
            LDUM = PLTICL('WHITE', RCOLOR)
            IWHITE = RCOLOR
            LDUM = PLTICL('BLACK', RCOLOR)
            IBLACK = RCOLOR
            IF (PLTICL (COLLST(I+2), RCOLOR)) THEN
              COLMAP(I) = RCOLOR
            ELSE
              COLMAP(I) = -1.0
            ENDIF
 100      CONTINUE
        ENDIF

        IF ((ICURDV .NE. LSTDV) .OR. (LSTBCK .NE. IDBCK)) THEN
          IF (IDBCK .EQ. 1) THEN
            RCOLOR = IBLACK
          ELSE IF (IDBCK .EQ. 2) THEN
            RCOLOR = IWHITE
          ELSE IF (COLMAP(IDBCK-2) .GE. 0.0) THEN
            RCOLOR = COLMAP(IDBCK-2)
          ELSE
            RCOLOR = IBLACK
            IDBCK = DEFBCK
          ENDIF

          IF (RCOLOR .GE. 0.0) THEN
            if (icurdv .eq. 2) then
              CALL PLTSTD (2, RCOLOR+0.01)
            else
              CALL PLTSTD (2, RCOLOR)
            endif
            IDBCKT = IDBCK
          ELSE
             BLACK = IBLACK
            CALL PLTSTD (2, BLACK)
            IDBCKT = DEFBCK
          ENDIF
        ENDIF

        IF (LSTBCK .NE. IDBCK) THEN
          LSTBCK = IDBCK
        ENDIF

        IF (ICURDV .NE. LSTDV) THEN
          LSTDV = ICURDV
        ENDIF

      ELSE IF (IFUNC .EQ. 1) THEN

        IF (FFEXST (IFLD, INTYP)) THEN

C              Check that next field has characters in it.

          IF (INTYP(IFLD) .GE. 0) THEN
            COLA = CFIELD(IFLD)
            IFLD = IFLD + 1
            CALL ABRSTR (COLF, COLA, COLLST)
            IF (COLF .EQ. ' ') THEN
              WRITE (*, 10000) COLA
10000         FORMAT (1X, A, ' not a valid color name.')
              GO TO 110
            ELSE
              IDCOL = LOCSTR (COLF, NCOLOR, COLLST)
              CALL FFADDC (COLF, INLINE)
              IF (IDCOL .GT. 0) THEN
                BCKGND = COLF
                IDBCK = IDCOL
                IDBCKT = IDCOL
              ELSE
                WRITE (*, 10000) COLA
                GO TO 110
              ENDIF
            ENDIF
          ELSE
            CALL PRTERR ('CMDERR',
     &        'Expected color name following BACKGROUND command')
            GOTO 110
          ENDIF
        ELSE

          CALL PRTERR ('CMDERR',
     &      'Expected color name following BACKGROUND command')
          GOTO 110

        ENDIF

      ELSE IF (IFUNC .EQ. 3) THEN

        IDBCK = DEFBCK
        IDBCKT = DEFBCK
        BCKGND = 'BLACK   '

      ENDIF

      RETURN

 110  CONTINUE
      RETURN 1

      END
