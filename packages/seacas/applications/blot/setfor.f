C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C============================================================================
      SUBROUTINE SETFOR (IFUNC, INLINE, IFLD, INTYP, CFIELD, *)
C============================================================================

C   --*** SETFOR ***  (BLOT) Process FOREGROUND command
C   --   Written by John Glick - 2/27/89
C   --
C   --Parameters:
C   --   IFUNC  - IN - = 1 if the call is to parse the FOREGROUND
C   --                     command.
C   --                 = 2 if the call is to reset the foreground color
C   --                     to the default color.
C   --   INLINE - IN/OUT - the parsed input lines for the log file
C   --   IFLD, INTYP, CFIELD, - IN/OUT - the free-field reader
C   --          index and character field.
C   --   * - return statement if command error; message is printed.

      include 'params.blk'
      INTEGER IFUNC
      CHARACTER*(*) INLINE(*)
      INTEGER IFLD, INTYP(*)
      CHARACTER*(*) CFIELD(*)

      include 'plcolr.blk'
      include 'plcol2.blk'

      LOGICAL FFEXST
      INTEGER IDCOL, LOCSTR

      CHARACTER*(mxstln) COLA, COLF
      include 'cmap-lst.blk'
      INTEGER DEFFOR
      SAVE DEFFOR

      DATA DEFFOR / 2 /

c      DATA FORGND /'WHITE   '/

C *****************************************************************

      IF (IFUNC .EQ. 1) THEN

         IF (FFEXST (IFLD, INTYP)) THEN

C              Check that next field has characters in it.

            IF (INTYP(IFLD) .GE. 0) THEN
               COLA = CFIELD(IFLD)
               IFLD = IFLD + 1
               CALL ABRSTR (COLF, COLA, COLLST)
               IF (COLF .EQ. ' ') THEN
                  WRITE (*, 10000) COLA
10000              FORMAT (1X, A, ' not a valid color name.')
                  GO TO 100
               ELSE
                  IDCOL = LOCSTR (COLF, NCOLOR, COLLST)
                  CALL FFADDC (COLF, INLINE)
                  IF (IDCOL .GT. 0) THEN
                     FORGND = COLF
                     IDFOR = IDCOL
                     IDFORT = IDCOL
                  ELSE
                     WRITE (*, 10000) COLA
                     GO TO 100
                  ENDIF
               ENDIF
            ELSE
               CALL PRTERR ('CMDERR',
     &            'Expected color name following FOREGROUND command')
               GOTO 100
            ENDIF
         ELSE

            CALL PRTERR ('CMDERR',
     &         'Expected color name following FOREGROUND command')
            GOTO 100

         ENDIF

      ELSE IF (IFUNC .EQ. 2) THEN

         IDFOR = DEFFOR
         IDFORT = DEFFOR
         FORGND = 'WHITE   '

      ENDIF

      RETURN

  100 CONTINUE
      RETURN 1

      END
