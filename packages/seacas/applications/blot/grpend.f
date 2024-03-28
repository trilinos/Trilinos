C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE GRPEND (MAYQUI, MAYTXT, NDONE, NTOTAL, GOBCK, *, *)
C=======================================================================

C   --*** GRPEND *** (GRPLIB) End a plot by requesting quit, hardcopy (PLT)
C   --   Written by Amy Gilkey - revised 01/22/88
C   --
C   --GRPEND is called at the end of a plot.  It may prompt the user for
C   --a response and take the appropriate alternate returns for special
C   --responses.
C   --
C   --If a single hardcopy plot is requested, the hardcopy device is
C   --selected and an internal flag is set and an alternate return is taken.
C   --If the internal flag is set upon entry to the routine, the terminal
C   --device is selected and the flag is reset.  Be sure to call this
C   --routine at the normal end of the each plot (even a single hardcopy plot).
C   --The internal flag is reset whenever the terminal device is selected.
C   --
C   --If the plot may be annotated with text, be sure that the color of the
C   --text to be displayed is set before this routine is called.
C   --
C   --If no response is possible and AUTO is set, the number of plots
C   --completed is displayed.
C   --
C   --Parameters:
C   --   MAYQUI - IN - true if QUIT is a possible alternative
C   --   MAYTXT - IN - true if TEXT is a possible alternative
C   --   NDONE - IN - the number of plots completed
C   --   NTOTAL - IN - the total number of plots
C   --   GOBCK  - in/out - on entry, true if go back supported,
C   --                   - on output, true if go back selected
C   --   * - the alternate return if hardcopy requested
C   --   * - the alternate return if quit requested
C   --
C   --Common Variables:
C   --   Sets and uses IHARD of /GRPCOM/
C   --   Uses ICURDV, DEVOK, DEVCOD, AUTOPL of /GRPCOM/

C   --Routines Called:
C   --   FREFLD - (SUPES) Free-field reader
C   --   SQZSTR - (STRLIB) Compress extra blanks
C   --   PLTBEL - (PLT) Ring bell
C   --   PLTFLU - (PLT) Flush buffer
C   --   PLTMOV - (PLT) Move cursor
C   --   PLTWAI - (PLT) Wait until key is pressed
C   --   GRABRT - (GRPLIB) Check for plot set abort
C   --   GRIKEY - (GRPLIB) Select cursor position
C   --   GRSDEV - (GRPLIB) Select device
C   --   GRSNAP - (GRPLIB) Handle device frame snapping

      PARAMETER (MXTEXT = 80)

      include 'grpcom.blk'

      LOGICAL MAYQUI, MAYTXT, GOBCK

      LOGICAL GRABRT
      LOGICAL REQWAI, ASKQUI, ASKHRD, ASKTXT, REPEAT
      CHARACTER*80 PROMPT
      CHARACTER REPLY

      LOGICAL LSTERR

      INTEGER NTEXT
      LOGICAL CTEXT(MXTEXT)
      REAL XTEXT(MXTEXT), YTEXT(MXTEXT)
      CHARACTER*80 TEXT(MXTEXT)
      SAVE NTEXT, CTEXT, XTEXT, YTEXT, TEXT

C   --If end of hardcopy plot, add any text

      IF (ISHARD) THEN
         DO 100 I = 1, NTEXT
            IF (CTEXT(I)) THEN
               CALL GRTEXC (XTEXT(I), YTEXT(I), TEXT(I))
            ELSE
               CALL GRTEXT (XTEXT(I), YTEXT(I), TEXT(I))
            END IF
  100    CONTINUE
      END IF

C   --End plot and snap frames

      CALL GRSNAP ('STOP', 0)

C   --Flush the buffer and ring bell and (on some devices) wait for key

      IF (GRABRT ()) RETURN 2
      CALL PLTMOV (0.0, 0.0)
      CALL PLTFLU
      CALL PLTBEL
      CALL PLTFLU
      REQWAI = (.NOT. AUTOPL(ICURDV)) .AND. (DEVCOD(ICURDV) .EQ. 'WAIT')
      IF (REQWAI) THEN
         CALL PLTWAI
         CALL PLTFLU
      END IF

C   --If end of hardcopy plot, reselect terminal device and return

      IF (ISHARD) THEN
         CALL GRSDEV (0)
         ISHARD = .FALSE.
         RETURN
      END IF

C   --Reset the text lines
      NTEXT = 0
      DX = 1.0 / 2
      DY = 0.75 / 2

C.. If a single plot, or the last plot in the set, don't prompt,
C   return immediately to prompt.
      ASKQUI = NDONE .LT. NTOTAL .AND. MAYQUI

      ASKHRD = DEVOK(2) .AND. (ICURDV .NE. 2)
C.. Only ask for text if the hardcopy or multiple plot prompting enabled
      ASKTXT = (ASKHRD .OR. ASKQUI) .AND. MAYTXT

      IF ((.NOT. AUTOPL(ICURDV))
     &   .AND. (ASKQUI .OR. ASKHRD .OR. ASKTXT)) THEN

C      --Check if user wants to quit or get hardcopy

         LSTERR = .FALSE.

         PROMPT = 'Enter'
         LPR = LENSTR (PROMPT) + 1
         IF (ASKQUI) THEN
            PROMPT(LPR+1:) = '"Q" to quit,'
            LPR = LENSTR (PROMPT) + 1
         END IF
         IF (gobck) THEN
            gobck = .false.
            if (ndone .gt. 1) then
               PROMPT(LPR+1:) = '"P" for previous,'
               LPR = LENSTR (PROMPT) + 1
            end if
         END IF
         IF (ASKHRD) THEN
            PROMPT(LPR+1:) = '"H" for hardcopy,'
            LPR = LENSTR (PROMPT) + 1
         END IF
         IF (ASKTXT) THEN
            PROMPT(LPR+1:) = '"T" for text,'
            LPR = LENSTR (PROMPT) + 1
         END IF
         PROMPT(LPR+1:) = '" " to continue'
         LPR = LENSTR (PROMPT) + 1

  110    CONTINUE

C      --Prompt for response, print trailing bell only if last response
C      --was in error

         IF (LSTERR) THEN
            PROMPT(LPR+1:LPR+1) = CHAR(7)
            LPR = LPR + 1
         END IF
         CALL FREFLD (0, 0, PROMPT(:LPR), 1,
     &      IOSTAT, NUMFLD, INTYP, REPLY, IDUM, RDUM)
         IF (LSTERR) LPR = LPR - 1

         LSTERR = .FALSE.
         REPEAT = .FALSE.

         IF (REPLY .EQ. 'Q') THEN
            IF (ASKQUI) THEN
C            --Quit
               RETURN 2
            END IF

         ELSE IF (REPLY .EQ. 'H') THEN
            IF (ASKHRD) THEN
C            --Request hardcopy, select hardcopy, set flag, and return
               CALL GRSDEV (2)
               ISHARD = .TRUE.
               RETURN 1
            ELSE
               REPEAT = .TRUE.
            END IF

         ELSE IF (REPLY .EQ. 'P') THEN
            GOBCK = .TRUE.
            RETURN

         ELSE IF (REPLY .EQ. 'T') THEN
            IF (ASKTXT .AND. (NTEXT .LT. MXTEXT)) THEN

C            --Input position and text string

               CALL GRIKEY (
     &            '   Select text position ("C" to center) ...',
     &            DX, DY, REPLY, *120)

C            --Write line to correct problem with cursor input followed
C            --by GETINP or FREFLD
               WRITE (*, *)

               IF (REPLY .EQ. 'C') THEN
                  CALL GETINP (0, 0, '   Text to center> ',
     &               TEXT(NTEXT+1), IOSTAT)
               ELSE
                  CALL GETINP (0, 0, '   Text> ',
     &               TEXT(NTEXT+1), IOSTAT)
               END IF
               WRITE (*, *)

               IF (TEXT(NTEXT+1) .NE. ' ') THEN
                  NTEXT = NTEXT + 1
                  CTEXT(NTEXT) = (REPLY .EQ. 'C')
                  XTEXT(NTEXT) = DX
                  YTEXT(NTEXT) = DY

C               --Display text string on device

                  IF (CTEXT(NTEXT)) THEN
                     CALL GRTEXC (XTEXT(NTEXT), YTEXT(NTEXT),
     &                  TEXT(NTEXT))
                  ELSE
                     CALL GRTEXT (XTEXT(NTEXT), YTEXT(NTEXT),
     &                  TEXT(NTEXT))
                  END IF
                  CALL PLTFLU
               END IF
            END IF

            REPEAT = .TRUE.

         ELSE IF (REPLY .NE. ' ') THEN
            LSTERR = .TRUE.
            REPEAT = .TRUE.
         END IF

         IF (REPEAT) GOTO 110
      ELSE
         if (gobck) gobck = .false.
      END IF

      IF (AUTOPL(ICURDV) .AND. (NTOTAL .GT. 0)) THEN
         WRITE (PROMPT, 10000, IOSTAT=IDUM) NDONE, NTOTAL
10000     FORMAT (4X, 'Plot ', I5, ' of ', I5)
         CALL SQZSTR (PROMPT, LSTR)
         WRITE (*, '(1X, A)') PROMPT(:LSTR)
      END IF

      RETURN
 120  CONTINUE
      return 2
      END
