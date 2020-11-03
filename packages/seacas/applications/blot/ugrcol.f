C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE UGRCOL (INDX, BLKCOL)
C=======================================================================

C   --*** UGRCOL *** (BLOT) Set color (PLT)
C   --   Written by John Glick - 11/29/88
C   --
C   --UGRCOL is similar to the GRCOLR routine of GRPLIB in that it
C   --sets the color of lines depending on the passed index.
C   --The colors are chosen differently than GRCOLR, however.
C   --If appropriate, the BLKCOL array is used to select the color.
C   --
C   --
C   --Parameters:
C   --   INDX - IN - the color index or block id:
C   --      If = -1, black is chosen as the color.
C   --      If = 0, the foreground color is chosen as the color.
C   --              If the foreground color is not available on
C   --              the current device, white is chosen.
C   --      If > 0, it is assumed to be an element block identifier.
C   --              If BLKCOL(INDX) > -2, then BLKCOL(INDX) is the
C   --              color index chosen (if it is a supported color on
C   --              the current device).
C   --   BLKCOL - IN - the user selected colors of the element blocks.
C   --                    BLKCOL(0) = 1 if the user defined material
C   --                                colors should be used in mesh plots.
C   --                              = -1 if program selected colors should
C   --                                not be used.
C   --                    BLKCOL(i) = the user selected color of element
C   --                               block i:
C   --                                  -2 - no color selected by user.
C   --                                  -1 - black
C   --                                   0 - white
C   --                                   1 - red
C   --                                   2 - green
C   --                                   3 - yellow
C   --                                   4 - blue
C   --                                   5 - cyan
C   --                                   6 - magenta
C   --
C   --Common Variables:
C   --   Uses ICURDV, MAXCOL, NUMCOL, MAPUSE of /GRPCOM/

C   --Routines Called:
C   --   PLTICL - (PLTLIB) Get color index from name
C   --   PLTSTD - (PLTLIB) Set device parameter
C   --      1 = (KCOLOR) set color
C   --   PLTSTG - (PLTLIB) Set graph parameter
C   --      6 = (KCOLIN) set line color for PLTGRH lines
C   --      44 = (KCOSYM) set symbol color for PLTGRH lines

      PARAMETER (KCOLOR=1)
      PARAMETER (KCOLIN=6, KCOSYM=44)

      COMMON /GRPCOC/ DEVNAM(2), DEVCOD(2)
      CHARACTER*3 DEVNAM
      CHARACTER*8 DEVCOD
      COMMON /GRPCOM/ ICURDV, ISHARD, DEVOK(2), TALKOK(2),
     &   NSNAP(2), IFONT(2), SOFTCH(2), AUTOPL(2),
     &   MAXCOL(2), NUMCOL(0:1,2), MAPALT(2), MAPUSE(2)
      LOGICAL ISHARD, DEVOK, TALKOK, SOFTCH, AUTOPL
      include 'dbnums.blk'
      include 'grcol.blk'

      COMMON /BCOLR/ BCOLCH
      LOGICAL BCOLCH

      include 'plcol2.blk'
      include 'params.blk'
      include 'cmap-lst.blk'
      INTEGER LSTD2
      SAVE LSTD2

      INTEGER INDX
      INTEGER BLKCOL(0:NELBLK)

      LOGICAL LDUM, PLTSTG, PLTSTD, PLTICL

      LOGICAL COLSPC(NCOLOR)
      INTEGER NUMFRE, FRECOL(NCOLOR)
      SAVE COLSPC, NUMFRE, FRECOL

      CHARACTER*8 BLCOLR, WHCOLR
      SAVE BLCOLR, WHCOLR

      DATA LSTD2 / 0 /
      DATA BLCOLR, WHCOLR / 'BLACK   ', 'WHITE   ' /

      IF (ICURDV .NE. LSTDEV) THEN
         LSTDEV = ICURDV
         LSTMAP = 0

C      --Reset last color
         LSTCOL = -999

C      --Set device black and white
         IF (MAXCOL(ICURDV) .GT. 0) THEN
            LDUM = PLTICL (WHCOLR, RWHITE)
            IWHITE = INT(RWHITE)
            LDUM = PLTICL (BLCOLR, RBLACK)
            IBLACK = INT(RBLACK)
         ELSE
            IBLACK = 0
            IWHITE = 1
         END IF

      ELSE IF (MAPUSE(ICURDV) .NE. LSTMAP) THEN

C      --Reset last color
         LSTCOL = -999
         LSTMAP = MAPUSE(ICURDV)

      ENDIF

      IF (ICURDV .NE. LSTD2) THEN
C           Get mapping of colors
         DO 100 I = 1, NCOLOR
            IF (PLTICL (COLLST(I+2), RCOLOR)) THEN
               COLMAP(I) = RCOLOR
            ELSE
               COLMAP(I) = -1
            ENDIF
  100    CONTINUE
      ENDIF

      IF (BCOLCH  .OR.  (ICURDV .NE. LSTD2)) THEN
         BCOLCH = .FALSE.
         DO 110 I = 1, NCOLOR
            COLSPC(I) = .FALSE.
  110    CONTINUE
         IF (BLKCOL(0) .EQ. 1) THEN
            DO 120 I = 1, NELBLK
               IF ((BLKCOL(I) .GT. 0)  .AND.
     &            (COLMAP(BLKCOL(I)) .NE. -1))
     &            COLSPC(BLKCOL(I)) = .TRUE.
  120       CONTINUE
         ENDIF
         NUMFRE = 0
         DO 130 I = 1, NCOLOR
            IF ((.NOT. COLSPC(I))  .AND.  (COLMAP(I) .NE. -1)) THEN
               NUMFRE = NUMFRE + 1
               FRECOL(NUMFRE) = I
            ENDIF
  130    CONTINUE
      ENDIF

      IF (ICURDV .NE. LSTD2) LSTD2 = ICURDV

      NCOL = NUMCOL(MAPUSE(ICURDV),ICURDV)

      IF (BLKCOL(0) .EQ. 1) THEN
         IF (INDX .EQ. -1) THEN
            ICOLOR = IBLACK
         ELSE IF (INDX .EQ. 0) THEN
            IDFORT = IDFOR
            IF (IDFOR .EQ. 1) THEN
               ICOLOR = IBLACK
            ELSE IF (IDFOR .EQ. 2) THEN
               ICOLOR = IWHITE
            ELSE IF (COLMAP(IDFOR-2) .GE. 0.0) THEN
               ICOLOR = INT(COLMAP(IDFOR-2))
            ELSE
               ICOLOR = IWHITE
               IDFORT = 2
            ENDIF
            COLFOR = REAL (ICOLOR)
         ELSE IF ((BLKCOL(INDX) .GT. 0)  .AND.
     &      (COLMAP(BLKCOL(INDX)) .NE. -1)) THEN
            ICOLOR = INT(COLMAP(BLKCOL(INDX)))
         ELSE IF (BLKCOL(INDX) .EQ. 0) THEN
            ICOLOR = IWHITE
         ELSE IF (BLKCOL(INDX) .EQ. -1) THEN
            ICOLOR = IBLACK
         ELSE
            IF (NUMFRE .GT. 0) THEN
               NBLNSP = 0
               DO 140 I = 1, INDX
                  IF (BLKCOL(I) .EQ. -2) THEN
                     NBLNSP = NBLNSP + 1
                  ELSE IF (BLKCOL(I) .GT. 0) THEN
                     IF (COLMAP(BLKCOL(I)) .EQ. -1) THEN
                        NBLNSP = NBLNSP + 1
                     ENDIF
                  ENDIF
  140          CONTINUE
               ICOLOR = INT(COLMAP (FRECOL(MOD (NBLNSP-1, NUMFRE)+1)))
            ELSE
               IDFORT = IDFOR
               IF (IDFOR .EQ. 1) THEN
                  ICOLOR = IBLACK
               ELSE IF (IDFOR .EQ. 2) THEN
                  ICOLOR = IWHITE
               ELSE IF (COLMAP(IDFOR-2) .GE. 0.0) THEN
                  ICOLOR = INT(COLMAP(IDFOR-2))
               ELSE
                  ICOLOR = IWHITE
                  IDFORT = 2
               ENDIF
               COLFOR = REAL (ICOLOR)
            ENDIF
         ENDIF

      ELSE
         IF ((INDX .GT. 0) .AND. (NCOL .GT. 0)) THEN
C ... This is a temporary kludge by GDS for Ultrix -
            ICOLOR = MOD (INDX-1, NCOL) +1
C            ICOLOR = MOD (INDX-1, NCOL)
            IF (MAPUSE(ICURDV) .EQ. 0) THEN
               IF (ICOLOR .GE. IBLACK) ICOLOR = ICOLOR + 1
               IF (ICOLOR .GE. IWHITE) ICOLOR = ICOLOR + 1
               IF (ICOLOR .GT. IBLACK .AND.
     &             ICOLOR .LT. IWHITE) ICOLOR = INT(COLMAP(ICOLOR))
            ELSE
               ICOLOR = ICOLOR + (6 + 2)
            END IF
         ELSE IF (INDX .EQ. -1) THEN
            ICOLOR = IBLACK
         ELSE
            IDFORT = IDFOR
            IF (IDFOR .EQ. 1) THEN
               ICOLOR = IBLACK
            ELSE IF (IDFOR .EQ. 2) THEN
               ICOLOR = IWHITE
            ELSE IF (COLMAP(IDFOR-2) .GE. 0.0) THEN
               ICOLOR = INT(COLMAP(IDFOR-2))
            ELSE
               ICOLOR = IWHITE
               IDFORT = 2
            ENDIF
            COLFOR = REAL (ICOLOR)
         END IF
      ENDIF

      IF (ICOLOR .NE. LSTCOL) THEN
         LDUM = PLTSTD (KCOLOR, REAL (ICOLOR))
         LDUM = PLTSTG (KCOLIN, REAL (ICOLOR))
         LDUM = PLTSTG (KCOSYM, REAL (ICOLOR))
         LSTCOL = ICOLOR
      END IF
      RETURN
      END

