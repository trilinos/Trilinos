C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DTLAB (A, NEWSET, NVWTIM, IPTIMS, TIMES,
     &   NAMECO, NAMES, IELBST, NMIN, NMAX, FMIN, FMAX,
     &   ISSNPS, IDNPS, ISSESS, IDESS, ISCNPS, ISCESS,
     &   LIDSP, BLKCOL, *)
C=======================================================================

C   --*** DTLAB *** (DETOUR) Label plot
C   --   Modified by John Glick - 11/9/88
C   --   Written by Amy Gilkey - revised 04/14/88
C   --   D. P. Flanagan, 07/27/82
C   --
C   --DTLAB calls MSLAB to draw the standard mesh plot label, then adds
C   --DETOUR-specific labeling.  MSLAB also draws the axes.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   NEWSET - IN - true iff starting a new plot set
C   --   NVWTIM - IN - the number of times on this plot
C   --   IPTIMS - IN - the plot time indices (starting with current plot)
C   --   TIMES - IN - the database times
C   --   NAMECO - IN - the coordinate names
C   --   NAMES - IN - the variable names
C   --   IELBST - IN - the element block status:
C   --      -1 = OFF, 0 = ON, but not selected, 1 = selected
C   --   NMIN, NMAX - IN - the number of variables values matching the
C   --      minimum and the maximum (for contour plots only)
C   --   FMIN, FMAX - IN - the minimum and maximum nodal value (contour only)
C   --   ISSNPS - IN - the indices of the selected node sets
C   --   IDNPS - IN - the node set IDs
C   --   ISSESS - IN - the indices of the selected side sets
C   --   IDESS - IN - the side set IDs
C   --   ISCNPS - IN/OUT - size = NUMNPS, set iff NEWSET
C   --   ISCESS - IN/OUT - size = NUMESS, set iff NEWSET
C   --   LIDSP(0:*)  - IN - the indices of the selected variables
C   --          whose values will be displayed on the plot legend.
C   --          LIDSP(0) = the number of variables in the list.
C   --          LIDSP(i) identifies the ith variable in the list.
C   --          If LIDSP(i) > 0, LIDSP(i) is the id of a history variable.
C   --          If LIDSP(i) < 0, -LIDSP(i) is the id of a global variable.
C   --          If LIDSP(i) = 0, TIME is to be displayed on the plot legend.
C   --   BLKCOL - IN/OUT - the user selected colors of the element blocks.
C   --                    BLKCOL(0) = 1 if the user defined material
C   --                                colors should be used in mesh plots.
C   --                              = -1 if program selected colors should
C   --                                be used.
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
C   --   * - return statement if the cancel function is active
C   --
C   --Common Variables:
C   --   Uses NDIM of /DBNUMS/
C   --   Uses IS3DIM of /D3NUMS/
C   --   Uses DOLEG of /LEGOPT/
C   --   Uses MODDET, MODTYP, IDTVAR, NNDVAR, NEDVAR of /DETOPT/
C   --   Uses NCNTR, LABINC, MAXMIN, MAXMAX of /CNTR/
C   --   Uses CHLSIZ of /LAYOUT/

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)

      include 'params.blk'
      include 'dbnums.blk'
      include 'dbnumgq.blk'
      include 'legopt.blk'
      include 'd3nums.blk'
      include 'mshopt.blk'
      include 'detopt.blk'
      include 'cntr.blk'
      include 'layout.blk'

      DIMENSION A(*)
      LOGICAL NEWSET
      INTEGER IPTIMS(*)
      REAL TIMES(*)
      CHARACTER*(MXSTLN) NAMECO(*)
      CHARACTER*(*) NAMES(*)
      INTEGER IELBST(NELBLK)
      INTEGER ISSNPS(NUMNPS,4)
      INTEGER ISSESS(NUMESS,4)
      INTEGER IDNPS(*)
      INTEGER IDESS(*)
      INTEGER ISCNPS(*)
      INTEGER ISCESS(*)
      INTEGER LIDSP(0:*)
      INTEGER BLKCOL(0:NELBLK)

      LOGICAL GRABRT
      INTEGER NUMMOD

      LOGICAL SOFTCH
      LOGICAL LABCON
      CHARACTER*164 STRING
      CHARACTER*20 RSTR(3)
      INTEGER INUM(3)
      REAL RNUM(3)

      REAL DLEGND(KTOP)
      SAVE DLEGND

      CHARACTER*(MXSTLN) VNAM(4)
      LOGICAL CONMOD, ISLINE, SPECTR, TWOCON
      SAVE DYVNAM, DYCNTR, DYMNMX,
     &   NVNAM, NLMAX, VNAM,
     &   CONMOD, ISLINE, SPECTR, NOVER, TWOCON, lcntr

C *** Standard mesh plot labeling ***

      CALL MSLAB (A, NEWSET, NVWTIM, IPTIMS, TIMES,
     &   NAMECO, NAMES, IELBST,
     &   ISSNPS, IDNPS, ISSESS, IDESS, ISCNPS, ISCESS,
     &   DLEGND, LIDSP, BLKCOL, *180)

C *** Plot legend ***

      IF (.NOT. DOLEG(1)) GOTO 170

C   --Set up labeling

      IF (NEWSET) THEN

C      --Get the number of variables and the number of lines needed
         NVNAM = MAX (NNDVAR, NEDVAR)
         IF (NVNAM .GT. 0) THEN
            IF (MOD (NVNAM, 3) .EQ. 0) THEN
               NLMAX = 3
            ELSE
               NLMAX = 2
            END IF
            LVNAM = 1 + (NVNAM-1) / NLMAX
         ELSE
            LVNAM = 0
         END IF

         CONMOD = (NUMMOD (MODDET, ' ', 'CONTOUR', ' ') .GE. 1)
     &      .OR. (NUMMOD (MODDET, ' ', 'ELEMCONT', ' ') .GE. 1)

         IF (CONMOD) THEN
C         --Get the number of lines for contours
            ISLINE = (NUMMOD (MODDET, MODTYP, 'CONTOUR', 'LINE') .GE. 1)
     &         .OR. (NUMMOD (MODDET, MODTYP, 'ELEMCONT', 'LINE') .GE. 1)
            IF (ISLINE) THEN
               SPECTR = .FALSE.
               LCNTR = NCNTR
            ELSE
               CALL GRGPAR ('SPECTRUM', 0, INUM, STRING)
               SPECTR = ((INUM(3) .GT. 0) .AND. (INUM(1) .GT. 10))
               IF (.NOT. SPECTR) THEN
                  LCNTR = NCNTR+1
               ELSE
                  LCNTR = 10+1
               END IF
            END IF

C         --Get the number of lines for min/max symbols
            IF ((MAXMIN .GT. 0) .OR. (MAXMAX .GT. 0)) THEN
               LNMNMX = 2
            ELSE
               LNMNMX = 0
            END IF
         ELSE
            LCNTR = 0
            LNMNMX = 0
         END IF

         NCENLN = LVNAM + MIN (LVNAM, 1) + LCNTR + MIN (LCNTR, 1)
     &      + LNMNMX + MIN (LNMNMX, 1) - 1

         IF (NCENLN .GT. 0) THEN
            DYLTOP = DLEGND(KTOP) - 1.5*CHLSIZ
            DYLBOT = DLEGND(KBOT) + 1.5*CHLSIZ
            CALL GRYCEN (CHLSIZ, DYLTOP, DYLBOT, NCENLN, NOVER)
            IF (NVNAM .GT. 0) THEN
               DYVNAM = DYLTOP
               DYLTOP = DYVNAM - (LVNAM+1) * CHLSIZ
            END IF
            IF (CONMOD) THEN
               DYCNTR = DYLTOP
               IF (NOVER .GT. 0) THEN
                  IF (.NOT. SPECTR) THEN
                     LCNTR = LCNTR - NOVER - 2
                  ELSE
                     LCNTR = LCNTR - NOVER
                     NOVER = 0
                  END IF
               END IF
               IF (NOVER .GT. 0) THEN
                  DYLTOP = DYLTOP - (LCNTR+2+1) * CHLSIZ
               ELSE
                  DYLTOP = DYLTOP - (LCNTR+1) * CHLSIZ
               END IF
               IF ((MAXMIN .GT. 0) .OR. (MAXMAX .GT. 0)) THEN
                  DYMNMX = DYLTOP
                  DYLTOP = DYLTOP -
     &               (MIN (MAXMIN, 1) + MIN (MAXMAX, 1) + 1) * CHLSIZ
               END IF
            END IF
         END IF
      END IF

C   --Get software character flag for current device
      CALL GRGPARD ('SOFTCHAR', 0, SOFTCH, STRING)

C   --Variable name(s)

      IF (NVNAM .GT. 0) THEN
         IF (NEWSET) THEN
            DO 100 IV = 1, NVNAM
               IF (IDTVAR(IV) .GT. 0) THEN
                  VNAM(IV) = NAMES(IDTVAR(IV))
               ELSE
                  VNAM(IV) = '0'
               END IF
  100       CONTINUE
         END IF

         DY = DYVNAM
         DO 110 ISTART = 1, NVNAM, NLMAX
            IEND = MIN (ISTART+NLMAX-1, NVNAM)
            WRITE (STRING, '(4 (A, 1X))') (VNAM(IV), IV=ISTART,IEND)
            CALL SQZSTR (STRING, LSTR)
            IF (ISTART .EQ. 1) THEN
               CALL GRTEXT (DLEGND(KLFT), DY, STRING(:LSTR))
            ELSE
               CALL GRTEXT (DLEGND(KLFT), DY, '  ' // STRING(:LSTR))
            END IF
            DY = DY - CHLSIZ
  110    CONTINUE
      END IF

C   --Contour labels

      IF (CONMOD) THEN
C      --Use the selected color table
         CALL GRCOLU ('ALTERNATE')

         IF (ISLINE) THEN
            LABCON = (LABINC .GE. 0)
            IF (LABCON) THEN
               DY = DYCNTR
               DO 120 I = 1, LCNTR
                  IF (GRABRT ()) RETURN 1
                  CALL GRCOLR (I)
                  CALL GRTEXT (DLEGND(KLFT), DY, CHAR(64+I) // ' =')
                  DY = DY - CHLSIZ
  120          CONTINUE
               IF (NOVER .GT. 0) THEN
                  IF (GRABRT ()) RETURN 1
                  CALL UGRCOL (0, BLKCOL)
                  CALL GRTEXT (DLEGND(KLFT), DY, ':')
                  DY = DY - CHLSIZ
                  IF (GRABRT ()) RETURN 1
                  I = NCNTR
                  CALL GRCOLR (I)
                  CALL GRTEXT (DLEGND(KLFT), DY, CHAR(64+I) // ' =')
               END IF
               IF (SOFTCH) THEN
                  CALL PLTXSE (DX, RDUM)
               ELSE
                  CALL PLTXHE (DX, RDUM)
               END IF
            END IF

         ELSE
            LABCON = .TRUE.
            DY = DYCNTR + 0.5*CHLSIZ
            DX = DLEGND(KLFT) + 1.5*CHLSIZ
            IF (.NOT. SPECTR) THEN
               BOXSIZ = CHLSIZ
               NC = LCNTR-1
            ELSE
               BOXSIZ = CHLSIZ * (LCNTR-1) / NCNTR
               NC = NCNTR
            END IF
            DO 130 I = 1, NC
               IF (GRABRT ()) RETURN 1
               CALL GRCOLR (I)
               CALL GRBOX ('S', DLEGND(KLFT), DX, DY-BOXSIZ, DY)
               DY = DY - BOXSIZ
  130       CONTINUE
            IF (NOVER .GT. 0) THEN
               DY = DY - CHLSIZ
               IF (GRABRT ()) RETURN 1
               CALL UGRCOL (0, BLKCOL)
               CALL GRTEXT (DLEGND(KLFT), DY, ':')
               IF (GRABRT ()) RETURN 1
               CALL GRCOLR (NCNTR)
               CALL GRBOX ('S', DLEGND(KLFT), DX, DY-CHLSIZ, DY)
            END IF
         END IF

         IF (LABCON) THEN
            IF (SOFTCH) THEN
               CALL PLTXSL (' 0', RDUM)
               CALL PLTXSL ('0', R)
            ELSE
               CALL PLTXHL (' 0', RDUM)
               CALL PLTXHL ('0', R)
            END IF
            DX = DX + (RDUM - R)
         ELSE
            DX = DLEGND(KLFT)
         END IF

         IF (LABCON) THEN
            CALL UGRCOL (0, BLKCOL)
         ENDIF

         IF (.NOT. SPECTR) THEN
            RNUM(2) = CNTRI (1)
            IF ((RNUM(2) .EQ. 0.0) .AND. (LCNTR .GT. 1))
     &         RNUM(2) = CNTRI (MIN (2, LCNTR))
            RNUM(3) = CNTRI (LCNTR)

            DY = DYCNTR
            DO 140 I = 1, LCNTR
               IF (GRABRT ()) RETURN 1
               RNUM(1) = CNTRI (I)
               IF (ISLINE .OR. .NOT. (
     &            (NOCMIN .AND. (I .EQ. 1)) .OR.
     &            (NOCMAX .AND. (I .GT. NCNTR)))) THEN
                  CALL NUMSTR (3, 4, RNUM, RSTR, LSTR)
                  CALL GRCALN (RSTR(1), LSTR)
                  IF (.NOT. LABCON) CALL GRCOLR (I)
                  CALL GRTEXT (DX, DY, RSTR(1)(:LSTR))
               END IF
               DY = DY - CHLSIZ
  140       CONTINUE
            IF (NOVER .GT. 0) THEN
               I = NCNTR
               IF (.NOT. ISLINE) THEN
                  IF (GRABRT ()) RETURN 1
                  RNUM(1) = CNTRI (I)
                  CALL NUMSTR (3, 4, RNUM, RSTR, LSTR)
                  CALL GRCALN (RSTR(1), LSTR)
                  IF (.NOT. LABCON) CALL GRCOLR (I)
                  CALL GRTEXT (DX, DY, RSTR(1)(:LSTR))
                  I = I + 1
               END IF
               DY = DY - CHLSIZ
               IF (GRABRT ()) RETURN 1
               RNUM(1) = CNTRI (I)
               CALL NUMSTR (3, 4, RNUM, RSTR, LSTR)
               CALL GRCALN (RSTR(1), LSTR)
               IF (.NOT. LABCON) CALL GRCOLR (I)
               CALL GRTEXT (DX, DY, RSTR(1)(:LSTR))
               DY = DY - CHLSIZ
            END IF

         ELSE
            RNUM(1) = CNTRI (1)
            RNUM(2) = (CNTRI (1) + CNTRI (NCNTR+1)) / 2
            RNUM(3) = CNTRI (NCNTR+1)
            CALL NUMSTR (3, 4, RNUM, RSTR, LSTR)
            DY = DYCNTR
            DO 150 I = 1, 3
               IF (GRABRT ()) RETURN 1
               CALL GRCALN (RSTR(I), LSTR)
               CALL GRTEXT (DX, DY, RSTR(I)(:LSTR))
               DY = DY - (LCNTR-1)*CHLSIZ/(3-1)
  150       CONTINUE
         END IF

         IF (.NOT. LABCON) THEN
            CALL UGRCOL (0, BLKCOL)
         ENDIF
C      --Use the standard color table
         CALL GRCOLU ('STANDARD')
      END IF

C   --Nodal variable min/max

      IF (CONMOD .AND. ((MAXMIN .GT. 0) .OR. (MAXMAX .GT. 0))) THEN

C      --Count the number of contour modes if more than one time

         IF (NEWSET) THEN
            IF (NVWTIM .GT. 1) THEN
               TWOCON = ((NUMMOD (MODDET, ' ', 'CONTOUR', ' ')
     &            + NUMMOD (MODDET, ' ', 'ELEMCONT', ' ')) .GT. 1)
            ELSE
               TWOCON = .FALSE.
            END IF
         END IF

         IF ((MAXMIN .GE. NMIN) .OR. (MAXMAX .GE. NMAX)) THEN
            IF (TWOCON) THEN
               RSTR(1) = 'MINIMUM'
               RSTR(2) = 'MAXIMUM'
            ELSE
               RNUM(1) = FMIN
               RNUM(2) = FMAX
               RSTR(1) = ' '
               RSTR(2) = ' '
               IF ((MAXMIN .GE. NMIN) .AND. (MAXMAX .GE. NMAX)) THEN
                  CALL NUMSTR (2, 4, RNUM, RSTR, LSTR)
               ELSE IF (MAXMIN .GE. NMIN) THEN
                  CALL NUMSTR1(4, RNUM(1), RSTR(1), LSTR)
               ELSE IF (MAXMAX .GE. NMAX) THEN
                  CALL NUMSTR1(4, RNUM(2), RSTR(2), LSTR)
               END IF
            END IF
            DY = DYMNMX
            DO 160 I = 1, 2
               IF (RSTR(I) .NE. ' ') THEN
                  IF (GRABRT ()) RETURN 1
                  CALL GRCALN (RSTR(I), LSTR)
#if NeedsDoubleEscape
                  IF (I .EQ. 1) THEN
                     CALL PLTXTS (DLEGND(KLFT), DY, '\\CI')
                  ELSE
                     CALL PLTXTS (DLEGND(KLFT), DY, '\\X')
                  END IF
                  CALL PLTXTS (DLEGND(KLFT), DY, '\\CS')
#else
                  IF (I .EQ. 1) THEN
                     CALL PLTXTS (DLEGND(KLFT), DY, '\CI')
                  ELSE
                     CALL PLTXTS (DLEGND(KLFT), DY, '\X')
                  END IF
                  CALL PLTXTS (DLEGND(KLFT), DY, '\CS')
#endif
                  CALL PLTXSE (DX, RDUM)
                  CALL GRTEXT (DX, DY, ' = ' // RSTR(I)(:LSTR))
                  DY = DY - CHLSIZ
               END IF
  160       CONTINUE
            DY = DY - CHLSIZ
         END IF
      END IF

  170 CONTINUE

C   --Flush buffer, so label is complete at this point
      CALL PLTFLU

      RETURN

  180 CONTINUE
      RETURN 1
      END
