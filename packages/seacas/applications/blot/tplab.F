C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE TPLAB (IPVAR, NCRV, NUMCRV, TIMLIM, NAMES,
     &   TXLAB, TYLAB, BLKCOL, MAPEL, MAPND, *)
C=======================================================================

C   --*** TPLAB *** (TPLOT) Label plot
C   --   Written by Amy Gilkey - revised 05/02/88
C   --
C   --TPLAB erases the display surface then labels the plot legend
C   --with the following:
C   --
C   --   BANNER - the title, creator, modifier, and drawn-by information.
C   --      Constant for entire run.  If QA only.
C   --   CAPTION - the plot caption.
C   --   LABEL - the plot variable names.
C   --
C   --It is assumed that time curves and X-Y curves are not overlaid.
C   --
C   --The vertical-to-horizontal aspect ratio is assumed to be 0.75.
C   --Label locations are set for a horizontal layout.
C   --
C   --Parameters:
C   --   IPVAR - IN - the /TPVARS/ index of the starting plot variable
C   --   NCRV - IN - the number of curves on the plot
C   --   NUMCRV - IN - true iff curves are to be numbered
C   --   TIMLIM - IN - the starting and ending time of a
C   --      variable-versus-variable plot
C   --   NAMES - IN - the variable names
C   --   TXLAB, TYLAB - OUT - labels for the X and Y axis
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
C   --   * - the return statement if the cancel function is active
C   --
C   --Common Variables:
C   --   Uses TITLE, CREATE, MODIFY, DRAW of /DBTITL/
C   --   Uses TIMPLT, ITVID, ITVNE of /TPVARS/
C   --   Uses DOQA, DOLEG, CAPTN of /LEGOPT/
C   --   Uses LINTYP, ISYTYP of /XYOPT/
C   --   Uses XLAB, YLAB of /XYLAB/
C   --   Uses CHLSIZ, DBORD0, DVIEW0 of /LAYOUT/

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4)

      PARAMETER (NUMSYM = 6, NUMLIN = 6)

      include 'params.blk'
      include 'dbnums.blk'
      include 'dbtitl.blk'
      include 'tpvars.blk'
      include 'legopt.blk'
      include 'xyopt.blk'
      include 'xylab.blk'
      include 'layout.blk'
      include 'debug.blk'

      LOGICAL NUMCRV
      REAL TIMLIM(2)
      CHARACTER*(*) NAMES(*)
      CHARACTER*(*) TXLAB, TYLAB
      INTEGER BLKCOL(0:NELBLK)
      INTEGER MAPEL(*), MAPND(*)

      LOGICAL GRABRT
      CHARACTER*80 LABSTR
      REAL DLEGND(KTOP)
      CHARACTER*20 RSTR(2)
      CHARACTER*10 FMT

      CHARACTER*3 SYMBOL(NUMSYM)
      CHARACTER*7 LTYPES(NUMLIN)
      SAVE SYMBOL, LTYPES

#if NeedsDoubleEscape
      DATA SYMBOL / '\\SQ', '\\DI', '\\CS', '\\X ', '\\TR', '\\CI' /
      DATA LTYPES / '\\SLINE ', '\\DLINE ', '\\DDLINE',
     &   '\\SDLINE', '\\LDLINE', '\\MDLINE' /
#else
      DATA SYMBOL / '\SQ', '\DI', '\CS', '\X ', '\TR', '\CI' /
      DATA LTYPES / '\SLINE ', '\DLINE ', '\DDLINE',
     &   '\SDLINE', '\LDLINE', '\MDLINE' /
#endif

C   --Set up axes labeling (needed before standard labeling)

      IF (XLAB .NE. ' ') THEN
         TXLAB = XLAB
      ELSE IF (TIMPLT) THEN
         TXLAB = 'TIME'
      ELSE
         N = IPVAR
         ID = ITVID(N)
         DO 100 NP = 2, NCRV
            N = N + 2
            IF (ITVID(N) .NE. ID) ID = 0
  100    CONTINUE
         IF (ID .NE. 0) THEN
            TXLAB = NAMES(ID)
         ELSE
            TXLAB = ' '
         END IF
      END IF

      IF (YLAB .NE. ' ') THEN
         TYLAB = YLAB
      ELSE
         N = IPVAR
         IF (.NOT. TIMPLT) N = N + 1
         ID = ITVID(N)
         DO 110 NP = 2, NCRV
            N = N + 1
            IF (.NOT. TIMPLT) N = N + 1
            IF (ITVID(N) .NE. ID) ID = 0
  110    CONTINUE
         IF (ID .NE. 0) THEN
            TYLAB = NAMES(ID)
         ELSE
            TYLAB = ' '
         END IF
      END IF

C *** Standard plot labeling ***

      CALL QALAB (DBORD0, DVIEW0, CHLSIZ,
     &   DOQA(2), DOAXIS(2), (TXLAB .NE. ' '),
     &   CAPTN(1,2), TITLE, CREATE, MODIFY, DRAW, DLEGND,
     &   BLKCOL,  *140)

C *** Plot legend ***

      IF (.NOT. DOLEG(2)) GOTO 130

C   --Find middle of the remaining legend area and set limits

      IF (TIMPLT) THEN
         NLINE = NCRV
      ELSE
         NLINE = 2*NCRV + 3
      END IF
      DTOP = DLEGND(KTOP) - CHLSIZ
      DBOT = DLEGND(KBOT)
      IF (DOQA(2)) THEN
         DTOP = DTOP - 2*CHLSIZ
         DBOT = DBOT + 2*CHLSIZ
      END IF
      CALL GRYCEN (CHLSIZ, DTOP, DBOT, NLINE, NOVER)
      IF (DOQA(2) .AND. (NOVER .GT. 0)) NOVER = NOVER - 2
C      --Allow 2 extra lines if QA (DBOT is invalid)

C   --Display plot item(s) (variable name and number) for each curve

      NC = NCRV
      IF (NOVER .GT. 0) THEN
         IF (TIMPLT) THEN
            NC = NCRV - NOVER - 1
         ELSE
            NC = NCRV - (NOVER+1)/2 - 1
         END IF
      END IF

      N = IPVAR
      DO 120 NP = 1, NC
         IF (GRABRT ()) RETURN 1
         CALL GRCOLR (NP)

         IF (NUMCRV) THEN
            IF (NP .LE. 1) THEN
               IF (NC .LT. 10) THEN
                  FMT = '(I1)'
               ELSE
                  FMT = '(I2)'
               END IF
            END IF
            WRITE (LABSTR, FMT, IOSTAT=IDUM) NP
            CALL GRCALN (LABSTR, L)
            CALL TPLABV (0, ITVID(N), NAMES(ITVID(N)), ITVNE(N),
     &         LABSTR(L+2:), MAPEL, MAPND)
         ELSE
            CALL TPLABV (0, ITVID(N), NAMES(ITVID(N)), ITVNE(N), LABSTR,
     *       MAPEL, MAPND)
         END IF

         IF (GRABRT ()) RETURN 1
         IF ((NCRV .LE. 1)
     &      .OR. ((ISYTYP .GE. 0) .AND. (LINTYP .GE. 0))) THEN
            CALL GRTEXT (DLEGND(KLFT), DTOP, LABSTR)
         ELSE
            CALL PLTXSE (DX, RDUM)
            IF (LINTYP .LT. 0) THEN
               ITYP = MOD (NP-1, NUMLIN) + 1
               CALL PLTXTS (DLEGND(KLFT), DTOP,
     &            LTYPES(ITYP)(:LENSTR(LTYPES(ITYP))))
               CALL PLTXSE (DX, RDUM)
            END IF
            IF (ISYTYP .LT. 0) THEN
               ISYM = MOD (NP-1, NUMSYM) + 1
               IF (LINTYP .GE. 0) THEN
                  CALL PLTXTS (DLEGND(KLFT), DTOP,
     &               SYMBOL(ISYM)(:LENSTR(SYMBOL(ISYM))))
                  CALL PLTXSE (DX, RDUM)
               ELSE
                  DXX = 0.5 * (DLEGND(KLFT) + DX)
     &               - 0.75 * 5.0/7.0*CHLSIZ
                  CALL PLTXTS (DXX, DTOP,
     &               SYMBOL(ISYM)(:LENSTR(SYMBOL(ISYM))))
               END IF
            END IF
            CALL GRTEXT (DX, DTOP, ' ' // LABSTR)
         END IF
         DTOP = DTOP - CHLSIZ
         N = N + 1

         IF (.NOT. TIMPLT) THEN
            IF (GRABRT ()) RETURN 1
            CALL TPLABV (0, ITVID(N), NAMES(ITVID(N)), ITVNE(N), LABSTR,
     *        MAPEL, MAPND)
            CALL GRTEXT (DLEGND(KLFT), DTOP, '  ' // LABSTR)
            DTOP = DTOP - CHLSIZ
            N = N + 1
         END IF
  120 CONTINUE

      CALL UGRCOL (0, BLKCOL)

      IF (NC .LT. NCRV) THEN
         IF (GRABRT ()) RETURN 1
         CALL GRTEXT (DLEGND(KLFT), DTOP, '...')
         DTOP = DTOP - CHLSIZ
      END IF

C   --Display times for variable-versus-variable plot

      DTOP = DTOP - CHLSIZ

      IF (.NOT. TIMPLT) THEN
         CALL NUMSTR (2, 4, TIMLIM, RSTR, LSTR)
         CALL PCKSTR (2, RSTR)
         CALL GRTEXT (DLEGND(KLFT), DTOP, 'TIMES ' // RSTR(1))
         DTOP = DTOP - CHLSIZ
         CALL GRTEXT (DLEGND(KLFT), DTOP, '  TO ' // RSTR(2))
         DTOP = DTOP - CHLSIZ
      END IF

  130 CONTINUE

C   --Flush buffer, so label is complete at this point
      CALL PLTFLU

      RETURN

  140 CONTINUE
      RETURN 1
      END
