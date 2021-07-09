C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SPLAB (A, NPTIMS, IPTIMS, TIMES, NENUM,
     &   IPVAR, NCURVE, NUMCRV, NAMES, TXLAB, TYLAB, LIDSP,
     &   BLKCOL,  MAPEL, MAPND, *)
C=======================================================================

C   --*** SPLAB *** (SPLOT) Label plot
C   --   Modified by John Glick - 11/9/88
C   --   Written by Amy Gilkey - revised 05/02/88
C   --
C   --SPLAB erases the display surface then labels the plot legend
C   --with the following:
C   --
C   --   BANNER - the title, creator, modifier, and drawn-by information.
C   --      Constant for entire run.  If QA only.
C   --   CAPTION - the plot caption.
C   --   LABEL - the node/element range, the plot variable names, the plot time.
C   --
C   --The vertical-to-horizontal aspect ratio is assumed to be 0.75.
C   --Label locations are set for a horizontal layout.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   NPTIMS - IN - the number of selected times
C   --   IPTIMS - IN - the selected time steps
C   --   TIMES - IN - the database times
C   --   NENUM - IN - the selected node/element numbers
C   --   IPVAR - IN - the /SPVARS/ index of the starting plot variable
C   --   NCURVE - IN - the number of curves on the plot
C   --   NUMCRV - IN - true iff curves are to be numbered
C   --   NAMES - IN - the variable names
C   --   TXLAB, TYLAB - OUT - labels for the X and Y axis
C   --   LIDSP(0:*)  - IN/OUT - the indices of the selected variables
C   --          whose values will be displayed on the plot legend.
C   --          ABS(LIDSP(0)) = the number of variables in the list.
C   --          SIGN(LIDSP(0)) specifies whether the variables in the
C   --                   list should have their values displayed on
C   --                   the plot legend.  If >0, they should;
C   --                   If <=0, they should not.
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
C   --   * - the return statement if the cancel function is active
C   --
C   --Common Variables:
C   --   Uses TITLE, CREATE, MODIFY, DRAW of /DBTITL/
C   --   Uses NODVAR, NNENUM of /SELNE/
C   --   Uses ISVID of /SPVARS/
C   --   Uses DOQA, DOLEG, CAPTN of /LEGOPT/
C   --   Uses LINTYP, ISYTYP, OVERLY, OVERTM of /XYOPT/
C   --   Uses XLAB, YLAB of /XYLAB/
C   --   Uses CHLSIZ, DBORD0, DVIEW0 of /LAYOUT/

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4)

      PARAMETER (NUMSYM = 6, NUMLIN = 6)

      include 'params.blk'
      include 'dbtitl.blk'
      include 'selne.blk'
      include 'spvars.blk'
      include 'legopt.blk'
      include 'xyopt.blk'
      include 'xylim.blk'
      include 'xylab.blk'
      include 'layout.blk'
      include 'dbnums.blk'

      common /debugc/ cdebug
      common /debugn/ idebug
      character*8 cdebug

      DIMENSION A(*)
      INTEGER IPTIMS(NPTIMS)
      REAL TIMES(*)
      INTEGER NENUM(NNENUM)
      LOGICAL NUMCRV
      CHARACTER*(*) NAMES(*)
      CHARACTER*(*) TXLAB, TYLAB
      INTEGER MAPEL(*), MAPND(*)

      LOGICAL GRABRT
      CHARACTER*64 LABSTR
      CHARACTER*20 RSTR(3)
      CHARACTER*20 FMT
      REAL RNUM(3)
      REAL DLEGND(KTOP)
      INTEGER LIDSP(0:*)
      INTEGER BLKCOL(0:NELBLK)

      LOGICAL FIRSTH, FIRSTG
      CHARACTER*(MXSTLN+1) DISPVL

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
      ELSE
         TXLAB = 'DISTANCE'
      END IF
      IF (YLAB .NE. ' ') THEN
         TYLAB = YLAB
      ELSE IF (.NOT. OVERLY) THEN
         TYLAB = NAMES(ISVID(IPVAR))
      ELSE
         TYLAB = ' '
      END IF

C *** Standard plot labeling ***

      CALL QALAB (DBORD0, DVIEW0, CHLSIZ,
     &   DOQA(2), DOAXIS(2), (TXLAB .NE. ' '),
     &   CAPTN(1,2), TITLE, CREATE, MODIFY, DRAW, DLEGND,
     &   BLKCOL,  *140)

C *** Plot legend ***

      IF (.NOT. DOLEG(2)) GOTO 130

C   --Find middle of the remaining legend area and set limits

      NUMDSP = LIDSP(0)
      IF (NUMDSP .GT. 0) THEN
         NDSPL = NUMDSP * NPTIMS
      ELSE
         NDSPL = 0
      ENDIF
      NLINE = 1 + 1 + NCURVE + 1 + NDSPL
      DTOP = DLEGND(KTOP) - CHLSIZ
      DBOT = DLEGND(KBOT)
      IF (DOQA(2)) THEN
         DTOP = DTOP - 2*CHLSIZ
         DBOT = DBOT + 2*CHLSIZ
      END IF
      CALL GRYCEN (CHLSIZ, DTOP, DBOT, NLINE, NOVER)
      IF (DOQA(2) .AND. (NOVER .GT. 0)) NOVER = NOVER - 2
C      --Allow 2 extra lines if QA (DBOT is invalid)

C   --Display node/element numbers

      IF (GRABRT ()) RETURN 1
      IF (NODVAR) THEN
        WRITE (RSTR(1), 10000, IOSTAT=IDUM) MAPND(NENUM(1)),
     *    MAPND(NENUM(NNENUM))
        CALL PCKSTR (1, RSTR(1))
        CALL GRTEXT (DLEGND(KLFT), DTOP, 'NODES ' // RSTR(1))
      ELSE
        WRITE (RSTR(1), 10000, IOSTAT=IDUM) MAPEL(NENUM(1)),
     *    MAPEL(NENUM(NNENUM))
        CALL PCKSTR (1, RSTR(1))
         CALL GRTEXT (DLEGND(KLFT), DTOP, 'ELEMENTS ' // RSTR(1))
      END IF
      DTOP = DTOP - CHLSIZ

C   --Display plot item(s) (variable name and number) for each curve

      DTOP = DTOP - CHLSIZ

      NCRV = NCURVE
      IF (NOVER .GT. 0) THEN
         IF (NCURVE .GT. 1) THEN
            NCRV = MIN (NCURVE, NLINE-4)
            IF (NCRV .LT. NCURVE) NCRV = NCRV - 1
         ELSE
            NCRV = 1
         ENDIF
      ENDIF
      NUSED = NCRV + 3
      IF (NCRV .LT. NCURVE) NUSED = NUSED + 1

      N = IPVAR
      DO 100 NP = 1, NCRV
         IF (GRABRT ()) RETURN 1
         IF (.NOT. OVERTM) CALL GRCOLR (NP)

         IF ((.NOT. NUMCRV) .OR. (.NOT. OVERLY)) THEN
            LABSTR = NAMES(ISVID(N))
         ELSE
            IF (NP .LE. 1) THEN
               IF (NCRV .LT. 10) THEN
                  FMT = '(I1, 1X, A)'
               ELSE
                  FMT = '(I2, 1X, A)'
               END IF
            END IF
            WRITE (LABSTR, FMT, IOSTAT=IDUM) NP, NAMES(ISVID(N))
            CALL GRCALN (LABSTR, L)
         END IF

         IF ((.NOT. OVERLY)
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
  100 CONTINUE

      IF (.NOT. OVERTM) CALL UGRCOL (0, BLKCOL)

      IF (NCRV .LT. NCURVE) THEN
         IF (GRABRT ()) RETURN 1
         CALL GRTEXT (DLEGND(KLFT), DTOP, '...')
         DTOP = DTOP - CHLSIZ
      END IF

C   --Show display variables for each curve

      DTOP = DTOP - CHLSIZ

      IF (NPTIMS .GT. 1) THEN
         RNUM(2) = TIMES(IPTIMS(1))
         IF ((RNUM(2) .EQ. 0.0) .AND. (NPTIMS .GT. 1))
     &      RNUM(2) = TIMES(IPTIMS(2))
         RNUM(3) = TIMES(IPTIMS(NPTIMS))
      END IF

      NLEFT = NLINE - NUSED
      IF (NDSPL .GT. NLEFT) THEN
         NDSP = NLEFT - 1
      ELSE
         NDSP = NDSPL
      ENDIF
      NTIM = NDSP / NUMDSP
      IF (MOD(NDSP, NUMDSP) .GT. 0) NTIM = NTIM + 1

      NP = 1
      NL = 1
  110 CONTINUE
      IF (NL .LE. NDSP) THEN
         IF (GRABRT ()) RETURN 1
         IF (OVERTM) CALL GRCOLR (NP)

         FIRSTH = .TRUE.
         FIRSTG = .TRUE.
         ND = 1
  120    CONTINUE
         IF ((ND .LE. NUMDSP)  .AND.  (NL .LE. NDSP)) THEN
            IF (LIDSP(ND) .EQ. 0) THEN
               DISPVL = 'TIME'
               VARVAL = TIMES(IPTIMS(NP))
            ELSE IF (LIDSP(ND) .GT. 0)THEN
               IF (FIRSTH) THEN
                  CALL DBVIX_BL ('H', 1, IXHV)
                  CALL MDRSRV ('SCRH', KVARH, NVARHI)
                  CALL MDSTAT (NERR, MEM)
                  IF (NERR .GT. 0)GOTO 130
                  CALL GETVAR (A, IXHV, 0, IPTIMS(NP), NVARHI,
     &               A(KVARH))
                  FIRSTH = .FALSE.
               ENDIF
               CALL DBVIX_BL ('H', LIDSP(ND), IDVAR)
               DISPVL = NAMES(IDVAR)
               VARVAL = A(KVARH+LIDSP(ND)-1)
            ELSE
               IF (FIRSTG) THEN
                  CALL DBVIX_BL ('G', 1, IXGV)
                  CALL MDRSRV ('SCRG', KVARG, NVARGL)
                  CALL MDSTAT (NERR, MEM)
                  IF (NERR .GT. 0)GOTO 130
                  CALL GETVAR (A, IXGV, 0, IPTIMS(NP), NVARGL,
     &               A(KVARG))
                  FIRSTG = .FALSE.
               ENDIF
               CALL DBVIX_BL ('G', -LIDSP(ND), IDVAR)
               DISPVL = NAMES(IDVAR)
               VARVAL = A(KVARG-LIDSP(ND)-1)
            ENDIF

            RNUM(1) = VARVAL
            CALL NUMSTR (MIN(NPTIMS,3), 4, RNUM, RSTR, LSTR)
            IF ((.NOT. NUMCRV) .OR. (.NOT. OVERTM)) THEN
               LABSTR = DISPVL(:LENSTR(DISPVL)) // '  ' // RSTR(1)
            ELSE
               IF (NP .LE. 1) THEN
                  IF (NTIM .LT. 10) THEN
                     FMT = '(I1, 1X, A, A)'
                  ELSE
                     FMT = '(I2, 1X, A, A)'
                  END IF
               END IF
               WRITE (LABSTR, FMT, IOSTAT=IDUM) NP,
     *           DISPVL(:LENSTR(DISPVL)) // '  ' // RSTR(1)
            END IF
            CALL GRCALN (LABSTR, LSTR)

            IF ((.NOT. OVERTM)
     &         .OR. ((ISYTYP .GE. 0) .AND. (LINTYP .GE. 0))) THEN
               CALL GRTEXT (DLEGND(KLFT), DTOP, LABSTR)
            ELSE
               IF (ISYTYP .LT. 0) THEN
                  ISYM = MOD (NP-1, NUMSYM) + 1
                  CALL PLTXTS (DLEGND(KLFT), DTOP,
     &               SYMBOL(ISYM)(:LENSTR(SYMBOL(ISYM))))
               ELSE IF (LINTYP .LT. 0) THEN
                  ITYP = MOD (NP-1, NUMLIN) + 1
                  CALL PLTXTS (DLEGND(KLFT), DTOP,
     &               LTYPES(ITYP)(:LENSTR(LTYPES(ITYP))))
               END IF
               CALL PLTXSE (DX, RDUM)
               CALL GRTEXT (DX, DTOP, ' ' // LABSTR)
            END IF
            DTOP = DTOP - CHLSIZ
            N = N + 1

            NL = NL + 1
            ND = ND + 1
            GO TO 120
         ENDIF

         IF (.NOT. FIRSTH) CALL MDDEL ('SCRH')
         IF (.NOT. FIRSTG) CALL MDDEL ('SCRG')
         NP = NP + 1
         GO TO 110
      ENDIF

      IF (OVERTM) CALL UGRCOL (0, BLKCOL)

      IF (NDSP .LT. NDSPL) THEN
         IF (GRABRT ()) RETURN 1
         CALL GRTEXT (DLEGND(KLFT), DTOP, '...')
         DTOP = DTOP - CHLSIZ
      END IF

  130 CONTINUE

C   --Flush buffer, so label is complete at this point
      CALL PLTFLU

      RETURN

  140 CONTINUE
      RETURN 1

10000  FORMAT (I6, '..', I6)
      END
