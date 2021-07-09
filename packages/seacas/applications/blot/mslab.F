C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MSLAB (A, NEWSET, NVWTIM, IPTIMS, TIMES,
     &   NAMECO, NAMES, IELBST,
     &   ISSNPS, IDNPS, ISSESS, IDESS, ISCNPS, ISCESS,
     &   DLEGND, LIDSP, BLKCOL, *)
C=======================================================================

C   --*** MSLAB *** (MESH) Label plot
C   --   Modified by John Glick - 11/9/88
C   --   Written by Amy Gilkey - revised 05/26/88
C   --   D. P. Flanagan, 07/27/82
C   --
C   --MSLAB erases the display surface then labels the graphic frame
C   --with the following graphics structures:
C   --
C   --   BANNER - the title, creator, modifier, and user information.
C   --      Constant for entire run.  If QA only.
C   --   LABEL - includes the variable name, contour labels, magnification
C   --      factor, plot time, and axis orientation indicator.  Except
C   --      for plot time and the min/max nodal values, this structure
C   --      is constant for a particular plot set.
C   --   CAPTION - the plot caption.
C   --   AXES - tick-marked with numeric labels and text labels.
C   --
C   --All view windows are outlined.
C   --
C   --The vertical-to-horizontal aspect ratio is assumed to be 0.75.
C   --Label locations are set for a horizontal layout.
C   --
C   --The lengths and minimum/maximums of the axes must be set before
C   --this routine is called.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   NEWSET - IN - true iff starting a new plot set
C   --   NVWTIM - IN - the number of times on this plot; 0 for no time
C   --   IPTIMS - IN - the plot time indices (starting with current plot)
C   --   TIMES - IN - the database times
C   --   NAMECO - IN - the coordinate names
C   --   NAMES - IN - the variable names
C   --   IELBST - IN - the element block status:
C   --      -1 = OFF, 0 = ON, but not selected, 1 = selected
C   --   ISSNPS - IN - the indices of the selected node sets
C   --   IDNPS - IN - the node set IDs
C   --   ISSESS - IN - the indices of the selected side sets
C   --   IDESS - IN - the side set IDs
C   --   ISCNPS - IN/OUT - size = NUMNPS, set iff NEWSET
C   --   ISCESS - IN/OUT - size = NUMESS, set iff NEWSET
C   --   DLEGND - OUT - the limits of the special legend area
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
C   --   * - return statement if the cancel function is active
C   --
C   --Common Variables:
C   --   Uses TITLE, CREATE, MODIFY, DRAW of /DBTITL/
C   --   Uses NDIM, NELBLK of /DBNUMS/
C   --   Uses IS3DIM of /D3NUMS/
C   --   Uses DOQA, DOLEG, DOAXIS, CAPTN of /LEGOPT/
C   --   Uses DFAC of /DEFORM/
C   --   Uses MSHDEF, NALVAR, DEADNP of /MSHOPT/
C   --   Uses XISSYM, YISSYM of /VIEWS/
C   --   Uses TICMSH of /MSHLIM/
C   --   Uses ROTMAT of /ROTOPT/
C   --   Uses CHLSIZ, DBORD0, DVIEW0, DVIEW, WVIEW of /LAYOUT/

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)

      PARAMETER (KXORIG=1, KYORIG=2, KXLENG=3, KYLENG=4)
      PARAMETER (KXNUMS=22, KXLABS=23, KYNUMS=47, KYLABS=48)
      PARAMETER (KSCALE=11)
      PARAMETER (KTICSZ=33)

      include 'params.blk'
      include 'dbtitl.blk'
      include 'dbnums.blk'
      include 'dbnumgq.blk'
      include 'd3nums.blk'
      include 'legopt.blk'
      include 'deform.blk'
      include 'mshopt.blk'
      include 'views.blk'
      include 'mshlim.blk'
      include 'rotopt.blk'
      include 'layout.blk'
      include 'layoud.blk'

      common /debugc/ cdebug
      common /debugn/ idebug
      character*8 cdebug

      DIMENSION A(*)
      LOGICAL NEWSET
      INTEGER IPTIMS(*)
      REAL TIMES(*)
      CHARACTER*(*) NAMECO(*)
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
      REAL DLEGND(KTOP)

      LOGICAL PLTGTG, PLTSTG, PLTSTG1, MPVIEW, MPORT2, LDUM
      LOGICAL GRABRT
      INTEGER NUMMOD, NDEFVW, IXVW

      LOGICAL SOFTCH
      LOGICAL XISCOP, YISCOP
      LOGICAL RGTVW, TOPVW
      CHARACTER*80 STRING
      CHARACTER*20 RSTR(3)
      REAL BUF(11)
      REAL X(4)
      REAL XPT(0:3), YPT(0:3), ZPT(0:3)
      REAL WAVIEW(KTOP), DAVIEW(KTOP)

      REAL ZERO3(3)
      SAVE ZERO3

      SAVE NNPS, NESS
      SAVE DTOP, DBOT, DYMAGN, DYELB, DYALIV, DYSETS,
     &   DYTIME, EYTIME, DYXYZ, EYXYZ

      CHARACTER*(MXSTLN+1) DISPVL
      REAL VARVAL
      INTEGER IDVAR
      INTEGER NUMDSP
      LOGICAL FIRSTH, FIRSTG

      DATA ZERO3 / 0.0, 0.0, 0.0 /

      RGTVW(I) = (I .EQ. 2) .OR. (I .EQ. 4)
      TOPVW(I) = (I .EQ. 1) .OR. (I .EQ. 2)

C *** Standard plot labeling ***

      CALL QALAB (DBORD0, DVIEW0, CHLSIZ,
     &   DOQA(1), DOAXIS(1), (.NOT. IS3DIM),
     &   CAPTN(1,1), TITLE, CREATE, MODIFY, DRAW, DLEGND,
     &   BLKCOL,  *220)

C *** Plot legend ***

      IF (.NOT. DOLEG(1)) GOTO 200

C   --Set up labeling

      IF (NEWSET) THEN
         DTOP = DLEGND(KTOP) - CHLSIZ
         IF (DOQA(1)) THEN
            DTOP = DTOP - 2*CHLSIZ
         END IF

         IF (NVWTIM .GT. 0) THEN
            DYMAGN = DTOP
            DTOP = DYMAGN - 2*CHLSIZ
         ELSE
            DYMAGN = 0.0
         END IF

         DYELB = DTOP
         DTOP = DYELB - 3*CHLSIZ

         IF (NALVAR .GT. 0) THEN
            DYALIV = DTOP
            IF (DEADNP) THEN
               DTOP = DYALIV - 3.5*CHLSIZ
            ELSE
               DTOP = DYALIV - 2*CHLSIZ
            END IF
         ELSE
            DYALIV = 0.0
         END IF

C      --Get the number of lines and the identifiers for node sets
C      --and side sets
         NNPS = 0
         NESS = 0
         DO 120 IVW = 1, NDEFVW (.FALSE.)
            IVIEW = IXVW (.FALSE., IVW)
            DO 100 IX = 1, NNPSET(IVIEW)
               I = ISSNPS(IX,IVIEW)
               IF (LOCINT (I, NNPS, ISCNPS) .LE. 0) THEN
                  NNPS = NNPS + 1
                  ISCNPS(NNPS) = I
               END IF
  100       CONTINUE
            DO 110 IX = 1, NESSET(IVIEW)
               I = ISSESS(IX,IVIEW)
               IF (LOCINT (I, NESS, ISCESS) .LE. 0) THEN
                  NESS = NESS + 1
                  ISCESS(NESS) = I
               END IF
  110       CONTINUE
  120    CONTINUE

         IF ((NNPS .GT. 0) .OR. (NESS .GT. 0)) THEN
            DYSETS = DTOP
            NL = MIN (NNPS, 1) + NNPS + MIN (NNPS, 1)
     &         + MIN (NESS, 1) + NESS + MIN (NESS, 1) - 1
            DTOP = DTOP - NL * CHLSIZ
         END IF

         DBOT = DLEGND(KBOT)
         IF (IS3DIM) THEN
            EYXYZ = DBOT + 0.5*CHLSIZ
            DYXYZ = EYXYZ + .1 - 0.5*CHLSIZ
            DBOT = DYXYZ + 2.5*CHLSIZ
         ELSE
            EYXYZ = 0.0
            DYXYZ = 0.0
         END IF
         IF (NVWTIM .GT. 0) THEN
            EYTIME = DBOT
            DYTIME = EYTIME + (NVWTIM-1) * CHLSIZ
            DBOT = DYTIME + CHLSIZ
         ELSE
            EYTIME = 0.0
            DYTIME = 0.0
         END IF
      END IF

C   --Set the top and bottom of the remaining legend area
      DLEGND(KTOP) = DTOP
      DLEGND(KBOT) = DBOT

C   --Get software character flag for current device
      CALL GRGPARD ('SOFTCHAR', 0, SOFTCH, STRING)

C   --Magnification factor

      IF (NVWTIM .GT. 0) THEN
         IF (GRABRT ()) RETURN 1
         IF ((NUMMOD (MSHDEF, ' ', 'DEFORM', ' ') .GE. 1)
     &      .AND. (DFAC .NE. 0.0)) THEN
            CALL NUMSTR (1, 4, DFAC, RSTR(1), LSTR)
            CALL GRTEXT (DLEGND(KLFT), DYMAGN,
     &         'MAGNIFIED BY ' // RSTR(1))
         ELSE
            CALL GRTEXT (DLEGND(KLFT), DYMAGN,
     &         'NO MAGNIFICATION')
         END IF
      END IF

C   --Active element blocks

      IF (GRABRT ()) RETURN 1
      CALL CNTELB (IELBST, NELBLK, NUMON, NUMSEL)
      CALL GRTEXT (DLEGND(KLFT), DYELB, 'ELEMENT BLOCKS ACTIVE:')
      WRITE (STRING, 10000, IOSTAT=IDUM) NUMSEL, NELBLK
10000  FORMAT (I5, ' OF ', I5)
      CALL SQZSTR (STRING, LSTR)
      CALL GRTEXT (DLEGND(KLFT), DYELB-CHLSIZ, '  ' // STRING(:LSTR))

C   --Birth/death variable indicator

      IF (NALVAR .GT. 0) THEN
         IF (GRABRT ()) RETURN 1
         CALL GRTEXT (DLEGND(KLFT), DYALIV,
     &      'BIRTH/DEATH: ' // NAMES(NALVAR))
         IF (DEADNP) THEN
            IF (GRABRT ()) RETURN 1
            DY = DYALIV - 1.5*CHLSIZ
#if NeedsDoubleEscape
            CALL PLTXTS (DLEGND(KLFT), DY+.4*CHLSIZ, '\\CDO')
#else
            CALL PLTXTS (DLEGND(KLFT), DY+.4*CHLSIZ, '\CDO')
#endif
            CALL PLTXSE (DX, RDUM)
            CALL GRTEXT (DX, DY, ' = DEAD NODE')
         END IF
      END IF

C   --Label the node sets

      DY = DYSETS
      IF (NNPS .GT. 0) THEN
         IF (DY .LT. DBOT) GOTO 150
         IF (GRABRT ()) RETURN 1
         CALL GRTEXT (DLEGND(KLFT), DY, 'NODE SETS:')
         DY = DY - CHLSIZ
         DO 130 IX = 1, NNPS
            IF (DY .LT. DBOT) GOTO 150
            INPS = ISCNPS(IX)
            IF (GRABRT ()) RETURN 1
            CALL GRCOLR (INPS)
#if NeedsDoubleEscape
            CALL PLTXTS (DLEGND(KLFT), DY, '\\CI')
#else
            CALL PLTXTS (DLEGND(KLFT), DY, '\CI')
#endif
            CALL PLTXSE (DX, RDUM)
            CALL INTSTR (1, 0, IDNPS(INPS), STRING, LSTR)
            CALL GRTEXT (DX, DY, ' = ID ' // STRING(:LSTR))
            DY = DY - CHLSIZ
  130    CONTINUE
         DY = DY - CHLSIZ
      END IF
      CALL UGRCOL (0, BLKCOL)

C   --Label the side sets

      IF (NESS .GT. 0) THEN
         IF (DY .LT. DBOT) GOTO 150
         IF (GRABRT ()) RETURN 1
         CALL GRTEXT (DLEGND(KLFT), DY, 'SIDE SETS:')
         DY = DY - CHLSIZ
         DO 140 IX = 1, NESS
            IF (DY .LT. DBOT) GOTO 150
            IESS = ISCESS(IX)
            IF (GRABRT ()) RETURN 1
            CALL GRCOLR (IESS)
            BOXSIZ = 0.8 * CHLSIZ
            DX = DLEGND(KLFT) + BOXSIZ
            CALL GRBOX (' ', DLEGND(KLFT), DX, DY, DY+BOXSIZ)
            CALL INTSTR (1, 0, IDESS(IESS), STRING, LSTR)
            CALL GRTEXT (DX, DY, ' = ID ' // STRING(:LSTR))
            DY = DY - CHLSIZ
  140    CONTINUE
         DY = DY - CHLSIZ
      END IF
      CALL UGRCOL (0, BLKCOL)

  150 CONTINUE
      CALL UGRCOL (0, BLKCOL)

C   --Plot time

      IF ((NVWTIM .GT. 0)  .AND.  (LIDSP(0) .GT. 0)) THEN

         IF (GRABRT ()) RETURN 1

         DY = EYTIME

C           Determine the number of variables whose value will be
C           displayed in the plot legend.

         NUMDSP = LIDSP(0)
         IF (NVWTIM .EQ. 4) THEN
            NUMDSP = MIN (NUMDSP, 1)
         ELSE IF (NVWTIM .EQ. 3) THEN
            NUMDSP = MIN (NUMDSP, 1)
         ELSE IF (NVWTIM .EQ. 1) THEN
            NUMDSP = MIN (NUMDSP, 4)
         ELSE
            NUMDSP = MIN (NUMDSP, 2)
         ENDIF

         DO 170 I = NVWTIM, 1, -1
            FIRSTH = .TRUE.
            FIRSTG = .TRUE.
            DO 160 J = NUMDSP, 1, -1
               IF (LIDSP(J) .EQ. 0) THEN
                  DISPVL = 'TIME'
                  VARVAL = TIMES(IPTIMS(I))
               ELSE IF (LIDSP(J) .GT. 0)THEN
                  IF (FIRSTH) THEN
                     CALL DBVIX_BL ('H', 1, IXHV)
                     CALL MDRSRV ('SCRH', KVARH, NVARHI)
                     CALL MDSTAT (NERR, MEM)
                     IF (NERR .GT. 0)GOTO 180
                     CALL GETVAR (A, IXHV, 0, IPTIMS(I), NVARHI,
     &                  A(KVARH))
                     FIRSTH = .FALSE.
                  ENDIF
                  CALL DBVIX_BL ('H', LIDSP(J), IDVAR)
                  DISPVL = NAMES(IDVAR)
                  VARVAL = A(KVARH+LIDSP(J)-1)
               ELSE
                  IF (FIRSTG) THEN
                     CALL DBVIX_BL ('G', 1, IXGV)
                     CALL MDRSRV ('SCRG', KVARG, NVARGL)
                     CALL MDSTAT (NERR, MEM)
                     IF (NERR .GT. 0)GOTO 180
                     CALL GETVAR (A, IXGV, 0, IPTIMS(I), NVARGL,
     &                  A(KVARG))
                     FIRSTG = .FALSE.
                  ENDIF
                  CALL DBVIX_BL ('G', -LIDSP(J), IDVAR)
                  DISPVL = NAMES(IDVAR)
                  VARVAL = A(KVARG-LIDSP(J)-1)
               ENDIF
               CALL NUMSTR (1, 4, VARVAL, RSTR(1), LSTR)
               CALL GRTEXT (DLEGND(KLFT), DY,
     &            DISPVL(:LENSTR(DISPVL)) // '  ' // RSTR(1)(:LSTR))
               DY = DY + CHLSIZ
  160       CONTINUE
            IF (.NOT. FIRSTH) CALL MDDEL ('SCRH')
            IF (.NOT. FIRSTG) CALL MDDEL ('SCRG')
  170    CONTINUE
      ENDIF

C   --Orientation indicator

  180 CONTINUE
      IF (IS3DIM) THEN
         CALL BL_ROTATE (1, 1, ROTMAT, ZERO3,
     &      0.0, 0.0, 0.0, XPT(0), YPT(0), ZPT(0))
         CALL BL_ROTATE (1, 1, ROTMAT, ZERO3,
     &      1.0, 0.0, 0.0, XPT(1), YPT(1), ZPT(1))
         CALL BL_ROTATE (1, 1, ROTMAT, ZERO3,
     &      0.0, 1.0, 0.0, XPT(2), YPT(2), ZPT(2))
         CALL BL_ROTATE (1, 1, ROTMAT, ZERO3,
     &      0.0, 0.0, 1.0, XPT(3), YPT(3), ZPT(3))
         X(KLFT) = MIN (XPT(0), XPT(1), XPT(2), XPT(3))
         X(KRGT) = MAX (XPT(0), XPT(1), XPT(2), XPT(3))
         X(KBOT) = MIN (YPT(0), YPT(1), YPT(2), YPT(3))
         X(KTOP) = MAX (YPT(0), YPT(1), YPT(2), YPT(3))
         CALL SQRLIM (X, X)
         DX = DLEGND(KLFT) + 2*CHLSIZ
         EX = DX + (DYXYZ - EYXYZ)
         LDUM = MPVIEW (DX, EX, EYXYZ, DYXYZ)
         LDUM = MPORT2 (X(KLFT), X(KRGT), X(KBOT), X(KTOP))
         CALL MP2PT (1, XPT(0), YPT(0), XPT(0), YPT(0), IDUM)
         DO 190 I = 1, 3
            CALL MP2PT (1, XPT(I), YPT(I), XPT(I), YPT(I), IDUM)
            CALL GRTEXC (XPT(I), YPT(I)-0.5*CHLSIZ, NAMECO(I))
            X(1) = 0
            X(2) = 0
            X(3) = 0
            X(I) = .9
            CALL BL_ROTATE (1, 1, ROTMAT, ZERO3,
     &         X(1), X(2), X(3), XPT(I), YPT(I), ZPT(I))
            CALL MP2PT (1, XPT(I), YPT(I), XPT(I), YPT(I), IDUM)
            CALL PLTARR (XPT(0), YPT(0), XPT(I), YPT(I), 0.5, .01)
  190    CONTINUE
      END IF

  200 CONTINUE

C *** Axes ***

C   --Reset coordinate system

      CALL MPINIT

      XISCOP = (.NOT. XISSYM) .AND. (MSHDEF(1) .NE. 'NONE')
      YISCOP = (.NOT. YISSYM) .AND. (MSHDEF(4) .NE. 'NONE')

      IF (.NOT. DOAXIS(1)) THEN
         LDUM = PLTSTG1 (KXNUMS, 0.0)
         LDUM = PLTSTG1 (KYNUMS, 0.0)
         LDUM = PLTSTG1 (KTICSZ, 0.0)
      END IF

C   --Compute tick-interval parameters

      IF (XISSYM) THEN
         WAVIEW(KLFT) = WVIEW(KLFT,2) - (WVIEW(KRGT,2) - WVIEW(KLFT,2))
         DAVIEW(KLFT) = DVIEW(KLFT,2) - (DVIEW(KRGT,2) - DVIEW(KLFT,2))
      ELSE
         WAVIEW(KLFT) = WVIEW(KLFT,2)
         DAVIEW(KLFT) = DVIEW(KLFT,2)
      END IF
      WAVIEW(KRGT) = WVIEW(KRGT,2)
      DAVIEW(KRGT) = DVIEW(KRGT,2)
      IF (YISSYM) THEN
         WAVIEW(KBOT) = WVIEW(KBOT,2) - (WVIEW(KTOP,2) - WVIEW(KBOT,2))
         DAVIEW(KBOT) = DVIEW(KBOT,2) - (DVIEW(KTOP,2) - DVIEW(KBOT,2))
      ELSE
         WAVIEW(KBOT) = WVIEW(KBOT,2)
         DAVIEW(KBOT) = DVIEW(KBOT,2)
      END IF
      WAVIEW(KTOP) = WVIEW(KTOP,2)
      DAVIEW(KTOP) = DVIEW(KTOP,2)

      WXATIC = TICMSH
      WYATIC = TICMSH
      CALL GRAPAR (.TRUE., WAVIEW, DAVIEW,
     &   WXALAB, WYALAB, WXAEND, WYAEND, WXATIC, WYATIC)
      WXAMIN = WAVIEW(KLFT)
      WXAMAX = WAVIEW(KRGT)
      WYAMIN = WAVIEW(KBOT)
      WYAMAX = WAVIEW(KTOP)

C   --Adjust for non-symmetric views where both min and max will be labeled,
C   --causing overlap; move first label in 1/2 tick

      IF (XISCOP) THEN
         IF (MAX (WXALAB - WXAMIN, WXAMAX - WXAEND) .LE. (0.2 * WXATIC))
     &      WXALAB = WXALAB + 0.5*WXATIC
      END IF
      IF (YISCOP) THEN
         IF (MAX (WYALAB - WYAMIN, WYAMAX - WYAEND) .LE. (0.2 * WYATIC))
     &      WYALAB = WYALAB + 0.5*WYATIC
      END IF

C   --Set up axis length

      DXALEN = DAVIEW(KRGT) - DAVIEW(KLFT)
      DYALEN = DAVIEW(KTOP) - DAVIEW(KBOT)

      LDUM = PLTSTG1 (KXLENG, DXALEN)
      LDUM = PLTSTG1 (KYLENG, DYALEN)

C   --Set up axis minimum, maximum, and tick intervals

      BUF(1) = 4
      BUF(2) = WXAMIN
      BUF(3) = WXALAB
      BUF(4) = WXAMAX
      BUF(5) = WXATIC
      BUF(6) = 0.0
      BUF(7) = WYAMIN
      BUF(8) = WYALAB
      BUF(9) = WYAMAX
      BUF(10)= WYATIC
      BUF(11)= 0.0
      LDUM = PLTSTG (KSCALE, BUF)

      LDUM = PLTGTG (KXNUMS, SZNUM)
      LDUM = PLTGTG (KXLABS, SZLAB)

      DO 210 IVW = 1, NDEFVW (.TRUE.)
         IVIEW = IXVW (.TRUE., IVW)

         IF ((IVIEW .EQ. 2)
     &      .OR. (XISCOP .AND. YISCOP)
     &      .OR. (XISCOP .AND. TOPVW (IVIEW))
     &      .OR. (YISCOP .AND. RGTVW (IVIEW))) THEN

            IF (GRABRT ()) RETURN 1

C         --Set up axes labeling

            IF (YISCOP .AND. TOPVW(IVIEW)) THEN
               LDUM = PLTSTG1 (KXNUMS, 0.0)
               LDUM = PLTSTG1 (KXLABS, 0.0)
            ELSE
               IF (DOAXIS(1)) LDUM = PLTSTG1 (KXNUMS, SZNUM)
               LDUM = PLTSTG1 (KXLABS, SZLAB)
            END IF
            IF (XISCOP .AND. RGTVW(IVIEW)) THEN
               LDUM = PLTSTG1 (KYNUMS, 0.0)
               LDUM = PLTSTG1 (KYLABS, 0.0)
            ELSE
               IF (DOAXIS(1)) LDUM = PLTSTG1 (KYNUMS, SZNUM)
               LDUM = PLTSTG1 (KYLABS, SZLAB)
            END IF

C         --Set the axis start

            IF (XISSYM .OR. (.NOT. RGTVW(IVIEW))) THEN
               DXASTA = DVIEW(KLFT,1)
            ELSE
               DXASTA = DVIEW(KLFT,2)
            END IF
            IF (YISSYM .OR. (.NOT. TOPVW(IVIEW))) THEN
               DYASTA = DVIEW(KBOT,4)
            ELSE
               DYASTA = DVIEW(KBOT,2)
            END IF

            LDUM = PLTSTG1 (KXORIG, DXASTA)
            LDUM = PLTSTG1 (KYORIG, DYASTA)

C         --Draw the axes

            if (dobox) then
              IF (.NOT. IS3DIM) THEN
                CALL PLTGPH (0., 0., 0, NAMECO(1), ' ', NAMECO(2), ' ')
              ELSE
                CALL PLTGPH (0., 0., 0, ' ', ' ', ' ', ' ')
              END IF
            END IF
          end if
  210 CONTINUE

C   --Flush buffer, so label is complete at this point
      CALL PLTFLU

      RETURN

  220 CONTINUE
      RETURN 1
      END
