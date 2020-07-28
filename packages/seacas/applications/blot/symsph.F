C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SYMSPH (LENE, LINK, NUMLNK, NUMATR, XN, YN, ZN, ATRIB,
     *  BLKCOL, FFLAG, IDELB, VARNP, MODDET, ISPSOR, RAD,
     *  IELBST, ISPBLK, HIDENP, *)
C=======================================================================

C   --*** SYMSPH *** (DETOUR) Plot element spheres
C   --   Written by Lee Taylor, 07/12/88
C   --   Amy Gilkey, 10/12/87
C   --   D. P. Flanagan, 12/08/83
C   --   Modified version 1.1a  November 1990  - R.J. Meyers
C   --           added color coded sphere capability
C   --
C   --SYMSPH drives the user symbol interface for element variables.
C   --It processes each element by element block, computing scaling factors
C   --and element information, then calling the user symbol routine.
C   --Only elements of selected element blocks are drawn.
C   --
C   --Parameters:
C   --   LENE - IN - the cumulative element counts by element block
C   --   LINK - IN - the master connectivity array
C   --   NUMLNK - IN - the number of nodes per element by block
C   --   NUMATR - IN - the number of attributes per element by block
C   --   XN, YN, ZN - IN - the nodal coordinates
C   --   ATRIB - IN - the element attributes array
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
C   --   FFLAG - IN - = 'NOFILL' if the polygons making the sphere
C   --                  are not to be filled.
C   --                = 'FILL' if they are to be filled.
C   --   VARNP - IN - function value array for painted spheres
C   --   MODDET - IN -  the DETOUR mode to check for contour mode
C   --   * - return statement if the cancel function is active
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)
      include 'mshlim.blk'
      include 'd3nums.blk'
      include 'dbnums.blk'
      include 'sphele.blk'

      INTEGER LENE(0:NELBLK),LINK(*)
      INTEGER NUMLNK(NELBLK),NUMATR(NELBLK)
      REAL XN(*), YN(*), ZN(*), ATRIB(*)
      INTEGER BLKCOL(0:NELBLK)
      INTEGER IDELB(*)
      CHARACTER*6 FFLAG
      REAL VARNP(*)
      CHARACTER*8 MODDET
      INTEGER ISPSOR(NUMEL)
      REAL RAD(NUMNP)
      INTEGER IELBST(NELBLK), ISPBLK(NUMNP)
      LOGICAL HIDENP(*)

      LOGICAL CONTOR

      include 'cntr.blk'
      DATA NPTSPH /0/

      XMIN = ZMMESH(KLFT)
      XMAX = ZMMESH(KRGT)
      YMIN = ZMMESH(KBOT)
      YMAX = ZMMESH(KTOP)

      NPTSPH = MIN(ABS(SPHPLT), NPTSPX)
      ANGINC = 8. * ATAN(1.) / NPTSPH
      ANGLE = 0.
      DO 100 N = 0, NPTSPH
        CS(N) = COS(ANGLE)
        SN(N) = SIN(ANGLE)
        ANGLE = ANGLE + ANGINC
 100  CONTINUE

      iel = 0
      IRAD = 1
      ILNK = 1
      do 104 ielb = 1, nelblk
        nel = lene(ielb) - lene(ielb-1)
        if (numlnk(ielb) .eq. 1 .AND. IELBST(IELB) .GT. 0) then
          DO 102 ISPH = ILNK, ILNK+NEL-1
            NODE = LINK(ISPH)
            IF (IS3DIM) THEN
               IF (HIDENP(NODE)) GOTO 102
            END IF
            if (numatr(ielb) .ge. 1) then
              RAD(NODE) = DEFRAD * ATRIB(IRAD)
            else
              RAD(NODE) = DEFRAD
            end if
            X = XN(NODE)
            Y = YN(NODE)
            R = RAD(NODE)
C ... Determine if portion of sphere is in viewing window
            if (x+r .ge. xmin .and. x-r .le. xmax .and.
     *        y+r .ge. ymin .and. y-r .le. ymax) then
              iel = iel + 1
              ISPSOR(IEL) = NODE
              ISPBLK(NODE) = ielb
            end if
            IRAD = IRAD + NUMATR(IELB)
 102      continue
        ELSE
          IRAD = IRAD + NEL * NUMATR(IELB)
        END IF
        ILNK = ILNK + NEL * NUMLNK(IELB)
 104  continue
      numsph = iel
C ... Sort elements from smallest Z coord to largest.
C     Sort is based on element center.
      if (is3dim) then
        call indexx (zn, ispsor, numsph, .FALSE.)
      end if

C CHECK TO SEE IF CONTOURING IS ON

      CONTOR = (MODDET .EQ. 'CONTOUR')

C IF CONTOURING MODE, SET ALTERNATE COLORS

      IF(CONTOR) THEN
         CALL GRCOLU('ALTERNATE')
      END IF

      ieblst = 0
      do 130 iel = 1, numsph
        node = ispsor(iel)
        xc = xn(node)
        yc = yn(node)
        rd = rad(node)
        ielb = ispblk(node)
        IF (CONTOR) THEN
          FUNVAL = VARNP(NODE)
          ICOL = INT( (FUNVAL-CMIN)/(CMAX-CMIN)*NCNTR + 1 )
          IF (NOCMIN .AND. (ICOL .LT. 1)) ICOL = 1
          IF (NOCMAX .AND. (ICOL .GT. NCNTR)) ICOL = NCNTR
          IF(ICOL .LT. 1) ICOL = -1
          IF(ICOL .GT. NCNTR) ICOL = -1
          CALL GRCOLR(ICOL)
        ELSE IF (IEBLST .NE. IELB .OR. FFLAG .EQ. 'FILL') THEN
          IEBLST = IELB
          CALL UGRCOL (IELB, BLKCOL)
        END IF

C            --Draw circle
        IF (NPTSPH .GT. 1) THEN
          DO 110 N = 0, NPTSPH
            XS(N) = XC + RD * CS(N)
            YS(N) = YC + RD * SN(N)
 110      CONTINUE

C ... Fill the sphere if 'FILL' specified OR if in CONTOR mode AND
C    the function value is within the the contour range.
C    If  we fill with black, then we obscure the spheres that are
C    behind this sphere.
          IF (NPTSPH .GT. 2 .AND. (FFLAG .EQ. 'FILL' .OR.
     *      (CONTOR .AND. ICOL .NE. -1))) THEN
            CALL MPD2PG (NPTSPH+1, XS, YS, 'S')
            CALL GRCOLR( -1 )
            CALL MPD2VC( NPTSPH, XS, YS, XS(1), YS(1) )
          ELSE
            CALL MPD2VC( NPTSPH, XS, YS, XS(1), YS(1) )
          ENDIF
        ELSE
#if NeedsDoubleEscape
          CALL MPD2SY (1, XC, YC, '\\CCS')
#else
          CALL MPD2SY (1, XC, YC, '\\CCS')
#endif
        END IF
 130  continue
      call pltflu

C IF CONTOURING MODE, SET STANDARD COLORS

      IF(CONTOR) THEN
         CALL GRCOLU('STANDARD')
      END IF

      RETURN

      END
