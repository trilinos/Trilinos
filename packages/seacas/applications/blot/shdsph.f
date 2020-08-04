C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SHDSPH (LENE, LINK, NUMLNK, NUMATR, XN, YN, ZN, ATRIB,
     *  BLKCOL, IDELB, ISPSOR, RAD,
     *  IELBST, ISPBLK, SHDCOL, ISHDCL, HIDENP, *)
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
      INTEGER ISPSOR(NUMEL)
      REAL RAD(NUMEL)
      INTEGER IELBST(NELBLK), ISPBLK(NUMEL)
C ... SHDCOL(1, *) = Red   Component
C     SHDCOL(2, *) = Green Component
C     SHDCOL(3, *) = Blue  Component
C     SHDCOL(4-7,*) - Future Use
      REAL SHDCOL(7,NELBLK)
C ... ISHDCL(1, *) = -1 if color not set, >0 if color set
C     ISHDCL(2, *) = Number of colors to use for this block (SET if 0)
C     ISHDCL(3, *) = Starting location in color map (SET)
      INTEGER ISHDCL(3,NELBLK)
      LOGICAL HIDENP(*)

      REAL XSI(4),YSI(4)
      LOGICAL FIRST

      include 'icrnbw.blk'
      include 'light.blk'
      include 'sphere.blk'
      DATA FIRST /.TRUE./

      XMIN = ZMMESH(KLFT)
      XMAX = ZMMESH(KRGT)
      YMIN = ZMMESH(KBOT)
      YMAX = ZMMESH(KTOP)

      if (.not. is3dim) return

      IF ( FIRST ) THEN
        FIRST = .FALSE.
      END IF

C ... Calculate surface normals
C ... On a unit sphere, normals are simply x,y,z coords of point.
      do 30 ip=1, npoly
        i = icon(1,ip)
        j = icon(2,ip)
        k = icon(3,ip)
        l = icon(4,ip)
        XMAG =  (xpts(i) + xpts(j) + xpts(k) + xpts(l)) / 4.0
        YMAG =  (ypts(i) + ypts(j) + ypts(k) + ypts(l)) / 4.0
        ZMAG =  1.0 - SQRT(XMAG**2 + YMAG**2)

C ... Determine dot product of surface normal and light vector.
        SMAG = 0.0
        SHDMX = 0.0
        DO 20 J = 1, NLIT
          SMAG = SMAG + max(0.0, (XMAG * LITE(5,J) +
     *      YMAG * LITE(6,J) + ZMAG * LITE(7,J)) * LITE(8,J))
 20     CONTINUE
        SHADE(IP) = SMAG
        SHDMX = MAX(SHDMX, SMAG)
 30   CONTINUE
      SHDMX = MIN(1.0, 1.0 / SHDMX)
      DO 40 ip = 1, npoly
        SHADE(IP) = SHADE(IP) * SHDMX
 40   continue

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
      call indexx (zn, ispsor, numsph, .FALSE.)

C ... Set colors for all blocks
      call setcol (nelblk, shdcol, ishdcl, ielbst, blkcol, idelb)

C ... Plotting of spheres starts here.
      do 130 iel = 1, numsph
        node   = ispsor(iel)
        xc     = xn(node)
        yc     = yn(node)
        rd     = rad(node)
        ielb   = ispblk(node)
        mincol = ishdcl(3,ielb)
        ncol   = ishdcl(2,ielb)

        do 120 ip = 1, npoly
          do 118 in = 1, 4
            xsi(in) = xc + rd * xpts(icon(in,ip))
            ysi(in) = yc + rd * ypts(icon(in,ip))
 118      continue
          smag = shade(ip)
          ICOL = MINCOL +
     *      min(NCOL-1, max(1, NINT(SHADE(IP) * DBLE(NCOL))))
          CALL GRCOLR (ICOL)
          CALL MPD2PG (4, XSI, YSI, 'S')
 120    continue
 130  continue
      call pltflu

      CALL GRCOLU('STANDARD')

      RETURN
      END
