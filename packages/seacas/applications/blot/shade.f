C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SHADE (LENF, NLNKF, LINKF, NXFAC, IXFAC,
     &   XN, YN, ZN, IELBST, BLKCOL, IDELB, SHDCOL, ISHDCL, IHIDOP, *)
C=======================================================================

C   --*** SOLID *** (DETOUR) Paint solid mesh (by index)
C   --   Modified by John Glick - 11/29/88
C   --   Written by Amy Gilkey - revised 10/27/87
C   --
C   --SOLID paints the mesh in the color of each elements element block.
C   --
C   --Parameters:
C   --   LENF - IN - the cumulative face counts by element block
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF - IN - the connectivity for all faces
C   --   NXFAC - IN - the number of ordered faces (if DOIXF)
C   --   IXFAC - IN - the indices of the ordered faces (if DOIXF)
C   --   XN, YN, ZN - IN - the nodal coordinates
C   --   IELBST - IN - the element block status (>0 if selected)
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
C   --   Uses NUMEL, NELBLK of /DBNUMS/
C   --   Uses IS3DIM of /D3NUMS/

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)
      include 'mshlim.blk'
      include 'icrnbw.blk'
      include 'light.blk'
      include 'dbnums.blk'
      include 'd3nums.blk'

      INTEGER LENF(0:NELBLK)
      INTEGER NLNKF(NELBLK)
      INTEGER LINKF(*)
      INTEGER IXFAC(*)
      REAL XN(*), YN(*), ZN(*)
      INTEGER IELBST(NELBLK)
      INTEGER BLKCOL(0:NELBLK)
      INTEGER IDELB(*)
C ... SHDCOL(1, *) = Red   Component
C     SHDCOL(2, *) = Green Component
C     SHDCOL(3, *) = Blue  Component
C     SHDCOL(4-7,*) - Future Use
      REAL SHDCOL(7,NELBLK)
C ... ISHDCL(1, *) = -1 if color not set, >0 if color set
C     ISHDCL(2, *) = Number of colors to use for this block (SET if 0)
C     ISHDCL(3, *) = Starting location in color map (SET)
      INTEGER ISHDCL(3,NELBLK)

C ... AMBIENT controls the range of hue values. The hue is scaled by:
C         HUE = (1 - AMBIENT) * HUE + AMBIENT
C ... It is somewhat similar to the ambient intensity.

      if (.not. is3dim) return

C ... Determine range of plot window
      XMIN = ZMMESH(KLFT)
      XMAX = ZMMESH(KRGT)
      YMIN = ZMMESH(KBOT)
      YMAX = ZMMESH(KTOP)

C ... Set colors for all blocks
      call setcol (nelblk, shdcol, ishdcl, ielbst, blkcol, idelb)

C ... Heres where we actually start plotting something.

      call excpus(time1)
      DO 100 IX = 1, NXFAC
        IFAC = IXFAC(IX)
        IELB = 0
        IXL = IDBLNK (IELB, IFAC, LENF, NLNKF)

        NNPF = NLNKF(IELB)
        CALL SHADEN (NNPF, LINKF(IXL), XN, YN, ZN,
     *    ishdcl(2,ielb), LITE, NLIT, ishdcl(3,ielb),
     *    shdcol(4,ielb), shdcol(5,ielb), shdcol(6,ielb),
     *    xmin, ymin, xmax, ymax)
 100  CONTINUE
      call excpus(time2)
      tottim = time2 - time1
      if (tottim .gt. 0.0) then
        write (*,900) nxfac, time2-time1, nxfac/(time2-time1)
      end if
 900  FORMAT (' Processed ', i6, ' faces in ', 1pe10.3, ' seconds. ',
     *  '(',1pe10.3,' faces/second)',/)

      CALL PLTFLU

      RETURN
      END

