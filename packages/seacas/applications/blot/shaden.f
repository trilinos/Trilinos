C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SHADEN (NLNKF, LINKF1, XN, YN, ZN, NCOL, LITE, NLIT,
     *  MINCOL, KDIFF, KSPEC, SPEXP, WXMIN, WYMIN, WXMAX, WYMAX)
C=======================================================================

C   --*** SOLIDF *** (DETOUR) Paint face
C   --   Written by Amy Gilkey - revised 09/24/85
C   --
C   --SOLIDF paints the a face of the mesh.
C   --
C   --Parameters:
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF1 - IN - the connectivity for the face
C   --   XN, YN, ZN - IN - the nodal coordinates
C   --   NCOL - IN - the number of colors
C   --   XLN, YLN, ZLN - IN - the light unit vector components
      INTEGER LINKF1(NLNKF)
      REAL XN(*), YN(*), ZN(*)
      REAL LITE(8,*)
      REAL KDIFF, KSPEC, SPEXP

      REAL XPTS(20), YPTS(20), ZPTS(20)
      SAVE LSTCOL
      DATA LSTCOL /0/

C... Coordinate System:

C    ^ Y
C    |
C    |
C    |
C    |         X
C    Z--------->

      XMAX = -1.0e30
      YMAX = -1.0e30
      XMIN =  1.0e30
      YMIN =  1.0e30
      DO 100 ILINK = 1, NLNKF
         XPTS(ILINK) = XN(LINKF1(ILINK))
         YPTS(ILINK) = YN(LINKF1(ILINK))
         ZPTS(ILINK) = ZN(LINKF1(ILINK))
         XMAX = MAX(XMAX, XPTS(ILINK))
         XMIN = MIN(XMIN, XPTS(ILINK))
         YMAX = MAX(YMAX, YPTS(ILINK))
         YMIN = MIN(YMIN, YPTS(ILINK))
  100 CONTINUE

C ... Determine if polygon is in window.
      IF (XMAX .LT. WXMIN .OR. XMIN .GT. WXMAX .OR.
     *    YMAX .LT. WYMIN .OR. YMIN .GT. WYMAX) RETURN

C ... Calculate surface normal
      if (nlnkf .eq. 4) then
        XMAG =  (YPTS(3) - YPTS(1)) * (ZPTS(4) - ZPTS(2)) -
     *          (ZPTS(3) - ZPTS(1)) * (YPTS(4) - YPTS(2))
        YMAG =  (ZPTS(3) - ZPTS(1)) * (XPTS(4) - XPTS(2)) -
     *          (XPTS(3) - XPTS(1)) * (ZPTS(4) - ZPTS(2))
        ZMAG =  (XPTS(3) - XPTS(1)) * (YPTS(4) - YPTS(2)) -
     *          (YPTS(3) - YPTS(1)) * (XPTS(4) - XPTS(2))
      else
C ... Newells method, "Procedural Elements for Computer Graphics", p209
        xmag = (ypts(nlnkf) - ypts(1)) * (zpts(nlnkf) - zpts(1))
        ymag = (zpts(nlnkf) - zpts(1)) * (xpts(nlnkf) - xpts(1))
        zmag = (xpts(nlnkf) - xpts(1)) * (ypts(nlnkf) - ypts(1))
        do 110 i=1, nlnkf-1
          j = i + 1
          xmag = xmag + (ypts(i) - ypts(j)) * (zpts(i) - zpts(j))
          ymag = ymag + (zpts(i) - zpts(j)) * (xpts(i) - xpts(j))
          zmag = zmag + (xpts(i) - xpts(j)) * (ypts(i) - ypts(j))
 110    continue
      end if
C ... If polygon facing away from viewer (-Z), return immediately
      if (zmag .lt. 0) return
      VMAG = SQRT(XMAG**2 + YMAG**2 + ZMAG**2)

C ... Normalize surface normal
      XMAG = XMAG / VMAG
      YMAG = YMAG / VMAG
      ZMAG = ZMAG / VMAG

C ... Determine dot product of surface normal and light vector.
      SMAG = 0.0
      DO 120 I = 1, NLIT
         SMAG = SMAG + max(0.0, (XMAG * LITE(5,I) + YMAG * LITE(6,I) +
     *                 ZMAG * LITE(7,I)) * LITE(8,I))
  120 CONTINUE

C ... Calculate Reflection Vector.
      if (SPEXP .GT. 0.0 .and. KSPEC .GT. 0.0) then
        RMAG = 0.0
        DO 160 I = 1, NLIT
          R = xmag * lite(5,i) + ymag*lite(6,i) + zmag*lite(7,i)
          x = 2 * xmag * R - lite(5,i)
          y = 2 * ymag * R - lite(6,i)
          z = 2 * zmag * R - lite(7,i)

C ... Normalize Reflection Vector
          vmag = sqrt(x**2 + y**2 + z**2)
          x = x / vmag
          y = y / vmag
          z = z / vmag

C ... NOTE: Sight vector is (0,0,1),

C ... Calculate dot product of Reflect and Sight
C     - sight = z, therefore, dot product = z
          RMAG = RMAG + max(0.0, Z)**SPEXP * LITE(8,I)
 160    CONTINUE

        SMAG = KDIFF * SMAG + KSPEC * RMAG
      end if

C ... Scale dot product (0 <= |smag| <= 1) to (0 < icol < ncol)
C     If smag <= 0, set icol = 1
C     maximum color = ncol-1
C     minimum color = 1

      if (ncol .eq. 1) then
        icol = mincol + 1
      else
        ICOL = MINCOL + min(NCOL-1, max(1, NINT(SMAG * DBLE(NCOL))))
      end if

      if (lstcol .ne. icol) then
        CALL GRCOLR (ICOL)
        lstcol = icol
      end if
      CALL MPD2PG (NLNKF, XPTS, YPTS, 'S')
      RETURN
      END
