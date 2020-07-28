C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE GRNCRV (LABSID, ICURVE, NPTS, XPTS, YPTS, DRAWLN)
C=======================================================================

C   --*** GRNCRV *** (GRPLIB) Label Curve (PLT)
C   --   Written by Amy Gilkey - revised 10/07/87
C   --
C   --GRNCRV finds the device coordinates of the first, last or middle
C   --point of a curve within the plot window and labels the curve with
C   --the line number.
C   --
C   --Parameters:
C   --   LABSID - IN - label curve on "FIRST", "MIDDLE" or "LAST" point
C   --   ICURVE - IN - the curve number
C   --   NPTS - IN - the number of points in the curve
C   --   XPTS, YPTS - IN - the X and Y points
C   --   DRAWLN - IN - true if vector, false if points

C   --Routines Called:
C   --   MP2VC - (PLTLIB) Get device coordinates of vector
C   --   MP2PT - (PLTLIB) Get device coordinates of point
C   --   PLTGTT - (PLTLIB) Get text parameter
C   --      2 = (KSCHSZ) software character size
C   --   INTSTR - (STRLIB) Convert integer to string

      PARAMETER (KSCHSZ=2)

      CHARACTER*8 LABSID
      INTEGER ICURVE
      INTEGER NPTS
      REAL XPTS(*), YPTS(*)
      LOGICAL DRAWLN

      LOGICAL LDUM, PLTGTT
      CHARACTER*2 STR2

C   --Find the nearest line or point within the window

      IF (LABSID .EQ. 'FIRST') THEN
         ISTART = 1
         IEND = NPTS
         INC = 1
      ELSE IF (LABSID .EQ. 'MIDDLE') THEN
         ISTART = NPTS / 2
         IEND = NPTS
         INC = 1
      ELSE
         ISTART = NPTS
         IEND = 1
         INC = -1
      END IF

      IF (DRAWLN) THEN
         DO 100 I = ISTART, IEND-INC, INC
            MASK = -1
            CALL MP2VC (1, XPTS(I), YPTS(I), XPTS(I+INC), YPTS(I+INC),
     &         DX, DY, DXX, DYY, MASK)
            IF (MOD (MASK,2) .NE. 0) GOTO 120
  100    CONTINUE

      ELSE
         DO 110 I = ISTART, IEND, INC
            MASK = -1
            CALL MP2PT (1, XPTS(I), YPTS(I), DX, DY, MASK)
            IF (MOD (MASK,2) .NE. 0) GOTO 120
  110    CONTINUE
      END IF

      GOTO 130

C   --Write the curve number (location is dependent on the placement)

  120 CONTINUE
C   --Get curve number
      CALL INTSTR (1, 0, ICURVE, STR2, LNUM)

      LDUM = PLTGTT (KSCHSZ, VCS)
      IF (DRAWLN) THEN
         DXO = .25
      ELSE
         DXO = .50
      END IF
      IF (LABSID .EQ. 'FIRST') THEN
         DX = DX - (LNUM+DXO) * VCS
      ELSE
         DX = DX + DXO*VCS
      END IF
      IF (LABSID .EQ. 'MIDDLE') THEN
         DY = DY + .002
      ELSE
         DY = DY - .5*VCS
      END IF

      CALL GRTEXT (DX, DY, STR2)

  130 CONTINUE
      RETURN
      END
