C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE GRAPAR (XYSAME, WVIEW, DVIEW,
     &   WXALAB, WYALAB, WXAEND, WYAEND, WXATIC, WYATIC)
C=======================================================================

C   --*** GRAPAR *** (GRPLIB) Determine axes parameters
C   --   Written by Amy Gilkey - revised 02/20/87
C   --
C   --GRAPAR sets up the axes parameters.  It determines "good" numbers
C   --for the axes numbers and sets the exponent and number of digits.
C   --It also sets the axis number and label number sizes.
C   --
C   --Parameters:
C   --   XYSAME - IN - true iff the X and Y axis have the same scale;
C   --      i.e., they are the same type although the values may differ
C   --   WVIEW - IN - the window corners (left, right, bottom, top)
C   --      in window (user) coordinates
C   --   DVIEW - IN - the window corners (left, right, bottom, top)
C   --      in device coordinates
C   --   WXALAB, WYALAB - OUT - the starting tick-mark of the X and Y axis
C   --   WXAEND, WXAEND - OUT - the ending tick-mark of the X and Y axis
C   --   WXATIC, WYATIC - IN/OUT - the X and Y axis tick-mark interval;
C   --      set only if equal zero or invalid

C   --Routines Called:
C   --   PLTGTG - (PLTLIB) Get graph parameter (see PLTSTG)
C   --   PLTSTG - (PLTLIB) Set graph parameter
C   --      1, 2 = (KXORIG, KYORIG) X, Y axis origin location
C   --      3, 4 = (KXLENG, KYLENG) X, Y axis length
C   --      11 = (KSCALE) axes parameters (see documentation)
C   --      22, 47 = (KXNUMS, KYNUMS) X, Y axis number size
C   --      23, 48 = (KXLABS, KYLABS) X, Y axis label size
C   --   PLTGTT - (PLTLIB) Get text parameter
C   --      2 = (KSCHSZ) software character size
C   --   GRAEXP - (GRPLIB) Set axis exponent

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4)

      PARAMETER (KSCHSZ=2)
      PARAMETER (KXORIG=1, KYORIG=2, KXLENG=3, KYLENG=4, KSCALE=11)
      PARAMETER (KXNUMS=22, KXLABS=23, KYNUMS=47, KYLABS=48)

      LOGICAL XYSAME
      REAL WVIEW(KTOP), DVIEW(KTOP)
      REAL WXALAB, WYALAB
      REAL WXAEND, WYAEND
      REAL WXATIC, WYATIC

      LOGICAL LDUM, PLTGTG, PLTSTG, PLTGTT
      LOGICAL NUMBX, NUMBY
      CHARACTER*8 NSTR
      REAL BUF(4)

      WXAMIN = WVIEW(KLFT)
      WXAMAX = WVIEW(KRGT)
      WYAMIN = WVIEW(KBOT)
      WYAMAX = WVIEW(KTOP)

      DXALEN = DVIEW(KRGT) - DVIEW(KLFT)
      DYALEN = DVIEW(KTOP) - DVIEW(KBOT)
      WXALEN = ABS (WXAMAX - WXAMIN)
      WYALEN = ABS (WYAMAX - WYAMIN)

C   --Set tick-interval

      TICFAC = (1.00 / 12) * 2
C      --TICFAC is divided by the axis length to get the length of 2 times
C      --a "nice" tick interval (in device coordinates); a "nice" tick interval
C      --is 12 ticks across the screen (horizontally)

      IF (WXATIC .NE. 0.0) THEN
         XN = WXALEN / WXATIC
         IF ((XN .LT. 1) .OR. (XN .GT. 20)) WXATIC = 0.0
      END IF
      IF (WXATIC .EQ. 0.0) THEN
         XTIC = WXALEN * (TICFAC / DXALEN)
         WRITE (NSTR, 10000) XTIC
         READ (NSTR, 10000) WXATIC
         WXATIC = .5 * WXATIC
      END IF

      IF (WYATIC .NE. 0.0) THEN
         XN = WYALEN / WYATIC
         IF ((XN .LT. 1) .OR. (XN .GT. 20)) WYATIC = 0.0
      END IF
      IF (WYATIC .EQ. 0.0) THEN
         YTIC = WYALEN * (TICFAC / DYALEN)
         WRITE (NSTR, 10000) YTIC
         READ (NSTR, 10000) WYATIC
         WYATIC = .5 * WYATIC
      END IF

      IF (XYSAME) THEN
         IF (WXATIC .LT. WYATIC) THEN
            WXATIC = WYATIC
         ELSE
            WYATIC = WXATIC
         END IF
      END IF

C   --Set axes starting and ending tick-marks

      WXALAB = WXATIC * AINT (WXAMIN/WXATIC)
      IF (WXALAB .LT. WXAMIN) WXALAB = WXALAB + WXATIC
      WXAEND = WXALAB + WXATIC * AINT ((WXAMAX-WXAMIN) / WXATIC)
      IF (WXAEND .GT. WXAMAX) WXAEND = WXAEND - WXATIC

      WYALAB = WYATIC * AINT (WYAMIN/WYATIC)
      IF (WYALAB .LT. WYAMIN) WYALAB = WYALAB + WYATIC
      WYAEND = WYALAB + WYATIC * AINT ((WYAMAX-WYAMIN) / WYATIC)
      IF (WYAEND .GT. WYAMAX) WYAEND = WYAEND - WYATIC

C   --Set axis numbering (exponent and number of decimal digits)

      LDUM = PLTGTG (KXNUMS, SZXNUM)
      LDUM = PLTGTG (KYNUMS, SZYNUM)
      NUMBX = (SZXNUM .GT. 0.0)
      NUMBY = (SZYNUM .GT. 0.0)

      BUF(1) = WXALAB
      BUF(2) = WXAEND
      BUF(3) = WYALAB
      BUF(4) = WYAEND

      IF (XYSAME .AND. NUMBX .AND. NUMBY) THEN
         CALL GRAEXP (' ', 4, BUF, WXATIC)
      ELSE
         IF (NUMBX) THEN
            CALL GRAEXP ('X', 2, BUF(1), WXATIC)
         ELSE
            CALL GRAEXP ('X', 0, BUF(1), WXATIC)
         END IF
         IF (NUMBY) THEN
            CALL GRAEXP ('Y', 2, BUF(3), WYATIC)
         ELSE
            CALL GRAEXP ('Y', 0, BUF(3), WYATIC)
         END IF
      END IF

C   --Set up label/numbering size

      LDUM = PLTGTG (KXNUMS, SZNUM)
      IF (SZNUM .LE. 0.0) LDUM = PLTGTG (KYNUMS, SZNUM)
      LDUM = PLTGTG (KXLABS, SZLAB)
      IF (SZLAB .LE. 0.0) LDUM = PLTGTG (KYLABS, SZLAB)

      IF ((SZNUM .GT. 0.0) .OR. (SZLAB .GT. 0.0)) THEN
C      --The size of the label/numbering is determined by the character size
C      --and the axis lengths
         LDUM = PLTGTT (KSCHSZ, VCS)
         R = VCS * 50.0 / (0.5 * (DXALEN + DYALEN))
         SZNUM = SZNUM * R
         SZLAB = SZLAB * R*(4.0/5.0)
C         --Factor by 4/5 to make label and numbering the same size

         LDUM = PLTGTG (KXNUMS, X)
         IF (X .GT. 0.0) LDUM = PLTSTG (KXNUMS, SZNUM)
         LDUM = PLTGTG (KXLABS, X)
         IF (X .GT. 0.0) LDUM = PLTSTG (KXLABS, SZLAB)
         LDUM = PLTGTG (KYNUMS, X)
         IF (X .GT. 0.0) LDUM = PLTSTG (KYNUMS, SZNUM)
         LDUM = PLTGTG (KYLABS, X)
         IF (X .GT. 0.0) LDUM = PLTSTG (KYLABS, SZLAB)
      END IF

      RETURN
10000  FORMAT (E8.1)
      END
