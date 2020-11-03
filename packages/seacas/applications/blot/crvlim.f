C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CRVLIM (AXIS, TIMPLT, MAXPTS, NPTS, NSPVAR, NEPVAR,
     &   PLTVAL)
C=======================================================================

C   --*** CRVLIM *** (TPLOT) Calculate min/max value for plot data
C   --   Written by Amy Gilkey - revised 11/06/87
C   --
C   --CRVLIM calculates the minimum and maximum of the plot data.
C   --
C   --Parameters:
C   --   AXIS - IN - the axis to scale ('X' or 'Y')
C   --   TIMPLT - IN - true if time plot Y data versus X-Y data
C   --   MAXPTS - IN - the maximum number of points on a curve (PLTVAL length)
C   --   NPTS - IN - the number of points on each curve
C   --   NSPVAR, NEPVAR - IN - the starting and ending plot variable
C   --      indices for min/max calculation (PLTVAL index, /TPVARS/ index)
C   --   PLTVAL - IN - the data array
C   --
C   --Common Variables:
C   --   Sets XMIN, XMAX, YMIN, YMAX of /XYLIM/

      include 'xylim.blk'

      CHARACTER AXIS
      LOGICAL TIMPLT
      INTEGER NPTS(*)
      REAL PLTVAL(MAXPTS,*)

      IF ((.NOT .TIMPLT) .AND. (AXIS .EQ. 'X')) THEN
         XMIN =  1.0E+30
         XMAX = -1.0E+30

         N = NSPVAR
  100    CONTINUE
         IF (N .LE. NEPVAR) THEN
            DO 110 I = 1, NPTS(N)
               XMIN = MIN (XMIN, PLTVAL(I,N))
               XMAX = MAX (XMAX, PLTVAL(I,N))
  110       CONTINUE
            N = N + 2
            GOTO 100
         END IF
      END IF

      IF (AXIS .EQ. 'Y') THEN
         YMIN =  1.0E+30
         YMAX = -1.0E+30

         N = NSPVAR
  120    CONTINUE
         IF (N .LE. NEPVAR) THEN
            IF (.NOT. TIMPLT) N = N + 1
            DO 130 I = 1, NPTS(N)
               YMIN = MIN (YMIN, PLTVAL(I,N))
               YMAX = MAX (YMAX, PLTVAL(I,N))
  130       CONTINUE
            N = N + 1
            GOTO 120
         END IF
      END IF

C ... Check for NaN (Not a Number).
C     NaN is defined to be not equal to any other number including itself
      if (ymin .ne. ymin .or. ymax .ne. ymax .or.
     *    xmin .ne. xmin .or. xmax .ne. xmax) then
        call prterr('WARNING',
     *    'The variable contains "NaN" (Not a Number) values.'//
     *    ' Check data.')
        if (xmin .ne. xmin) xmin = -1.0e30
        if (xmax .ne. xmax) xmax =  1.0e30
        if (ymin .ne. ymin) ymin = -1.0e30
        if (ymax .ne. ymax) ymax =  1.0e30
      end if

      RETURN
      END
