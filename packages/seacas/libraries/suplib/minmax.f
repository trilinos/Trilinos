C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MINMAX (NPTS, PTS, VMIN, VMAX)
C=======================================================================

C   --*** MINMAX *** (ETCLIB) Calculate min/max value
C   --   Written by Amy Gilkey - revised 11/10/87
C   --
C   --MINMAX calculates the minimum and maximum of the data.
C   --
C   --Parameters:
C   --   NPTS - IN - the number of points
C   --   PTS - IN - the points
C   --   VMIN, VMAX - OUT - the maximum and maximum value of the points

      INTEGER NPTS
      REAL PTS(*)
      REAL VMIN, VMAX

      VMIN =  1.0E+30
      VMAX = -1.0E+30
      DO 10 I = 1, NPTS
         VMIN = MIN (VMIN, PTS(I))
         VMAX = MAX (VMAX, PTS(I))
   10 CONTINUE

C ... Check for NaN (Not a Number).
C     NaN is defined to be not equal to any other number including itself
      IF (VMIN .NE. VMIN .OR. VMAX .NE. VMAX) THEN
        CALL PRTERR('WARNING',
     *    'Data in subroutine MINMAX contains "NaN" values. Check Data')
        if (vmin .ne. vmin) vmin = -1.0e30
        if (vmax .ne. vmax) vmax =  1.0e30
      END IF

      RETURN
      END
