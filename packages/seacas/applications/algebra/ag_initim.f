C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE INITIM (MAXTIM, NSTEPS, TIMES, TMIN, TMAX,
     &                   DELT, NINTV, NPTIMS, IPTIMS)
C=======================================================================

C   --*** INITIM *** (TIMSEL) Initialize time step parameters
C   --   Written by Amy Gilkey - revised 11/12/87
C   --
C   --INITIM initializes the time step parameters.  The initialization
C   --is dependent on the MAXTIM parameter.
C   --
C   --Parameters:
C   --   MAXTIM - IN - the initialization parameter
C   --      =0 - select all times
C   --      >0 - select NINTV interval times (up to number of steps - 1)
C   --           with interval starting at TMIN + offset
C   --      <0 - select -NINTV interval times (up to number of steps)
C   --           with interval starting at TMIN
C   --   NSTEPS - IN -  the number of time steps
C   --   TIMES  - IN -  the database times
C   --   TMIN   - OUT - the minimum selected time; set to minimum time
C   --   TMAX   - OUT - the maximum selected time; set to maximum time
C   --   DELT   - OUT - the interval between selected times
C   --                  (<0 = selected times)
C   --                  set to 0 (for all times) or interval
C   --   NINTV  - IN/OUT - the number of times between tmin and tmax to select
C   --                     (negative for zero interval);
C   --   NPTIMS - OUT - the number of selected times
C   --   IPTIMS - OUT - the selected time step numbers

C   --Routines Called:
C   --   MINMAX - (ETCLIB) Find minimum and maximum value
C   --   MINMXL - (ETCLIB) Find minimum and maximum value of selected values
C   --   NUMEQL - (ETCLIB) Count the number of equal values

      INTEGER MAXTIM
      INTEGER NSTEPS
      REAL TIMES(*)
      REAL TMIN
      REAL TMAX
      INTEGER IPTIMS(*)

       CALL MINMAX (NSTEPS, TIMES, TMIN, TMAX)

       NSTEPX = NSTEPS

      IF (MAXTIM .NE. 0) THEN
         IF (MAXTIM .GT. 0) THEN
            NINTV = MIN (MAXTIM, NSTEPX-1)
         ELSE
            NINTV = - MIN (MAXTIM, NSTEPX)
         END IF
         DELT = 999.0

      ELSE IF (MAXTIM .EQ. 0) THEN
         DELT = 0.0
         NINTV = 0
      END IF

      NPTIMS = -999
      CALL CALTIM (TMIN, TMAX, DELT, NINTV,
     &             NSTEPS, TIMES, NPTIMS, IPTIMS)

      RETURN
      END
