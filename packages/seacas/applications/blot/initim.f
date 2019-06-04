C Copyright(C) 2009-2017 National Technology & Engineering Solutions of
C Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.
C     * Neither the name of NTESS nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

C $Log: initim.f,v $
C Revision 1.2  2009/03/25 12:36:45  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:03:44  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:52:35  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE INITIM (MAXTIM, WHONLY, NSTEPS, TIMES, WHOTIM,
     &   TMIN, TMAX, DELT, NINTV, NPTIMS, IPTIMS)
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
C   --   WHONLY - IN - true iff only whole times may be selected
C   --   NSTEPS - IN - the number of time steps
C   --   TIMES - IN - the database times
C   --   WHOTIM - IN - true iff TIMES(i) is whole (versus history) time step
C   --   TMIN - OUT - the minimum selected time; set to minimum time
C   --   TMAX - OUT - the maximum selected time; set to maximum time
C   --   DELT - OUT - the interval between selected times
C   --      (<0 = selected times); set to 0 (for all times) or interval
C   --   NINTV - IN/OUT - the number of times between tmin and tmax to select
C   --      (negative for zero interval);
C   --   NPTIMS - OUT - the number of selected times
C   --   IPTIMS - OUT - the selected time step numbers

C   --Routines Called:
C   --   MINMAX - (etclib) Find minimum and maximum value
C   --   MINMXL - (etclib) Find minimum and maximum value of selected values
C   --   NUMEQL - (etclib) Count the number of equal values

      LOGICAL WHONLY
      REAL TIMES(*)
      LOGICAL WHOTIM(*)
      INTEGER IPTIMS(*)

      IF (.NOT. WHONLY) THEN
         CALL MINMAX (NSTEPS, TIMES, TMIN, TMAX)
      ELSE
         CALL MINMXL (NSTEPS, WHOTIM, TIMES, TMIN, TMAX)
      END IF
      IF (.NOT. WHONLY) THEN
         NSTEPX = NSTEPS
      ELSE
         NSTEPX = NUMEQL (.TRUE., NSTEPS, WHOTIM)
      END IF

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
      CALL CALTIM (WHONLY, TMIN, TMAX, DELT, NINTV,
     &   NSTEPS, TIMES, WHOTIM, NPTIMS, IPTIMS)
      RETURN

      END
