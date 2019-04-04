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

C $Log: caltim.f,v $
C Revision 1.2  2009/03/25 12:36:42  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 19:55:14  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:47:51  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE CALTIM (WHONLY, TMIN, TMAX, DELT, NINTV,
     &   NSTEPS, TIMES, WHOTIM, NPTIMS, IPTIMS)
C=======================================================================

C   --*** CALTIM *** (TIMSEL) Calculate the selected times
C   --   Written by Amy Gilkey - revised 11/10/87
C   --
C   --CALTIM calculates the selected times based on the time
C   --parameters.  The calculation is done only if the number of selected
C   --times is less than 0.
C   --
C   --Parameters:
C   --   WHONLY - IN - true iff only whole times may be selected
C   --   TMIN - IN - the minimum selected time
C   --   TMAX - IN - the maximum selected time
C   --   DELT - IN - the interval between selected times
C   --      0 < - user selected times mode
C   --      0 = - all available times mode
C   --      0 > - uniform time interval mode (dependent on NINTV if <> 0)
C   --   NINTV - IN - in user selected times mode, the number of times;
C   --      in uniform time interval mode, negative if zero interval
C   --   NSTEPS - IN - the number of time steps
C   --   TIMES - IN - the database times
C   --   WHOTIM - IN - true iff TIMES(i) is whole (versus history) time step
C   --   NPTIMS - IN/OUT - the number of selected times; the
C   --      calculation is done only if this is negative
C   --   IPTIMS - IN/OUT - the selected time steps

C   --Routines Called:
C   --   LOCREA - (ETCLIB) Locate closest value
C   --   NUMEQL - (ETCLIB) Count the number of equal values

      LOGICAL WHONLY
      REAL TMIN, TMAX, DELT
      INTEGER NINTV
      INTEGER NSTEPS
      REAL TIMES(*)
      LOGICAL WHOTIM(*)
      INTEGER NPTIMS
      INTEGER IPTIMS(*)

      LOGICAL DOSORT

      IF (NPTIMS .GE. 0) RETURN

      IF (.NOT. WHONLY) THEN
         NSTEPX = NSTEPS
      ELSE
         NSTEPX = NUMEQL (.TRUE., NSTEPS, WHOTIM)
      END IF

      DOSORT = .FALSE.

      IF (NSTEPX .EQ. 0) THEN

C      --Check for no time steps

         NPTIMS = 0
         TMIN = 0.0
         TMAX = 0.0
         DELT = 0.0
         NINTV = 0

      ELSE IF (DELT .GT. 0.0) THEN

C      --Select database time steps between TMIN and TMAX at interval DELT

         IF (TMIN .GT. TMAX) THEN
            NPTIMS = 0

         ELSE IF (TMIN .EQ. TMAX) THEN
            IF (.NOT. WHONLY) THEN
               IPTIMS(1) = LOCREA (TMIN, NSTEPS, TIMES)
            ELSE
               IPTIMS(1) = LOCRL (TMIN, NSTEPS, WHOTIM, TIMES)
            END IF
            NPTIMS = 1

         ELSE
            IF (NINTV .GT. 0) THEN
               IS = 1
               IE = NINTV
               DELT = (TMAX - TMIN) / NINTV
            ELSE IF (NINTV .LT. 0) THEN
               IS = 0
               IE = - NINTV - 1
               IF (NINTV .NE. -1) DELT = (TMAX - TMIN) / (-NINTV-1)
            ELSE
               IS = 0
               TIMMAX = TMAX + 0.01*DELT
               IE = INT ((TIMMAX - TMIN) / DELT)
            END IF
            ILAST = 0
            N = 0

            DO 100 I = IS, IE
               T = TMIN + DELT*I
               IF (.NOT. WHONLY) THEN
                  ITHIS = LOCREA (T, NSTEPS, TIMES)
               ELSE
                  ITHIS = LOCRL (T, NSTEPS, WHOTIM, TIMES)
               END IF
               IF (ILAST .NE. ITHIS) THEN
                  IWHERE = LOCINT (ITHIS, N, IPTIMS)
                  IF (IWHERE .LE. 0) THEN
                     N = N + 1
                     IPTIMS(N) = ITHIS
                     IF (ILAST .GE. ITHIS) DOSORT = .TRUE.
                     ILAST = ITHIS
                  END IF
               END IF
  100       CONTINUE
            NPTIMS = N
         END IF

      ELSE IF (DELT .EQ. 0.0) THEN

C      --Select all the database time steps between TMIN and TMAX

         N = 0
         IF (.NOT. WHONLY) THEN
            ISTIM = LOCREA (TMIN, NSTEPS, TIMES)
            IETIM = LOCREA (TMAX, NSTEPS, TIMES)
         ELSE
            ISTIM = LOCRL (TMIN, NSTEPS, WHOTIM, TIMES)
            IETIM = LOCRL (TMAX, NSTEPS, WHOTIM, TIMES)
         END IF
         TIMMIN = TIMES(ISTIM)
         TIMMAX = TIMES(IETIM)
         DO 110 ITHIS = 1, NSTEPS
            IF ((TIMES(ITHIS) .GE. TIMMIN)
     &         .AND. (TIMES(ITHIS) .LE. TIMMAX)) THEN
               N = N + 1
               IPTIMS(N) = ITHIS
            END IF
  110    CONTINUE
         NPTIMS = N

      ELSE

C      --Select the database time steps nearest a specified list of times
C      --(in IPTIMS)

         ILAST = 0
         DO 120 I = 1, NINTV
            ITHIS = IPTIMS(I)
            IF (ILAST .GE. ITHIS) DOSORT = .TRUE.
            ILAST = ITHIS
  120    CONTINUE
         NPTIMS = NINTV
      END IF

C   --Sort the selected time steps numbers (in IPTIMS)

      IF (DOSORT) THEN
         ILAST = 0
         N = 0
         DO 140 I = 1, NPTIMS
            ITHIS = IPTIMS(I)
            IXTHIS = I
            DO 130 J = I+1, NPTIMS
               IF (IPTIMS(J) .LT. ITHIS) THEN
                  ITHIS = IPTIMS(J)
                  IXTHIS = J
               END IF
  130       CONTINUE
            IPTIMS(IXTHIS) = IPTIMS(I)
            IF (ILAST .NE. ITHIS) THEN
               N = N + 1
               IPTIMS(N) = ITHIS
               ILAST = ITHIS
            END IF
  140    CONTINUE
         NPTIMS = N
      END IF

      RETURN
      END
