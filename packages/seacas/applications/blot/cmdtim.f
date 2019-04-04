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

C $Log: cmdtim.f,v $
C Revision 1.2  2009/03/25 12:36:43  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 19:56:25  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:48:33  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE CMDTIM (INLINE,
     &   VERB, IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &   NSTEPS, TIMES, WHOTIM, WHONLY, TMIN, TMAX, DELT, NINTV,
     &   NPTIMS, IPTIMS)
C=======================================================================

C   --*** CMDTIM *** (TIMSEL) Process time step parameter command
C   --   Written by Amy Gilkey - revised 05/18/88
C   --
C   --CMDTIM processes a time step parameter command.  The commands are:
C   --   TMIN - sets the minimum selected time TMIN
C   --   TMAX - sets the maximum selected time TMAX
C   --   DELTIME - sets the selected time interval DELT
C   --   ZINTV - sets the number of selected times NINTV
C   --      with interval starting at TMIN
C   --   NINTV - sets the number of selected times NINTV
C   --      with interval starting at TMIN + offset
C   --   ALLTIMES - sets DELT = 0 for all times
C   --   TIMES - selects time steps by time
C   --   STEPS - selects time steps by number
C   --
C   --Parameters:
C   --   INLINE - IN/OUT - the parsed input lines for the log file
C   --   VERB - IN/OUT - the time step parameter command verb; set for SHOW
C   --   IFLD - IN/OUT - the field number
C   --   INTYP - IN - the input types from the free field reader
C   --   CFIELD - IN - the character fields
C   --   IFIELD - IN - the integer fields
C   --   RFIELD - IN - the real fields
C   --   NSTEPS - IN - the number of time steps
C   --   TIMES - IN - the database times
C   --   WHOTIM - IN - true iff TIMES(i) is whole (versus history) time step
C   --   WHONLY - IN - true iff only whole times may be selected
C   --   TMIN - IN/OUT - the minimum selected time
C   --   TMAX - IN/OUT - the maximum selected time
C   --   DELT - IN/OUT - the interval between selected times
C   --      (<0 = selected times)
C   --   NINTV - IN/OUT - the number of times between tmin and tmax to select
C   --      (negative for zero interval)
C   --   NPTIMS - IN/OUT - the number of selected times
C   --   IPTIMS - IN/OUT - the selected time step numbers

C   --Routines Called:
C   --   LOCREA - (ETCLIB) Find closest value
C   --   MINMAX - (ETCLIB) Find minimum and maximum value
C   --   MINMXL - (ETCLIB) Find minimum and maximum value of selected values
C   --   NUMEQL - (ETCLIB) Count the number of equal values

      CHARACTER*(*) INLINE(*)
      CHARACTER*(*) VERB
      INTEGER INTYP(*)
      CHARACTER*(*) CFIELD(*)
      INTEGER IFIELD(*)
      REAL RFIELD(*)
      REAL TIMES(*)
      LOGICAL WHOTIM(*)
      LOGICAL WHONLY
      INTEGER IPTIMS(*)

      LOGICAL FFEXST, FFMATC
      LOGICAL ISTIME
      CHARACTER*5 STR5

      IF (.NOT. WHONLY) THEN
         NSTEPX = NSTEPS
      ELSE
         NSTEPX = NUMEQL (.TRUE., NSTEPS, WHOTIM)
      END IF

      IF ((VERB .EQ. 'TMIN') .OR. (VERB .EQ. 'TMAX')) THEN
         CALL FFADDC (VERB, INLINE(1))
         IF (.NOT. WHONLY) THEN
            CALL MINMAX (NSTEPS, TIMES, TIMMIN, TIMMAX)
         ELSE
            CALL MINMXL (NSTEPS, WHOTIM, TIMES, TIMMIN, TIMMAX)
         END IF
         IF (VERB .EQ. 'TMIN') THEN
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &         'minimum time', TIMMIN, X, *120)
            TMIN = X
            CALL FFADDR (TMIN, INLINE(1))
         ELSE
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &         'maximum time', TIMMAX, X, *120)
            TMAX = X
            CALL FFADDR (TMAX, INLINE(1))
         END IF
         IF (TMIN .GT. TMAX) CALL PRTERR ('CMDWARN',
     &      'Minimum time is greater than maximum time')
         IF (DELT .LT. 0) THEN
            NINTV = 0
            DELT = 0.0
         END IF

      ELSE IF ((VERB .EQ. 'NINTV') .OR. (VERB .EQ. 'ZINTV')) THEN
         CALL FFADDC (VERB, INLINE(1))
         IF (VERB .EQ. 'NINTV') THEN
            NT = MIN (10, NSTEPX-1)
         ELSE
            NT = MIN (10, NSTEPX)
         END IF
         CALL FFINTG (IFLD, INTYP, IFIELD,
     &      'number of times', NT, N, *120)
         NINTV = N
         CALL FFADDI (NINTV, INLINE(1))
         IF (VERB .EQ. 'NINTV') THEN
            NINTV = N
         ELSE
            NINTV = -N
         END IF
         DELT = 999.0

      ELSE IF (VERB .EQ. 'DELTIME') THEN
         CALL FFADDC (VERB, INLINE(1))
         NT = - MIN (10, NSTEPX)
         IF (NT .LT. -1) THEN
            DDELT = (TMAX - TMIN) / (-NT-1)
            IF (DDELT .LE. 0.) DDELT = 1.0
         ELSE
            DDELT = 0.0
         END IF
         CALL FFREAL (IFLD, INTYP, RFIELD,
     &      'delta time', DDELT, X, *120)
         DELT = X
         CALL FFADDR (DELT, INLINE(1))
         NINTV = 0

      ELSE IF (VERB .EQ. 'ALLTIMES') THEN
         CALL FFADDC (VERB, INLINE(1))
         DELT = 0.0
         NINTV = 0

      ELSE IF ((VERB .EQ. 'TIMES') .OR. (VERB .EQ. 'STEPS')) THEN
         CALL FFADDC (VERB, INLINE(1))
         ISTIME = (VERB .EQ. 'TIMES')
         VERB = 'DELTIME'
         IF (DELT .GE. 0.0) THEN
            NINTV = NPTIMS
            DELT = -1.0
         END IF

C      --Reset is assumed, unless ADD is first parameter
         IF (.NOT. FFMATC (IFLD, INTYP, CFIELD, 'ADD', 1)) NINTV = 0

  100    CONTINUE
         IF (FFEXST (IFLD, INTYP)) THEN
            IF (ISTIME) THEN
               CALL FFREAL (IFLD, INTYP, RFIELD,
     &            'time, ignored', 0.0, T, *110)
               CALL FFADDR (T, INLINE(1))
               IF (.NOT. WHONLY) THEN
                  ISTEP = LOCREA (T, NSTEPS, TIMES)
               ELSE
                  ISTEP = LOCRL (T, NSTEPS, WHOTIM, TIMES)
               END IF

            ELSE
               CALL FFINTG (IFLD, INTYP, IFIELD,
     &            'step, ignored', 0, ISTEP, *110)
               IF ((ISTEP .LE. 0) .OR. (ISTEP .GT. NSTEPS)) THEN
                  CALL INTSTR (1, 0, ISTEP, STR5, LSTR)
                  CALL PRTERR ('CMDERR', 'Step ' // STR5(:LSTR)
     &               // ' does not exist, ignored')
                  GOTO 110
               ELSE IF (WHONLY) THEN
                  IF (.NOT. WHOTIM(ISTEP)) THEN
                     CALL INTSTR (1, 0, ISTEP, STR5, LSTR)
                     CALL PRTERR ('CMDERR', 'Step ' // STR5(:LSTR)
     &                  // ' is a history time step, ignored')
                     GOTO 110
                  END IF
               END IF
               CALL FFADDI (ISTEP, INLINE(1))
            END IF

            IF (LOCINT (ISTEP, NINTV, IPTIMS) .LE. 0) THEN
               NINTV = NINTV + 1
               IPTIMS(NINTV) = ISTEP
            END IF
  110       CONTINUE
            GOTO 100
         END IF
      END IF

      NPTIMS = -999
      CALL CALTIM (WHONLY, TMIN, TMAX, DELT, NINTV,
     &   NSTEPS, TIMES, WHOTIM, NPTIMS, IPTIMS)

  120 CONTINUE
      RETURN
      END
