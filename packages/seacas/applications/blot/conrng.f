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

C $Log: conrng.f,v $
C Revision 1.3  2009/03/25 12:36:43  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.2  2003/10/02 17:01:51  gdsjaar
C Fixed the setting of an artificial range for a constant negative
C value. It was incorrectly making max < min and then iterating until
C overflow. Changed to subtract and add the absolute value of the
C constant value.
C
C Removed call to initialize random file since there is none with
C exodusII
C
C Removed saving of memory pointers in rndvar; instead just find them
C when reading an element variable.
C
C Revision 1.1  1994/04/07 19:57:13  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:49:00  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE CONRNG (ISLINE, FMINI, FMAXI, NCNTR, DELC, CMIN, CMAX)
C=======================================================================

C   --*** CONRNG *** (DETOUR) Calculate contour range
C   --   Written by Amy Gilkey - revised 06/09/87
C   --
C   --CONRNG calculates a new contour range from the minimum and maximum
C   --values.  Note that the contour interval and range are rounded to
C   --get "nice" numbers.  For the line contour, the contour minimum and
C   --maximum should fall between the actual value range (about one-half
C   --an interval in).  For the paint contour, the contour range must
C   --cover the actual value range.
C   --
C   --Parameters:
C   --   ISLINE - IN - true if line contour, else paint contour
C   --   FMINI, FMAXI - IN - the minimum and maximum value in the contour plot
C   --   NCNTR - IN/OUT - the number of contours defined
C   --   DELC - OUT - the contour interval
C   --   CMIN - OUT - the minimum contour value
C   --   CMAX - OUT - the maximum contour value

      LOGICAL ISLINE

      CHARACTER*8 STRING

      IF (NCNTR .LE. 0) NCNTR = 6

C   --Check for minimum equal maximum

      IF (FMINI .EQ. FMAXI) THEN
         CALL PRTERR ('CMDWARN', 'Contour variable does not vary'
     &      // ' - an artificial range is supplied')
         IF (FMINI .EQ. 0.0) THEN
            FMIN = -1.0
            FMAX = 1.0
         ELSE
            FMIN = FMINI - .05 * ABS(FMINI)
            FMAX = FMINI + .05 * ABS(FMINI)
         END IF
      ELSE
         FMIN = FMINI
         FMAX = FMAXI
      END IF

      IF (ISLINE) THEN

C      --Compute decimal scale factor

         DELC = (FMAX - FMIN) / NCNTR
         WRITE (STRING, '(E8.1)') DELC
         STRING(4:4) = '1'
         READ (STRING, '(E8.1)') P10

C      --Compute rounded increment and offset

         DELC = .5 * P10 * NINT (2. * DELC / P10)
         CMIN = .5 *
     &      (P10 * NINT ((FMAX + FMIN) / P10) - (NCNTR-1) * DELC)
         CMAX = CMIN + (NCNTR-1) * DELC

      ELSE
         FRNG = FMAX - FMIN
  100    CONTINUE

C      --Compute decimal scale factor

         DELC = FRNG / NCNTR
         WRITE (STRING, '(E8.1)') DELC
         STRING(4:4) = '1'
         READ (STRING, '(E8.1)') P10

C      --Compute rounded increment and offset

         DELC = .25 * P10 * NINT (4. * DELC / P10)
         CMIN = .5 * (P10 * NINT ((FMAX + FMIN) / P10) - NCNTR * DELC)
         IF ((FMIN .GE. 0.0) .AND. (CMIN .LT. 0.0)) CMIN = 0.0
         CMAX = CMIN + NCNTR * DELC

C      --Adjust and recompute, if needed

         IF ((CMIN .GT. FMIN) .OR. (CMAX .LT. FMAX)) THEN
            IF (CMIN .GT. FMIN) FRNG = FRNG + .25*DELC
            IF (CMAX .LT. FMAX) FRNG = FRNG + .25*DELC
            GOTO 100
         END IF

      END IF

      RETURN
      END
