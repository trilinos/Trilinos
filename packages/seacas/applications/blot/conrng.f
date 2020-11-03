C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

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
