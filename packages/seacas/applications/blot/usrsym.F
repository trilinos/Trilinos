C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE USRSYM (ISYTYP, IS3DIM, X0, Y0, Z0, EVAR, SCALE)
C=======================================================================

C   --*** USRSYM *** (DETOUR) Plot element symbol (slip/angle version)
C   --   Written by Amy Gilkey - revised 06/15/87
C   --   D. P. Flanagan, 12/08/83
C   --
C   --If the element variable name begins with 'CRACK', USRSYM displays
C   --an 'X' to indicate a sliding joint or an 'O' to indicate an opening
C   --joint (Slip plot).
C   --
C   --If the element variable name begins with 'ANGLE', USRSYM draws a
C   --line segment at a specified angle (Angle plot). A '*' is drawn to
C   --indicate rubblization.
C   --
C   --Parameters:
C   --   ISYTYP - IN - the symbol type (as in MODTYP of /DETOPT/)
C   --   IS3DIM - IN - true iff 3D versus 2D
C   --   X0, Y0, Z0 - IN - the element center coordinates
C   --   EVAR - IN - the element variable
C   --      Crack angle for angle plot:
C   --         >360 indicates no crack
C   --      Crack state indicator for slip plot:
C   --         -1 = joint sliding
C   --          0 = joint stationary
C   --         +1 = joint opening
C   --   SCALE - IN - the scale factor for crack angle,
C   --      default is device to window conversion

      PARAMETER (DTR = .0174532925)

      common /debugc/ cdebug
      common /debugn/ idebug
      character*8 cdebug

      CHARACTER*(*) ISYTYP
      LOGICAL IS3DIM

      IF (ISYTYP .EQ. 'ANGLE') THEN

C      --Angle plot

         IF (EVAR .GT. 360.0) THEN
            RETURN

         ELSE IF (EVAR .LT. -360.0) THEN

C         --Indicate multiple failures

#if NeedsDoubleEscape
            CALL MPD2SY (1, X0, Y0, '\\CX')
            CALL MPD2SY (1, X0, Y0, '\\CCS')
#else
            CALL MPD2SY (1, X0, Y0, '\CX')
            CALL MPD2SY (1, X0, Y0, '\CCS')
#endif
         ELSE

C         --Indicate crack angle

            SX0 = .01 * SCALE * COS(EVAR*DTR)
            SY0 = .01 * SCALE * SIN(EVAR*DTR)
            CALL MPD2VC (1, X0-SX0, Y0-SY0, X0+SX0, Y0+SY0)
         END IF

      ELSE IF (ISYTYP .EQ. 'CRACK') THEN

C      --Slip plot

         IF (EVAR .EQ. 0.0) THEN
            RETURN

         ELSE IF (EVAR .LT. 0.0) THEN

C         --Indicate sliding joint

#if NeedsDoubleEscape
            CALL MPD2SY (1, X0, Y0, '\\CX')
#else
            CALL MPD2SY (1, X0, Y0, '\CX')
#endif
         ELSE

C         --Indicate opening joint

#if NeedsDoubleEscape
            CALL MPD2SY (1, X0, Y0, '\\CCI')
#else
            CALL MPD2SY (1, X0, Y0, '\CCI')
#endif
         END IF

      else if (cdebug .eq. 'SYMBOL') then
#if NeedsDoubleEscape
         call mpd2sy (1, x0, y0, '\\cx')
#else
         call mpd2sy (1, x0, y0, '\cx')
#endif

      END IF

      RETURN
      END
