C Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software.
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
C     * Neither the name of Sandia Corporation nor the names of its
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

C $Log: usrsym.F,v $
C Revision 1.3  2009/03/25 12:36:49  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.2  1997/06/26 17:51:37  caforsy
C Substituted if NeedsDoubleEscape for ifdef sun.  This define is needed on cougar, sun and ibm platforms
C
C Revision 1.1  1994/04/07 20:17:22  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
CRevision 1.3  1993/07/08  21:57:17  gdsjaar
CReplaced SUN with sun in all ifdef statements.
C
CRevision 1.2  1990/12/14  08:59:27  gdsjaar
CAdded RCS Id and Log to all files
C
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
