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

C $Log: plshow.f,v $
C Revision 1.4  2009/03/25 12:36:46  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.3  1994/06/13 18:29:05  gdsjaar
C Fixed background and foreground color setting. (I think)
C
C Revision 1.2  1994/06/13  17:11:35  gdsjaar
C Fixed background, foreground, and softcharacters to check for full
C string rather than truncated at 8 characters.
C
c Revision 1.1  1994/04/07  20:07:16  gdsjaar
c Initial checkin of ACCESS/graphics/blotII2
c
c Revision 1.3  1993/09/24  17:32:42  gdsjaar
c Added an outline off/on command to toggle drawing of the view window outline
c
c Revision 1.2  1990/12/14  08:54:52  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE PLSHOW (SHOTYP, XYTYPE, MESHOK, TIMES, WHOTIM, IPTIMS)
C=======================================================================

C   --*** PLSHOW *** (BLOT) Display general parameter information
C   --   Written by Amy Gilkey - revised 04/26/88
C   --
C   --PLSHOW displays the general plot parameters.
C   --
C   --The SHOW options with the items they display are:
C   --   TMIN     - time step selection parameters
C   --   TMAX     -
C   --   NINTV    -
C   --   ZINTV    -
C   --   DELTIME  -
C   --   ALLTIMES -
C   --   STEPS    - the selected time step times
C   --   TIMES    -
C   --   QA       - plot QA information flag
C   --   AXIS     - plot axis numbering flag
C   --   LEGEND   - plot legend flag
C   --   CAPTION  - plot caption
C   --   SOFTCHAR - software vs hardware characters selected
C   --   FONT     - font to use
C   --   COLOR    - number of colors to use and maximum number of colors
C   --   SPECTRUM -
C   --   SNAP     - number of frames to snap for a movie
C   --   AUTO     - automatic vs user-directed plotting
C   --   BACKGROUND - set background color
C   --   FOREGROUND - set foreground color
C   --Parameters:
C   --   SHOTYP - IN - the expanded SHOW option string
C   --   XYTYPE - IN - true iff current program is an XY curve versus mesh plot
C   --   MESHOK - IN - true iff mesh can be displayed
C   --   TIMES - IN - the database times
C   --   WHOTIM - IN - true iff whole (versus history) time step
C   --   IPTIMS - IN - the selected time steps
C   --
C   --Common Variables:
C   --   Uses DOQA, DOLEG, DOAXIS, CAPTN of /LEGOPT/
C   --   Uses NPTIMS, TMIN, TMAX, DELT, NINTV of /TIMES/

      COMMON /LEGOPT/ DOQA(2), DOLEG(2), DOAXIS(2), DOBOX
      LOGICAL DOQA, DOLEG, DOAXIS, DOBOX
      COMMON /LEGOPC/ CAPTN(3,2)
      CHARACTER*80 CAPTN
      COMMON /TIMES/  NPTIMS, NPTIMW, TMIN, TMAX, DELT, NINTV,
     &   NTIMIN, WHONLY, HISTOK
      LOGICAL WHONLY, HISTOK

      include 'plcolr.blk'

      CHARACTER*(*) SHOTYP
      LOGICAL XYTYPE, MESHOK
      REAL TIMES(*)
      LOGICAL WHOTIM(*)
      INTEGER IPTIMS(*)

C *** Time step control ***

      IF ((SHOTYP .EQ. 'TMIN') .OR. (SHOTYP .EQ. 'TMAX')
     &   .OR. (SHOTYP .EQ. 'DELTIME')
     &   .OR. (SHOTYP .EQ. 'NINTV') .OR. (SHOTYP .EQ. 'ZINTV')
     &   .OR. (SHOTYP .EQ. 'ALLTIMES')
     &   .OR. (SHOTYP .EQ. 'STEPS') .OR. (SHOTYP .EQ. 'TIMES')) THEN
         CALL SHOTSP (WHONLY, TMIN, TMAX, DELT, NINTV, NPTIMS)

         IF (((SHOTYP .EQ. 'STEPS') .OR. (SHOTYP .EQ. 'TIMES'))
     &      .AND. (NPTIMS .GT. 0)) THEN
            CALL SHPTIM (WHONLY, NPTIMS, IPTIMS, TIMES, WHOTIM)
         END IF

      ELSE IF (SHOTYP .EQ. 'HISTORY') THEN
         IF (HISTOK) THEN
            IF (WHONLY) THEN
               WRITE (*, 10000) 'Select whole times steps only'
            ELSE
               WRITE (*, 10000) 'Select history and whole time steps'
            END IF
         ELSE
            WRITE (*, 10000) 'Select whole times steps only'
         END IF

C *** Display options ***

      ELSE IF ((SHOTYP .EQ. 'QA') .OR. (SHOTYP .EQ. 'LEGEND')
     &   .OR. (SHOTYP .EQ. 'AXIS') .OR. (SHOTYP .EQ. 'CAPTION') .OR.
     *     (SHOTYP .EQ. 'OUTLINE')) THEN
         IF (XYTYPE) CALL SHOLEG
     &      (SHOTYP, XYTYPE, DOQA(2), DOLEG(2), DOAXIS(2), CAPTN(1,2),
     *    DOBOX)
         IF (MESHOK) CALL SHOLEG
     &      (SHOTYP, .FALSE., DOQA(1), DOLEG(1), DOAXIS(1), CAPTN(1,1),
     *     DOBOX)

      ELSE IF ((SHOTYP .EQ. 'SOFTCHARACTERS') .OR. (SHOTYP .EQ. 'FONT')
     &   .OR. (SHOTYP .EQ. 'COLOR') .OR. (SHOTYP .EQ. 'SPECTRUM')
     &   .OR. (SHOTYP .EQ. 'SNAP') .OR. (SHOTYP .EQ. 'AUTO')) THEN
         CALL SHODEV (SHOTYP)

      ELSE IF (SHOTYP .EQ. 'BACKGROUND') THEN
         WRITE (*, 10000) 'Background color is ', BCKGND

      ELSE IF (SHOTYP .EQ. 'FOREGROUND') THEN
         WRITE (*, 10000) 'Foreground color is ', FORGND

      END IF

      RETURN
10000  FORMAT (1X, 5A)
      END
