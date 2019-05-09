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

C $Log: xyshow.f,v $
C Revision 1.2  2009/03/25 12:36:49  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:18:15  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.3  1994/01/14  19:47:47  gdsjaar
c Added capability to write both grafaid and xmgr neutral files.
c
c Revision 1.2  1990/12/14  08:59:39  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE XYSHOW (SHOTYP)
C=======================================================================

C   --*** XYSHOW *** (XYPLOT) Display XY plot parameter information
C   --   Written by Amy Gilkey - revised 12/24/87
C   --
C   --XYSHOW displays the XY plot parameters.
C   --
C   --The SHOW options with the items they display are:
C   --   RATIOXY  - the X to Y axis length ratio
C   --   XSCALE   - the X axis minimum and maximum if user-defined
C   --   YSCALE   - the Y axis minimum and maximum if user-defined
C   --   XTICK    - the X axis tick interval
C   --   YTICK    - the Y axis tick interval
C   --   XLABEL   - the X axis label
C   --   YLABEL   - the Y axis label
C   --   ACURVE   - the neutral file curve name, number and increment
C   --   NCURVE   -
C   --   GRID     - the grid/no grid option
C   --   LINES    - the line, symbol options
C   --   SYMBOLS  -
C   --   CRVNUM   - the curve numbering option
C   --   OVERLAY  - the overlay / separate curves or times option
C   --   SAMESCAL - the curve scaling option
C   --   NORMAL   -
C   --
C   --Parameters:
C   --   SHOTYP - IN - the expanded SHOW option string
C   --
C   --Common Variables:
C   --   Uses DOGRID, LINTYP, ISYTYP, LABSID, OVERLY, OVERTM, IAXSCA of /XYOPT/
C   --   Uses ASPECT, IXSCAL, IYSCAL, XMIN, XMAX, YMIN, YMAX, XTICK, YTICK
C   --      of /XYLIM/
C   --   Uses XLAB, YLAB of /XYLAB/
C   --   Uses NEUOPN, NUMCRV, INCCRV, CRVNAM of /NEUTR./

      include 'xyopt.blk'
      include 'xylim.blk'
      include 'xylab.blk'
      include 'neutr.blk'

      CHARACTER*(*) SHOTYP

      CHARACTER*80 STRING
      CHARACTER*20 RSTR(2)
      REAL RNUM(2)
      CHARACTER*5 STRI1, STRI2

      IF (SHOTYP .EQ. 'RATIOXY') THEN
         IF (ASPECT .GE. .99) THEN
            CALL NUMSTR (1, 3, ASPECT, RSTR(1), LSTR)
         ELSE
            CALL NUMSTR (1, 2, ASPECT, RSTR(1), LSTR)
         END IF
         WRITE (*, 10000) 'X to Y aspect ratio = ', RSTR(1)(:LSTR)

      ELSE IF (SHOTYP .EQ. 'XSCALE') THEN
         IF (IXSCAL .EQ. 'SET') THEN
            RNUM(1) = XMIN
            RNUM(2) = XMAX
            CALL NUMSTR (2, 4, RNUM, RSTR, LSTR)
            WRITE (STRING, '(5A)') 'X axis scaling: ',
     &         RSTR(1)(:LSTR), ' to ', RSTR(2)(:LSTR)
            CALL SQZSTR (STRING, LSTR)
            WRITE (*, 10000) STRING(:LSTR)
         ELSE
            WRITE (*, 10000) 'X axis will be automatically scaled'
         END IF

      ELSE IF (SHOTYP .EQ. 'YSCALE') THEN
         IF (IYSCAL .EQ. 'SET') THEN
            RNUM(1) = YMIN
            RNUM(2) = YMAX
            CALL NUMSTR (2, 4, RNUM, RSTR, LSTR)
            WRITE (STRING, '(5A)') 'Y axis scaling: ',
     &         RSTR(1)(:LSTR), ' to ', RSTR(2)(:LSTR)
            CALL SQZSTR (STRING, LSTR)
            WRITE (*, 10000) STRING(:LSTR)
         ELSE
            WRITE (*, 10000) 'Y axis will be automatically scaled'
         END IF

      ELSE IF (SHOTYP .EQ. 'XTICK') THEN
         IF (XTICK .EQ. 0.0) THEN
            WRITE (*, 10000) 'X axis tick interval automatically scaled'
         ELSE
            CALL NUMSTR (1, 4, XTICK, RSTR(1), LSTR)
            WRITE (*, 10000) 'X axis tick interval = ', RSTR(1)(:LSTR)
         END IF

      ELSE IF (SHOTYP .EQ. 'YTICK') THEN
         IF (YTICK .EQ. 0.0) THEN
            WRITE (*, 10000) 'Y axis tick interval automatically scaled'
         ELSE
            CALL NUMSTR (1, 4, YTICK, RSTR(1), LSTR)
            WRITE (*, 10000) 'Y axis tick interval = ', RSTR(1)(:LSTR)
         END IF

      ELSE IF (SHOTYP .EQ. 'XLABEL') THEN
         IF (XLAB .EQ. ' ') THEN
            WRITE (*, 10000) 'X axis label is not defined'
         ELSE
            WRITE (*, 10000) 'X axis label: ', XLAB(:LENSTR(XLAB))
         END IF

      ELSE IF (SHOTYP .EQ. 'YLABEL') THEN
         IF (YLAB .EQ. ' ') THEN
            WRITE (*, 10000) 'Y axis label is not defined'
         ELSE
            WRITE (*, 10000) 'Y axis label: ', YLAB(:LENSTR(YLAB))
         END IF

      ELSE IF ((SHOTYP .EQ. 'ACURVE')
     &   .OR. (SHOTYP .EQ. 'NCURVE')) THEN
         WRITE (*, 10000) 'Neutral file curve name = ', CRVNAM
         CALL INTSTR (1, 0, NUMCRV, STRI1, L1)
         CALL INTSTR (1, 0, INCCRV, STRI2, L2)
         WRITE (*, 10000) '   Curve number = ',
     &      STRI1(:L1), ', incremented by ', STRI2(:L2)

      ELSE IF (SHOTYP .EQ. 'GRID') THEN
         IF (DOGRID) THEN
            WRITE (*, 10000) 'Draw grid on plot'
         ELSE
            WRITE (*, 10000) 'No grid on plot'
         END IF

      ELSE IF ((SHOTYP .EQ. 'LINES') .OR. (SHOTYP .EQ. 'SYMBOLS')) THEN
         IF (ISYTYP .LE. -1) THEN
            STRING = 'symbols'
         ELSE IF (ISYTYP .EQ. 0) THEN
            STRING = ' '
         ELSE
            CALL INTSTR (1, 0, ISYTYP, STRI1, LSTR)
            STRING = 'symbol ' // STRI1(:LSTR)
         END IF
         LSTR = LENSTR(STRING)
         IF (LINTYP .LE. -1) THEN
            IF (ISYTYP .EQ. 0) THEN
               WRITE (*, 10000) 'Plot with varying line'
            ELSE
               WRITE (*, 10000) 'Plot with varying line and ',
     &            STRING(:LSTR)
            END IF
         ELSE IF (LINTYP .EQ. 0) THEN
            WRITE (*, 10000) 'Plot with ', STRING(:LSTR)
         ELSE
            IF (ISYTYP .EQ. 0) THEN
               WRITE (*, 10000) 'Plot with solid line'
            ELSE
               WRITE (*, 10000) 'Plot with solid line and ',
     &            STRING(:LSTR)
            END IF
         END IF

      ELSE IF (SHOTYP .EQ. 'CRVNUM') THEN
         IF (LABSID .EQ. 'NONE') THEN
            WRITE (*, 10000) 'Do not number curves'
         ELSE IF (LABSID .EQ. 'FIRST') THEN
            WRITE (*, 10000) 'Number curves at first point'
         ELSE IF (LABSID .EQ. 'MIDDLE') THEN
            WRITE (*, 10000) 'Number curves at middle point'
         ELSE IF (LABSID .EQ. 'LAST') THEN
            WRITE (*, 10000) 'Number curves at last point'
         END IF

      ELSE IF (SHOTYP .EQ. 'OVERLAY') THEN
         IF (OVERLY) THEN
            WRITE (*, 10000) 'Overlay curves for all variables'
         ELSE IF (OVERTM) THEN
            WRITE (*, 10000) 'Overlay curves for all times'
         ELSE
            WRITE (*, 10000) 'Separate curves'
         END IF

      ELSE IF ((SHOTYP .EQ. 'SAMESCAL')
     &   .OR. (SHOTYP .EQ. 'NORMAL')) THEN
         IF (IAXSCA .EQ. 'ALL') THEN
            WRITE (*, 10000) 'Axis scaling defaults to',
     &         ' same scale for all curves on all plots'
         ELSE IF (IAXSCA .EQ. 'PLOT') THEN
            WRITE (*, 10000) 'Axis scaling defaults to',
     &         ' same scale for all curves on each plot'
         ELSE IF (IAXSCA .EQ. 'CURVE') THEN
            WRITE (*, 10000) 'Axis scaling defaults to',
     &         ' normalized scale for each curve'
         END IF

      END IF

      RETURN
10000  FORMAT (1X, 10A)
      END
