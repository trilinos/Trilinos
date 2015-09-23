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

C $Log: xyaxis.f,v $
C Revision 1.3  2009/03/25 12:36:49  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.2  1999/03/09 19:41:44  gdsjaar
C Fixed missing parameter definition of MSHBOR in blotII2.f
C
C Cleaned up parameters and common blocks in other routines.
C
C Revision 1.1  1994/04/07 20:17:50  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.3  1993/09/24  17:32:45  gdsjaar
c Added an outline off/on command to toggle drawing of the view window outline
c
c Revision 1.2  1990/12/14  08:59:35  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE XYAXIS (IXCRV, DOGRID, XLAB, YLAB, BLKCOL,  *)
C=======================================================================

C   --*** XYAXIS *** (XYPLOT) Draw the plot axis
C   --   Written by Amy Gilkey - revised 03/27/87
C   --
C   --XYAXIS draws the plot window with numbered and labeled axes.
C   --The device-to-user mapping is also done.
C   --
C   --Parameters:
C   --   IXCRV - IN - the plot type flag
C   --      0 if all curves on plot have the same scale
C   --      1 if first plot of a multi-scale plot
C   --      n if other plot of a multi-scale plot
C   --   DOGRID - IN - true iff a grid is to be drawn (IXCRV <= 1)
C   --   XLAB, YLAB - IN - labels for the X and Y axis
C   --   BLKCOL - IN/OUT - the user selected colors of the element blocks.
C   --                    BLKCOL(0) = 1 if the user defined material
C   --                                colors should be used in mesh plots.
C   --                              = -1 if program selected colors should
C   --                                be used.
C   --                    BLKCOL(i) = the user selected color of element
C   --                               block i:
C   --                                  -2 - no color selected by user.
C   --                                  -1 - black
C   --                                   0 - white
C   --                                   1 - red
C   --                                   2 - green
C   --                                   3 - yellow
C   --                                   4 - blue
C   --                                   5 - cyan
C   --                                   6 - magenta
C   --
C   --   * - the return statement if the cancel function is active
C   --
C   --Common Variables:
C   --   Uses DOGRID of /LEGOPT/
C   --   Uses ASPECT, IXSCAL, IYSCAL, XMIN, XMAX, YMIN, YMAX, XTICK, YTICK
C   --      of /XYLIM/
C   --   Uses CHLSIZ, DBORD0, DVIEW0 of /LAYOUT/

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4)

      PARAMETER (KSZAXE=28, KSZTIC=33, KSZGRI=34, KSZLIN=29)
      PARAMETER (KONGRI=15)
      PARAMETER (KXNUMS=22, KXLABS=23, KYNUMS=47, KYLABS=48)

      include 'dbnums.blk'
      include 'legopt.blk'
      include 'xylim.blk'
      include 'layout.blk'
      include 'plcol2.blk'

      LOGICAL DOGRID
      CHARACTER*(*) XLAB, YLAB

      LOGICAL GRABRT
      LOGICAL PLTGTG, PLTSTG, MPVIEW, MPORT2, LDUM
      REAL WVIEW(KTOP), DVIEW(KTOP), XVIEW(KTOP)
      CHARACTER*80 TXLAB, TYLAB
      INTEGER BLKCOL(0:NELBLK)

C   --Reset coordinate system

      CALL MPINIT

C   --Set up the viewport

      CALL GRVIEW (ASPECT, DVIEW0, DVIEW)
      WVIEW(KLFT) = XMIN
      WVIEW(KRGT) = XMAX
      WVIEW(KBOT) = YMIN
      WVIEW(KTOP) = YMAX

      LDUM = MPVIEW (DVIEW(KLFT), DVIEW(KRGT),
     &   DVIEW(KBOT), DVIEW(KTOP))
      LDUM = MPORT2 (WVIEW(KLFT), WVIEW(KRGT),
     &   WVIEW(KBOT), WVIEW(KTOP))

C   --Set up the axis labels

      IF (IXCRV .LE. 1) THEN
         TXLAB = XLAB
         TYLAB = YLAB
      ELSE
         TXLAB = ' '
         TYLAB = ' '
      END IF

      IF (GRABRT ()) RETURN 1
C   --NOTE* this is the last chance to cancel; PLT variables will be
C   --set and they need to be reset before exit

C   --Set up for grid drawing or not

      IF (IXCRV .LE. 1) THEN
         IF (DOGRID) THEN
            LDUM = PLTSTG (KONGRI, 1.0)
            LDUM = PLTGTG (KSZGRI, SZGRID)
            LDUM = PLTGTG (KSZTIC, SZTIC)
            LDUM = PLTSTG (KSZTIC, SZGRID)
         ELSE
            LDUM = PLTSTG (KONGRI, 0.0)
         END IF
      ELSE
         LDUM = PLTSTG (KONGRI, 0.0)
      END IF

      IF (.NOT. DOAXIS(2)) THEN
C      --Turn off the axis numbering (and the exponent on the label)
         LDUM = PLTSTG (KXNUMS, 0.0)
         LDUM = PLTSTG (KYNUMS, 0.0)
      END IF

C           Set colors for the axes and labels
      CALL UGRCOL (0, BLKCOL)
      LDUM = PLTSTG (10, COLFOR)
      LDUM = PLTSTG (16, COLFOR)
      LDUM = PLTSTG (17, COLFOR)
      LDUM = PLTSTG (43, COLFOR)
      LDUM = PLTSTG (44, COLFOR)
      LDUM = PLTSTG (45, COLFOR)
      LDUM = PLTSTG (46, COLFOR)

C   --For first plot of multi-scale plot, label axis 0 to 1

      IF (IXCRV .EQ. 1) THEN
         TXTICK = XTICK
         TYTICK = YTICK

         IF (IXSCAL .EQ. 'CURVE') THEN
            XVIEW(KLFT) = 0.0
            XVIEW(KRGT) = 1.0
            CALL EXPMAX (' ', XVIEW(KLFT), XVIEW(KRGT))
            IF (TXTICK .EQ. 0.0) TXTICK = 0.20
            XVIEW(KBOT) = WVIEW(KBOT)
            XVIEW(KTOP) = WVIEW(KTOP)
         ELSE IF (IYSCAL .EQ. 'CURVE') THEN
            XVIEW(KLFT) = WVIEW(KLFT)
            XVIEW(KRGT) = WVIEW(KRGT)
            XVIEW(KBOT) = 0.0
            XVIEW(KTOP) = 1.0
            CALL EXPMAX (' ', XVIEW(KBOT), XVIEW(KTOP))
            IF (TYTICK .EQ. 0.0) TYTICK = 0.20
         END IF

         CALL GRAXES (.FALSE., XVIEW, DVIEW,
     &      TXTICK, TYTICK, TXLAB, TYLAB)
         TXLAB = ' '
         TYLAB = ' '
         IF (DOGRID) LDUM = PLTSTG (KONGRI, 0.0)
      END IF

C   --Turn off everything for multi-scale plots (GRAXES must be called to
C   --set up scaling for PLTCUR routine)

      IF (IXCRV .GE. 1) THEN
         LDUM = PLTSTG (KXNUMS, 0.0)
         LDUM = PLTSTG (KXLABS, 0.0)
         LDUM = PLTSTG (KYNUMS, 0.0)
         LDUM = PLTSTG (KYLABS, 0.0)
         LDUM = PLTGTG (KSZTIC, SZTIC)
         LDUM = PLTSTG (KSZTIC, 0.0)
         LDUM = PLTGTG (KSZAXE, SZAXES)
         LDUM = PLTSTG (KSZAXE, 0.0)
      END IF

C   --Draw the axes

      CALL GRAXES (.FALSE., WVIEW, DVIEW, XTICK, YTICK, TXLAB, TYLAB)

C   --Reset PLT variables

      IF (IXCRV .LT. 1) THEN
         IF (DOGRID) LDUM = PLTSTG (KSZTIC, SZTIC)
      ELSE
         LDUM = PLTSTG (KSZTIC, SZTIC)
         LDUM = PLTSTG (KSZAXE, SZAXES)
      END IF

C   --Flush buffer, so label is complete at this point
      CALL PLTFLU

      RETURN
      END
