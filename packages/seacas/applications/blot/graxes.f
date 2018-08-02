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

C $Log: graxes.f,v $
C Revision 1.2  2009/03/25 12:36:44  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:02:03  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:51:23  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE GRAXES (XYSAME, WVIEW, DVIEW, WXATIC, WYATIC,
     &   TXLAB, TYLAB)
C=======================================================================

C   --*** GRAXES *** (GRPLIB) Draw axes (PLT)
C   --   Written by Amy Gilkey - revised 02/20/87
C   --
C   --GRAXES sets up and draws the axes for the plot.  It determines
C   --"good" numbers for the axes numbers (including the exponents) and
C   --sets the label/numbering size.
C   --
C   --Parameters:
C   --   XYSAME - IN - true iff the X and Y axis have the same scale;
C   --      i.e., they are the same type although the values may differ
C   --   WVIEW - IN - the window corners (left, right, bottom, top)
C   --      in window (user) coordinates
C   --   DVIEW - IN - the window corners (left, right, bottom, top)
C   --      in device coordinates
C   --   WXATIC, WYATIC - IN - the X and Y axis tick-mark interval;
C   --      default if equal zero or invalid
C   --   TXLAB, TYLAB - IN - the X and Y axis labels

C   --Routines Called:
C   --   PLTGPH - (PLTLIB) Draw the axes with labels
C   --   PLTSTG - (PLTLIB) Set graph parameter
C   --      1, 2 = (KXORIG, KYORIG) X, Y axis origin location
C   --      3, 4 = (KXLENG, KYLENG) X, Y axis length
C   --      11 = (KSCALE) axes parameters (see documentation)
C   --      22, 47 = (KXNUMS, KYNUMS) X, Y axis number size
C   --      23, 48 = (KXLABS, KYLABS) X, Y axis label size
C   --   GRAPAR - (GRPLIB) Select axis parameters

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4)
      PARAMETER (KXORIG=1, KYORIG=2, KXLENG=3, KYLENG=4, KSCALE=11)
      PARAMETER (KXNUMS=22, KXLABS=23, KYNUMS=47, KYLABS=48)

      LOGICAL XYSAME
      REAL WVIEW(KTOP), DVIEW(KTOP)
      REAL WXATIC, WYATIC
      CHARACTER*80 TXLAB, TYLAB

      LOGICAL LDUM, PLTSTG
      REAL BUF(11)

C   --Set device axis start and length

      DXAST = DVIEW(KLFT)
      DYAST = DVIEW(KBOT)
      DXALEN = DVIEW(KRGT) - DVIEW(KLFT)
      DYALEN = DVIEW(KTOP) - DVIEW(KBOT)

      LDUM = PLTSTG (KXORIG, DXAST)
      LDUM = PLTSTG (KYORIG, DYAST)
      LDUM = PLTSTG (KXLENG, DXALEN)
      LDUM = PLTSTG (KYLENG, DYALEN)

C   --Set axis minimum and maximum and tick intervals (and exponents and
C   --numbering size)

      TXATIC = WXATIC
      TYATIC = WYATIC
      CALL GRAPAR (XYSAME, WVIEW, DVIEW,
     &   WXALAB, WYALAB, WXAEND, WYAEND, TXATIC, TYATIC)

      BUF(1) = 4
      BUF(2) = WVIEW(KLFT)
      BUF(3) = WXALAB
      BUF(4) = WVIEW(KRGT)
      BUF(5) = TXATIC
      BUF(6) = 0.0
      BUF(7) = WVIEW(KBOT)
      BUF(8) = WYALAB
      BUF(9) = WVIEW(KTOP)
      BUF(10)= TYATIC
      BUF(11)= 0.0
      LDUM = PLTSTG (KSCALE, BUF)

C   --Draw the axes

      CALL PLTGPH (0., 0., 0, TXLAB, ' ', TYLAB, ' ')

      RETURN
      END
