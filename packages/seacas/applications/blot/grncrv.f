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

C $Log: grncrv.f,v $
C Revision 1.2  2009/03/25 12:36:44  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:02:36  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:51:47  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE GRNCRV (LABSID, ICURVE, NPTS, XPTS, YPTS, DRAWLN)
C=======================================================================

C   --*** GRNCRV *** (GRPLIB) Label Curve (PLT)
C   --   Written by Amy Gilkey - revised 10/07/87
C   --
C   --GRNCRV finds the device coordinates of the first, last or middle
C   --point of a curve within the plot window and labels the curve with
C   --the line number.
C   --
C   --Parameters:
C   --   LABSID - IN - label curve on "FIRST", "MIDDLE" or "LAST" point
C   --   ICURVE - IN - the curve number
C   --   NPTS - IN - the number of points in the curve
C   --   XPTS, YPTS - IN - the X and Y points
C   --   DRAWLN - IN - true if vector, false if points

C   --Routines Called:
C   --   MP2VC - (PLTLIB) Get device coordinates of vector
C   --   MP2PT - (PLTLIB) Get device coordinates of point
C   --   PLTGTT - (PLTLIB) Get text parameter
C   --      2 = (KSCHSZ) software character size
C   --   INTSTR - (STRLIB) Convert integer to string

      PARAMETER (KSCHSZ=2)

      CHARACTER*8 LABSID
      INTEGER ICURVE
      INTEGER NPTS
      REAL XPTS(*), YPTS(*)
      LOGICAL DRAWLN

      LOGICAL LDUM, PLTGTT
      CHARACTER*2 STR2

C   --Find the nearest line or point within the window

      IF (LABSID .EQ. 'FIRST') THEN
         ISTART = 1
         IEND = NPTS
         INC = 1
      ELSE IF (LABSID .EQ. 'MIDDLE') THEN
         ISTART = NPTS / 2
         IEND = NPTS
         INC = 1
      ELSE
         ISTART = NPTS
         IEND = 1
         INC = -1
      END IF

      IF (DRAWLN) THEN
         DO 100 I = ISTART, IEND-INC, INC
            MASK = -1
            CALL MP2VC (1, XPTS(I), YPTS(I), XPTS(I+INC), YPTS(I+INC),
     &         DX, DY, DXX, DYY, MASK)
            IF (MOD (MASK,2) .NE. 0) GOTO 120
  100    CONTINUE

      ELSE
         DO 110 I = ISTART, IEND, INC
            MASK = -1
            CALL MP2PT (1, XPTS(I), YPTS(I), DX, DY, MASK)
            IF (MOD (MASK,2) .NE. 0) GOTO 120
  110    CONTINUE
      END IF

      GOTO 130

C   --Write the curve number (location is dependent on the placement)

  120 CONTINUE
C   --Get curve number
      CALL INTSTR (1, 0, ICURVE, STR2, LNUM)

      LDUM = PLTGTT (KSCHSZ, VCS)
      IF (DRAWLN) THEN
         DXO = .25
      ELSE
         DXO = .50
      END IF
      IF (LABSID .EQ. 'FIRST') THEN
         DX = DX - (LNUM+DXO) * VCS
      ELSE
         DX = DX + DXO*VCS
      END IF
      IF (LABSID .EQ. 'MIDDLE') THEN
         DY = DY + .002
      ELSE
         DY = DY - .5*VCS
      END IF

      CALL GRTEXT (DX, DY, STR2)

  130 CONTINUE
      RETURN
      END
