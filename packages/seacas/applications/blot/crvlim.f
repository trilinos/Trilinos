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

C $Log: crvlim.f,v $
C Revision 1.3  2009/03/25 12:36:43  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.2  1997/11/11 21:44:42  gdsjaar
C Added check for NaN (not a number) in the curve limit determination.
C Prints a warning message and sets the min/max to be +/-1e30.
C
C Previous behavior was to hang...
C
C Revision 1.1  1994/04/07 19:57:25  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:49:08  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE CRVLIM (AXIS, TIMPLT, MAXPTS, NPTS, NSPVAR, NEPVAR,
     &   PLTVAL)
C=======================================================================

C   --*** CRVLIM *** (TPLOT) Calculate min/max value for plot data
C   --   Written by Amy Gilkey - revised 11/06/87
C   --
C   --CRVLIM calculates the minimum and maximum of the plot data.
C   --
C   --Parameters:
C   --   AXIS - IN - the axis to scale ('X' or 'Y')
C   --   TIMPLT - IN - true if time plot Y data versus X-Y data
C   --   MAXPTS - IN - the maximun number of points on a curve (PLTVAL length)
C   --   NPTS - IN - the number of points on each curve
C   --   NSPVAR, NEPVAR - IN - the starting and ending plot variable
C   --      indices for min/max calculation (PLTVAL index, /TPVARS/ index)
C   --   PLTVAL - IN - the data array
C   --
C   --Common Variables:
C   --   Sets XMIN, XMAX, YMIN, YMAX of /XYLIM/

      include 'xylim.blk'

      CHARACTER AXIS
      LOGICAL TIMPLT
      INTEGER NPTS(*)
      REAL PLTVAL(MAXPTS,*)

      IF ((.NOT .TIMPLT) .AND. (AXIS .EQ. 'X')) THEN
         XMIN =  1.0E+30
         XMAX = -1.0E+30

         N = NSPVAR
  100    CONTINUE
         IF (N .LE. NEPVAR) THEN
            DO 110 I = 1, NPTS(N)
               XMIN = MIN (XMIN, PLTVAL(I,N))
               XMAX = MAX (XMAX, PLTVAL(I,N))
  110       CONTINUE
            N = N + 2
            GOTO 100
         END IF
      END IF

      IF (AXIS .EQ. 'Y') THEN
         YMIN =  1.0E+30
         YMAX = -1.0E+30

         N = NSPVAR
  120    CONTINUE
         IF (N .LE. NEPVAR) THEN
            IF (.NOT. TIMPLT) N = N + 1
            DO 130 I = 1, NPTS(N)
               YMIN = MIN (YMIN, PLTVAL(I,N))
               YMAX = MAX (YMAX, PLTVAL(I,N))
  130       CONTINUE
            N = N + 1
            GOTO 120
         END IF
      END IF

C ... Check for NaN (Not a Number). 
C     NaN is defined to be not equal to any other number including itself
      if (ymin .ne. ymin .or. ymax .ne. ymax .or.
     *    xmin .ne. xmin .or. xmax .ne. xmax) then
        call prterr('WARNING',
     *    'The variable contains "NaN" (Not a Number) values.'//
     *    ' Check data.')
        if (xmin .ne. xmin) xmin = -1.0e30
        if (xmax .ne. xmax) xmax =  1.0e30
        if (ymin .ne. ymin) ymin = -1.0e30
        if (ymax .ne. ymax) ymax =  1.0e30
      end if

      RETURN
      END
