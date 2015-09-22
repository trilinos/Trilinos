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

C=======================================================================
      SUBROUTINE MINMAX (NPTS, PTS, VMIN, VMAX)
C=======================================================================
C$Id: minmax.f,v 1.3 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: minmax.f,v $
CRevision 1.3  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.2  1997/11/11 21:45:59  gdsjaar
CAdded check for NaN (not a number) in the min/max limit determination.
CPrints a warning message and sets the min/max to be +/-1e30.
C
CRevision 1.1.1.1  1990/08/14 16:15:37  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:15:36  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:38  gdsjaar
c Initial revision
c 

C   --*** MINMAX *** (ETCLIB) Calculate min/max value
C   --   Written by Amy Gilkey - revised 11/10/87
C   --
C   --MINMAX calculates the minimum and maximum of the data.
C   --
C   --Parameters:
C   --   NPTS - IN - the number of points
C   --   PTS - IN - the points
C   --   VMIN, VMAX - OUT - the maximum and maximum value of the points

      INTEGER NPTS
      REAL PTS(*)
      REAL VMIN, VMAX

      VMIN =  1.0E+30
      VMAX = -1.0E+30
      DO 10 I = 1, NPTS
         VMIN = MIN (VMIN, PTS(I))
         VMAX = MAX (VMAX, PTS(I))
   10 CONTINUE

C ... Check for NaN (Not a Number). 
C     NaN is defined to be not equal to any other number including itself
      IF (VMIN .NE. VMIN .OR. VMAX .NE. VMAX) THEN
        CALL PRTERR('WARNING',
     *    'Data in subroutine MINMAX contains "NaN" values. Check Data')
        if (vmin .ne. vmin) vmin = -1.0e30
        if (vmax .ne. vmax) vmax =  1.0e30
      END IF

      RETURN
      END
