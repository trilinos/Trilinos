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
      SUBROUTINE MINMXL (NPTS, PTOK, PTS, VMIN, VMAX)
C=======================================================================
C$Id: minmxl.f,v 1.2 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: minmxl.f,v $
CRevision 1.2  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:15:39  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:15:38  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:38  gdsjaar
c Initial revision
c 

C   --*** MINMXL *** (ETCLIB) Calculate min/max value in selected list
C   --   Written by Amy Gilkey - revised 11/10/87
C   --
C   --MINMXL calculates the minimum and maximum of the data selected by
C   --a logical value.
C   --
C   --Parameters:
C   --   NPTS - IN - the number of points
C   --   ISOK - IN - use point i iff ISOK(i)
C   --   PTS - IN - the points
C   --   VMIN, VMAX - OUT - the maximum and maximum value of the points

      INTEGER NPTS
      LOGICAL PTOK(*)
      REAL PTS(*)
      REAL VMIN, VMAX

      VMIN =  1.0E+30
      VMAX = -1.0E+30
      DO 10 I = 1, NPTS
         IF (PTOK(I)) THEN
            VMIN = MIN (VMIN, PTS(I))
            VMAX = MAX (VMAX, PTS(I))
         END IF
   10 CONTINUE

      RETURN
      END
