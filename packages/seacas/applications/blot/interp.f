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

C $Log: interp.f,v $
C Revision 1.2  2009/03/25 12:36:45  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:03:56  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:52:45  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      LOGICAL FUNCTION INTERP_BL (CNTR, F1, F2, PSI)
C=======================================================================

C   --*** INTERP_BL *** (DETOUR) Compute interception point
C   --   Written by Amy Gilkey - revised 11/21/85
C   --   D. P. Flanagan, 3/25/83
C   --
C   --INTERP_BL tests if a contour value falls within an interval.
C   --A degenerate interval fails regardless of the contour value. If the
C   --test is passed, the nomalized interval coordinate of the contour
C   --value is computed.  The returned function value is true only if the
C   --interval test is passed.
C   --
C   --Parameters:
C   --   CNTR - IN - the contour value
C   --   F1 - IN - the interval origin value
C   --   F2 - IN - the interval terminus value
C   --   PSI - OUT - normalized interval coordinate

C   --Test if coordinate lies within the interval

      INTERP_BL = (F2 .NE. F1) .AND.
     &   (((F1 .LE. CNTR) .AND. (CNTR .LE. F2))
     &   .OR. ((F2 .LE. CNTR) .AND. (CNTR .LE. F1)))

C   --Compute normalized interval coordinate

      IF (INTERP_BL) PSI = (CNTR - F1) / (F2 - F1)
      RETURN
      END
