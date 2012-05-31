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

C $Log: expmax.f,v $
C Revision 1.2  2009/03/25 12:36:44  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:00:35  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:49:56  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE EXPMAX (LABSID, VMIN, VMAX)
C=======================================================================

C   --*** EXPMAX *** (BLOT) Expands min/max values
C   --   Written by Amy Gilkey - revised 10/07/87
C   --
C   --EXPMAX expands the minimum and maximum values by 2.5% each.
C   --It also expands the appropriate limit by 2.5% for numbering.
C   --
C   --Parameters:
C   --   LABSID - IN - if 'FIRST', expand for numbering on left side;
C   --      if 'LAST', expand for numbering on right side
C   --   VMIN/VMAX - IN/OUT - the minimum and maximum axis values

      CHARACTER*(*) LABSID
      REAL VMIN, VMAX

      IF (VMIN .EQ. VMAX) THEN
         IF (VMIN .EQ. 0.0) THEN
            VMIN = -1.0
            VMAX = 1.0
         ELSE
            VRNG = 0.025 * ABS(VMIN)
            VMIN = VMIN - VRNG
            VMAX = VMAX + VRNG
         END IF

      ELSE
         VRNG = 0.025 * (VMAX - VMIN)
         VMIN = VMIN - VRNG
         VMAX = VMAX + VRNG
         IF (LABSID .EQ. 'FIRST') VMIN = VMIN - 0.025 * (VMAX - VMIN)
         IF (LABSID .EQ. 'LAST') VMAX = VMAX + 0.025 * (VMAX - VMIN)
      END IF

      RETURN
      END
