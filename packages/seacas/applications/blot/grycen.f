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

C $Log: grycen.f,v $
C Revision 1.2  2009/03/25 12:36:45  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:03:00  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:52:07  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE GRYCEN (CHLSIZ, TOPLIN, BOTLIN, NUMLIN, NUMOVR)
C=======================================================================

C   --*** GRYCEN *** (GRPLIB) Find center line of text area
C   --   Written by Amy Gilkey - revised 02/20/87
C   --
C   --GRYCEN finds the center of a section of screen and returns the
C   --top and bottom line coordinates for the given number of lines.
C   --
C   --Parameters:
C   --   CHLSIZ - IN - the size of a line of text
C   --   TOPLIN, BOTLIN - IN/OUT - the device coordinates of the bottom of
C   --      the top and bottom lines of text
C   --   NUMLIN - IN/OUT - the number of lines requested, reduced by
C   --      NUMOVR if too long
C   --   NUMOVR - OUT - if the number of lines requested is greater than
C   --      the number of lines allowed, NUMOVR is the number of lines
C   --      deleted

      MAXLIN = INT((TOPLIN - BOTLIN) / CHLSIZ + 1 + 0.25)
      IF (NUMLIN .GT. MAXLIN) THEN
         NUMOVR = NUMLIN - MAXLIN
         NUMLIN = MAXLIN
      ELSE
         NUMOVR = 0
      END IF
      CEN = 0.5 * (TOPLIN + BOTLIN)
      TOPLIN = CEN + 0.5 * (NUMLIN-1) * CHLSIZ
      BOTLIN = TOPLIN - (NUMLIN-1) * CHLSIZ

      RETURN
      END
