C Copyright(C) 2009-2017 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
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
C
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

C=======================================================================
      SUBROUTINE FFADDV (IRANGE, LINE)
C=======================================================================
C$Id: ffaddv.f,v 1.2 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: ffaddv.f,v $
CRevision 1.2  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:14:18  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:14:17  gdsjaar
c Initial revision
c
c Revision 1.1  90/08/09  13:39:23  gdsjaar
c Initial revision
c

C   --*** FFADDV *** (FFLIB) Add integer range to line
C   --   Written by Amy Gilkey - revised 11/16/87
C   --
C   --FFADDV adds an integer range (as a character string) to a line.
C   --
C   --Parameters:
C   --   IRANGE - IN - the integer range to add;
C   --      IRANGE(1) TO IRANGE(2) BY IRANGE(3)
C   --   LINE - IN/OUT - the line being built

      INTEGER IRANGE(3)
      CHARACTER*(*) LINE

      CALL FFADDI (IRANGE(1), LINE)
      IF (IRANGE(2) .NE. IRANGE(1)) THEN
         CALL FFADDC ('TO', LINE)
         CALL FFADDI (IRANGE(2), LINE)
      END IF
      IF (IRANGE(3) .NE. 1) THEN
         CALL FFADDC ('BY', LINE)
         CALL FFADDI (IRANGE(3), LINE)
      END IF

      RETURN
      END
