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
      SUBROUTINE CPYSTR (LEN, IFROM, ITO)
C=======================================================================
C$Id: cpystr.f,v 1.2 2009/03/25 12:46:01 gdsjaar Exp $
C$Log: cpystr.f,v $
CRevision 1.2  2009/03/25 12:46:01  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:12:29  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:12:28  gdsjaar
c Initial revision
c
c Revision 1.1  90/08/09  13:39:07  gdsjaar
c Initial revision
c

C   --*** CPYSTR *** (STRLIB) Copy all strings in list
C   --   Written by Amy Gilkey - revised 03/15/88
C   --
C   --CPYSTR copies all the strings in a list to another list.
C   --
C   --Parameters:
C   --   LEN - IN - the number of strings in the list
C   --   IFROM - IN - the input list
C   --   ITO - OUT - the copied list

      INTEGER LEN
      CHARACTER*(*) IFROM(*), ITO(*)

      DO 100 I = 1, LEN
         ITO(I) = IFROM(I)
  100 CONTINUE

      RETURN
      END
