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
      SUBROUTINE NUM2IX (NSETS, NINSET, IXSET)
C=======================================================================
C$Id: num2ix.f,v 1.2 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: num2ix.f,v $
CRevision 1.2  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:15:45  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:15:44  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:39  gdsjaar
c Initial revision
c 

C   --*** NUM2IX *** (ETCLIB) Change number-in-set to set indices
C   --   Written by Amy Gilkey - revised 10/23/87
C   --
C   --NUM2IX creates the set indices given the number in each set.
C   --
C   --Parameters:
C   --   NSETS - IN - the number of sets
C   --   NINSET - IN - the number in each set
C   --   IXSET - OUT - the set indices; set i = IXSET(i-1)+1 .. IXSET(i)

      INTEGER NINSET(*)
      INTEGER IXSET(0:*)

      IXSET(0) = 0
      DO 10 I = 1, NSETS
         IXSET(I) = IXSET(I-1) + NINSET(I)
   10 CONTINUE

      RETURN
      END
