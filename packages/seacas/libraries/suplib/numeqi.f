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
      INTEGER FUNCTION NUMEQI (INT, LENLST, INTLST)
C=======================================================================
C$Id: numeqi.f,v 1.2 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: numeqi.f,v $
CRevision 1.2  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:15:48  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:15:47  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:40  gdsjaar
c Initial revision
c 

C   --*** NUMEQI *** (ETCLIB) Count number of occurances of integer in list
C   --   Written by Amy Gilkey - revised 11/10/87
C   --
C   --NUMEQI returns the number of times the given integer occurs in a list
C   --of integers.
C   --
C   --Parameters:
C   --   INT - IN - the integer to be counted
C   --   LENLST - IN - the number of integers in the list
C   --   INTLST - IN - the list of integers to be searched

      INTEGER INT
      INTEGER LENLST
      INTEGER INTLST(*)

      NUMEQI = 0
      DO 10 I = 1, LENLST
         IF (INT .EQ. INTLST(I)) NUMEQI = NUMEQI + 1
   10 CONTINUE

      RETURN
      END
