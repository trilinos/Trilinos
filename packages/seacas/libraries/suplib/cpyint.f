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
      SUBROUTINE CPYINT (LEN, IFROM, ITO)
C=======================================================================
C   --*** CPYINT *** (ETCLIB) Copy all integers in list
C   --   Written by Amy Gilkey - revised 11/03/87
C   --
C   --CPYINT copies all the integers in a list to another list.
C   --
C   --Parameters:
C   --   LEN - IN - the number of integers in the list
C   --   IFROM - IN - the input list
C   --   ITO - OUT - the copied list

      INTEGER LEN
      INTEGER IFROM(*), ITO(*)

      DO 100 I = 1, LEN-7,8
         ITO(I+0) = IFROM(I+0)
         ITO(I+1) = IFROM(I+1)
         ITO(I+2) = IFROM(I+2)
         ITO(I+3) = IFROM(I+3)
         ITO(I+4) = IFROM(I+4)
         ITO(I+5) = IFROM(I+5)
         ITO(I+6) = IFROM(I+6)
         ITO(I+7) = IFROM(I+7)
  100 CONTINUE
      do 110 J = I, LEN
         ITO(J) = IFROM(J)
 110  continue

      RETURN
      END
