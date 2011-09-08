C Copyright (c) 2007 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
C retains certain rights in this software.
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
C 

C************************************************************************
      subroutine rwpval(ndbi, ndbo, flag, nump, nval, names, values,
     *  namlen)
C************************************************************************
C     This subroutine with read properties from an input file and
C     write properties to an output text file
C     ndbi   - IN - file id to read from
C     ndbo   - IN - file id to write to
C     flag   - IN - flag indicating what type of properties to read/write
C     nump   - IN - number of properties
C     nval   - IN - number of values
C     names  - IN - names of the properties
C     values - IN - values of the properties

      integer ndbi, ndbo, flag, nump, nval
      character*(namlen) names(*)
      integer values(*)

      if (nump .eq. 0) return

      call exgpn (ndbi, flag, names, ierr)

      do 100 i = 1, nump
         call exgpa (ndbi, flag, names(i), values, ierr)
         write (ndbo,1000) '! Property Name: '
         write (ndbo,1000) names(i)
         write (ndbo,1000) '! Property Value(s): '
         write (ndbo,1010) (values(j), j = 1, nval)
 100  continue

 1000 format (A)
 1010 format (7(I10, 1X))
      return
      end
