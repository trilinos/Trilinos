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
      subroutine rwpval1(ntxt, ndb, flag, nump, nval, id, names,
     &  values, *)
C************************************************************************
C     This subroutine with read properties from an input text file and
C     write properties to an output exodusII file
C     ntxt   - IN - file id to read from
C     ndb    - IN - file id to write to
C     flag   - IN - flag indicating what type of properties to read/write
C     nump   - IN - number of properties
C     nval   - IN - number of values
C     names  - IN - names of the properties
C     values - IN - values of the properties

      include 'exodusII.inc'

      integer ntxt, ndb, flag, nump, nval
      integer id(*)
      character*(mxstln) names(*)
      integer values(*)

      if (nump .eq. 0) return

      do 100 i = 1, nump
C ... Skip comment line
        read (ntxt, *, end=130, err=130) 
        read (ntxt, '(A)', end=120, err=120) names(i)
C ... Skip comment line
        read (ntxt, *, end=130, err=130) 
        read (ntxt, *, end=130, err=130) (values(j), j = 1, nval)
        
C ... Write the property (unless it is 'ID')
        if (names(i) .ne. 'ID                              ') then
          do 90 j = 1, nval
            call expp(ndb, flag, id(j), names(i), values(j), ierr)
 90       continue
        end if
 100  continue
      return

 120  continue
 130  continue
      return 1
      end
