C    Copyright(C) 2008 Sandia Corporation.  Under the terms of Contract
C    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C    certain rights in this software
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C    * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C              
C    * Redistributions in binary form must reproduce the above
C      copyright notice, this list of conditions and the following
C      disclaimer in the documentation and/or other materials provided
C      with the distribution.
C                            
C    * Neither the name of Sandia Corporation nor the names of its
C      contributors may be used to endorse or promote products derived
C      from this software without specific prior written permission.
C                                                    
C    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C    
C=======================================================================
      SUBROUTINE SELBLK (NUMSEL, IXSEL, NELBLK, LISBLK, NUMELB, NUMLNK,
     *  LINK, ISCR, NUMNP, EBTYPE)
C=======================================================================

      include 'params.blk'
      
      INTEGER IXSEL(*)
      INTEGER LISBLK(0:*)
      INTEGER NUMELB(*)
      INTEGER NUMLNK(*)
      INTEGER LINK(*)
      INTEGER ISCR(*)
      LOGICAL SELECT
      CHARACTER*(MXSTLN) EBTYPE(*)
      CHARACTER*40 STRA
      CHARACTER*132 MSG

      do 80 i=1, numnp
        iscr(i) = 0
 80   continue
      
      islnk = 1
      do 100 ielb = 1, nelblk
        select = .false.
        do 90 ix = 1, lisblk(0)
          if (ielb .eq. lisblk(ix)) then
            select = .true.
          end if
 90     continue
        if (ebtype(ielb) .eq. 'nsided' .or.
     *    ebtype(ielb) .eq. 'NSIDED') THEN
          numnod = numlnk(ielb)
        else
          numnod = numlnk(ielb) * numelb(ielb)
        end if
        if (select) then
          call selblk1(ielb, numnod, link(islnk), iscr)
        end if
        ISLNK = ISLNK + numnod
 100  CONTINUE
C
      numsel = 0
      do 120 i=1, numnp
        if (iscr(i) .gt. 0) then
          numsel = numsel + 1
          ixsel(numsel) = i
        end if
 120  continue

      write (stra, 10000) numsel
      call pckstr(1, stra)
      MSG = STRA(:lenstr(stra)) // ' nodes selected'
      call prterr('CMDSPEC', MSG)
10000 format(I10)
      return
      end

      subroutine selblk1(ielb, numnod, link, iscr)

      integer link(*)
      integer iscr(*)
      
      do i=1, numnod
        node = link(i)
        iscr(node) = iscr(node) + 1
      end do
      return
      end
