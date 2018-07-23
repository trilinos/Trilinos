C Copyright(C) 2011-2017 National Technology & Engineering Solutions of
C Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C * Redistributions of source code must retain the above copyright
C    notice, this list of conditions and the following disclaimer.
C           
C * Redistributions in binary form must reproduce the above
C   copyright notice, this list of conditions and the following
C   disclaimer in the documentation and/or other materials provided
C   with the distribution.
C                         
C * Neither the name of NTESS nor the names of its
C   contributors may be used to endorse or promote products derived
C   from this software without specific prior written permission.
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

c=======================================================================
      subroutine muntt (nelblk0, nelblk, nvarel, itt0, itt, istat)
C=======================================================================
C   --   ISTAT - IN - the status of each block:
C   --      0 = same
C   --      - = delete
C   --      n = combine with block n

      implicit none
      integer nelblk0, nelblk, nvarel
      integer istat(*)
      logical itt0 (nelblk0, nvarel)
      logical itt  (nelblk,  nvarel)
      
      integer ielb, ielb0, ivar
      
      ielb = 0
      do 10 ielb0 = 1, nelblk0
         if (istat(ielb0) .eq. 0) then
            ielb = ielb + 1
            do 20 ivar = 1, nvarel
               itt(ielb,ivar) = itt0(ielb0, ivar)
C               write (*,*) ielb0, ivar,itt0(ielb0, ivar)
 20         continue
         else if (istat(ielb0) .lt. 0) then
C ... Zero out the original truth table since we don't need to 
C     read those variables off of the input database...
            do 30 ivar = 1, nvarel
               itt0(ielb0,ivar) = .false.
 30         continue
            
         else
            call prterr('PROGRAM',
     *       'Element block combination does not work if'
     *       //' there are element variables on the database.')
            stop
         end if
 10   continue
      if (ielb .ne. nelblk) then
         call prterr('PROGRAM', 'Problem in muntt')
         stop
      end if
      return
      end
