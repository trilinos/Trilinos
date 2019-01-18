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

C=======================================================================
      SUBROUTINE SMOGS2 (xscr, yscr, isbnd, x, y, numel,
     $     link, nlink, numnp, nscr)
C=======================================================================

C***********************************************************************
C
C  SUBROUTINE SMOGS  =  MESH SMOOTHING BY LAPLACE-S USING GAUSS-SEIDEL
C
C***********************************************************************

      real x(*), y(*), xscr(*), yscr(*)
      integer numel, nlink, link(nlink, *)
      logical isbnd(*)
      integer nscr(*)

      if (nlink .ne. 4) return
      do 20 iel = 1, numel
         n1 = link(1, iel)
         n2 = link(2, iel)
         n3 = link(3, iel)
         n4 = link(4, iel)

         if (.not. isbnd(n1)) then
            xscr(n1) = xscr(n1) + x(n2) + x(n4)
            yscr(n1) = yscr(n1) + y(n2) + y(n4)
            nscr(n1) = nscr(n1) + 2
         end if

         if (.not. isbnd(n2)) then
            xscr(n2) = xscr(n2) + x(n1) + x(n3)
            yscr(n2) = yscr(n2) + y(n1) + y(n3)
            nscr(n2) = nscr(n2) + 2
         end if

         if (.not. isbnd(n3)) then
            xscr(n3) = xscr(n3) + x(n2) + x(n4)
            yscr(n3) = yscr(n3) + y(n2) + y(n4)
            nscr(n3) = nscr(n3) + 2
         end if

         if (.not. isbnd(n4)) then
            xscr(n4) = xscr(n4) + x(n1) + x(n3)
            yscr(n4) = yscr(n4) + y(n1) + y(n3)
            nscr(n4) = nscr(n4) + 2
         end if

 20   continue
      return
      end
