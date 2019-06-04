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
      SUBROUTINE SMOGS3 (xscr, yscr, zscr, isbnd, x, y, z, numel,
     $     link, nlink, numnp, nscr)
C=======================================================================

C***********************************************************************
C
C  SUBROUTINE SMOGS  =  MESH SMOOTHING BY LAPLACE-S USING GAUSS-SEIDEL
C
C***********************************************************************

      real x(*), y(*), z(*), xscr(*), yscr(*), zscr(*)
      integer numel, nlink, link(nlink, *)
      logical isbnd(*)
      integer nscr(*)
      integer map(3,8)

      data map /2,4,5,  3,1,6,  4,2,7,  1,3,8,
     $          1,8,6,  2,5,7,  3,6,8,  4,7,5/

      if (nlink .ne. 8) return
      do 20 iel = 1, numel
         do 10 ino = 1, nlink

            ni = link(ino, iel)
            if (.not. isbnd(ni)) then
               n1 = link( map(1, ino), iel)
               n2 = link( map(2, ino), iel)
               n3 = link( map(3, ino), iel)
               xscr(ni) = xscr(ni) + x(n1) + x(n2) + x(n3)
               yscr(ni) = yscr(ni) + y(n1) + y(n2) + y(n3)
               zscr(ni) = zscr(ni) + z(n1) + z(n2) + z(n3)
               nscr(ni) = nscr(ni) + 3
            end if
 10      continue
 20   continue
      return
      end
