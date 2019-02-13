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
      SUBROUTINE DONRM2 (NODES, ISBND, COSIN, NUMEL, LINK, X, Y,
     $     NLINK, NUMNP)
C=======================================================================
      integer NODES(*)
      logical isbnd(*)
      REAL X(*), Y(*), COSIN(2,*)
      INTEGER LINK(NLINK, *)
      parameter (tol = 1.0e-6)
C
      do 10 inod = 1, numnp
         nodes(inod) = 0
 10   continue

      do 30 ilnk = 1, nlink
         do 20 iel = 1, numel
            nodes(link(ilnk, iel)) = 1
 20      continue
 30   continue

      do 60 iel = 1, numel
         do 50 iseg=1,nlink
            if (iseg .eq. nlink) then
               isegp1 = 1
            else
               isegp1 = iseg + 1
            end if

            XI = x( link(ISEG,iel) )
            YI = y( link(ISEG,iel) )
C
            XJ = x( link(ISEGp1,iel) )
            YJ = y( link(isegp1,iel) )
C
            DX = XI - XJ
            DY = YI - YJ
            RMAG = SQRT ( DX**2 + DY**2)
C
            cosin(1, link(iseg, iel)) = cosin(1, link(iseg, iel)) -
     $           dy / rmag
            cosin(2, link(iseg, iel)) = cosin(2, link(iseg, iel)) +
     $           dx / rmag

            cosin(1, link(isegp1, iel)) = cosin(1, link(isegp1, iel)) -
     $           dy / rmag
            cosin(2, link(isegp1, iel)) = cosin(2, link(isegp1, iel)) +
     $           dx / rmag
 50      continue
 60   continue

C
      do 70 inod = 1, numnp
         if (nodes(inod) .ne. 0) then
            if (abs(cosin(1, inod)) .gt. tol .or.
     $           abs(cosin(2, inod)) .gt. tol) then
               nodes(inod) = -1
               isbnd(inod) = .TRUE.
            end if
         end if
 70   continue

      RETURN
      END
