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
      subroutine dobnd (x, y, z, numelb, nlink, idelb, link,
     $     cosin, nodes, isbnd, numnp, ndim, nelblk)
C=======================================================================
      real x(*), y(*), z(*)
      integer numelb(*), nlink(*), idelb(*), link(*)
      real cosin(ndim,*)
      integer nodes(*)
      logical isbnd(*)

      ielnk = 0
      iel1  = 0
      do 5 i=1, numnp
         isbnd(i) = .FALSE.
 5    continue

      do 100 iblk = 1, nelblk

      do 10 i=1, numnp
         cosin(1,i) = 0.0
         cosin(2,i) = 0.0
         if (ndim .eq. 3) cosin(3,i) = 0.0
 10   continue

      islnk = ielnk+1
      ielnk = ielnk + nlink(iblk) * numelb(iblk)

C ... ISBND(node) = .TRUE. if boundary node at return from donrm*

      if (ndim .eq. 2) then
         call donrm2( nodes, isbnd, cosin, numelb(iblk), link(islnk),
     $        x, y, nlink(iblk), numnp)
      else
         call donrm3( nodes, isbnd, cosin, numelb(iblk), link(islnk),
     $        x, y, z, nlink(iblk), numnp)
      end if
      iel1 = iel1 + numelb(iblk)
 100  continue
      return
      end
