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
      subroutine dbonam(ndbout, ndim, nameco,
     &  nvargl, ixgv, nvarnp, ixnv, nvarel, ixev,
     *  nvarns, ixnsv, nvarss, ixssv, names)
C=======================================================================
      include 'gp_namlen.blk'

      character*(maxnam) nameco(ndim)
      character*(maxnam) names(*)

      if (ndim .gt. 0) then
        call expcon(ndbout, nameco, ierr)
      end if
      if (nvargl .gt. 0) then
         call expvp(ndbout, 'g', nvargl, ierr)
         call expvan(ndbout, 'g', nvargl, names(ixgv), ierr)
      end if

      if (nvarnp .gt. 0) then
         call expvp(ndbout, 'n', nvarnp, ierr)
         call expvan(ndbout, 'n', nvarnp, names(ixnv), ierr)
      end if

      if (nvarel .gt. 0) then
         call expvp(ndbout, 'e', nvarel, ierr)
         call expvan(ndbout, 'e', nvarel, names(ixev), ierr)
      end if

      if (nvarns .gt. 0) then
         call expvp(ndbout, 'm', nvarns, ierr)
         call expvan(ndbout, 'm', nvarns, names(ixnsv), ierr)
      end if

      if (nvarss .gt. 0) then
         call expvp(ndbout, 's', nvarss, ierr)
         call expvan(ndbout, 's', nvarss, names(ixssv), ierr)
      end if
      return
      end
