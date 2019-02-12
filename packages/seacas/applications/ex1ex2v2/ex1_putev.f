C Copyright (C) 2009-2017 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
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
C     * Neither the name of NTESS nor the names of its
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

      subroutine putev (idexo, nwstep, nelblk, nvarel, numelb,
     &     varel, idelb, isevok, ierr)
C======================================================================
C -- *** PUTEV *** Put Element Variables in regular netCDF file
C --
C --PUTEV calls the exodus II interface routine that writes the
C --      element variable values into the regular netCDF file.
C --
C --Parameters:
C --   idexo  - IN  - EXODUS file ID returned from a previous call to
C --                  EXCRE or EXOPEN.
C --   nwstep - IN  - The time step number.
C --   nelblk - IN  - The number of element blocks.
C --   nvarel - IN  - The number of element variables.
C --   numelb - IN  - An array containing the number of elements per
C --                  element block.
C --   varel  - IN  - An array containing the element variables.
C --   idelb  - IN  - Array of element block IDs
C --   isevok - IN  - Element variable truth table
C --   ierr   - OUT - Returned error code.  If no errors occurred, 0
C --                  is returned.

      include 'exodusII.inc'
      integer numelb(*), nwstep, nelblk, nvarel, idelb(*), ierr
      integer isevok(nvarel,*)
      real varel(*)

      ielo = 1
      do 200 ielb = 1, nelblk
        do 100 ivar = 1, nvarel
          if (isevok(ivar,ielb) .ne. 0) then
             call expev (idexo, nwstep, ivar, idelb(ielb), numelb(ielb),
     &            varel(ielo), ierr)
             if (ierr .lt. 0) then
      		call exerr ('putev','Error calling expev', exlmsg)
             endif
             ielo = ielo + numelb(ielb)
          endif
100     continue
200   continue
      return
      end
