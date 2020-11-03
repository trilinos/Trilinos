C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

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
