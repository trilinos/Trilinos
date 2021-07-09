C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

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
