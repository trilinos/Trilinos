C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details
      subroutine putss(ndbout, numess, idess, neess, ndess, ixeess,
     *  ixdess, lteess, ltsess, facess, *)

      integer idess(*), neess(*), ndess(*)
      integer ixeess(*), ixdess(*)
      integer lteess(*), ltsess(*)
      real facess(*)

      if (numess .gt. 0) then

        call expcss(ndbout, idess, neess, ndess, ixeess, ixdess,
     &    lteess, ltsess, facess, ierr)
        if (ierr .gt. 0) return 1
      end if
      return
      end

