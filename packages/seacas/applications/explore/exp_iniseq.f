C    Copyright(C) 2024 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      subroutine iniseq(icnt, map)
      integer map(*)
      do i=1, icnt
        map(i) = i
      end do
      return
      end
