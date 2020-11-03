C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      subroutine fixmap(numel, map)
      integer map(numel)

      do 10 i=1, numel
        if (map(i) .eq. 0) map(i) = i
 10   continue
      return
      end
