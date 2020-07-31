C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details
      subroutine chkss(nelblk, numelb, iblock)

      INTEGER NUMELB(*)
      INTEGER iblock(*)

      itot = 0
      do 10 i=1, nelblk
         iblock(i) = itot + numelb(i)
         itot = itot + numelb(i)
 10   continue
      return
      end

