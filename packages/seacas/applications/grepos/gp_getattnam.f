C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details
      subroutine getattnam(ndb, nelblk, id, numatr, names)
      include 'gp_namlen.blk'
      integer ndb, nelblk
      integer id(*)
      integer numatr(*)
      character*(maxnam) names(*)

      ibeg = 1
      do iel = 1, nelblk
        if (numatr(iel) .gt. 0) then
          call exgean(ndb, id(iel), numatr(iel), names(ibeg), ierr)
        end if
        ibeg = ibeg + numatr(iel)
      end do

      return
      end

      subroutine putattnam(ndb, nelblk, id, numatr, names)
      include 'gp_namlen.blk'
      integer ndb, nelblk
      integer id(*)
      integer numatr(*)
      character*(maxnam) names(*)

      ibeg = 1
      do iel = 1, nelblk
        if (numatr(iel) .gt. 0) then
          call expean(ndb, id(iel), numatr(iel), names(ibeg), ierr)
        end if
        ibeg = ibeg + numatr(iel)
      end do

      return
      end
