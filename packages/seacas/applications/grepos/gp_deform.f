C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

      subroutine deform(x, y, z, numnp, ndim, dx, dy, dz, ndb, idefst)
      real x(*), y(*), z(*), dx(*), dy(*), dz(*)

      if (numnp .le. 0) return

C ... Read the displacements from the database.
C     Assume they are the first 'ndim' nodal variables...
      call exgnv (ndb, idefst, 1, numnp, dx, ierr)
      if (ndim .ge. 2) then
        call exgnv (ndb, idefst, 2, numnp, dy, ierr)
      end if
      if (ndim .eq. 3) then
        call exgnv (ndb, idefst, 3, numnp, dz, ierr)
      end if

C ... Deform the variables...
      do i=1, numnp
        x(i) = x(i) + dx(i)
      end do

      if (ndim .ge. 2) then
        do i=1, numnp
          y(i) = y(i) + dy(i)
        end do
      end if

      if (ndim .eq. 3) then
        do i=1, numnp
          z(i) = z(i) + dz(i)
        end do
      end if

      return
      end
