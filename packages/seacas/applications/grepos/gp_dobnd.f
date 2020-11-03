C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      subroutine dobnd (x, y, z, numelb, nlink, idelb, link,
     $     cosin, nodes, isbnd, numnp, ndim, nelblk)
C=======================================================================
      real x(*), y(*), z(*)
      integer numelb(*), nlink(*), idelb(*), link(*)
      real cosin(ndim,*)
      integer nodes(*)
      logical isbnd(*)

      ielnk = 0
      iel1  = 0
      do 5 i=1, numnp
         isbnd(i) = .FALSE.
 5    continue

      do 100 iblk = 1, nelblk

      do 10 i=1, numnp
         cosin(1,i) = 0.0
         cosin(2,i) = 0.0
         if (ndim .eq. 3) cosin(3,i) = 0.0
 10   continue

      islnk = ielnk+1
      ielnk = ielnk + nlink(iblk) * numelb(iblk)

C ... ISBND(node) = .TRUE. if boundary node at return from donrm*

      if (ndim .eq. 2) then
         call donrm2( nodes, isbnd, cosin, numelb(iblk), link(islnk),
     $        x, y, nlink(iblk), numnp)
      else
         call donrm3( nodes, isbnd, cosin, numelb(iblk), link(islnk),
     $        x, y, z, nlink(iblk), numnp)
      end if
      iel1 = iel1 + numelb(iblk)
 100  continue
      return
      end
