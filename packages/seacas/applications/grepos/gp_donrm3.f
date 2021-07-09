C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C     1 - 12x15  2 - 26x21  3 - 37x32  4 - 48x43
C       - 14x12    - 23x26    - 34x37    - 41x48
C       - 15x14    - 21x23    - 32x34    - 43x41

C     5 - 51x56  6 - 62x67  7 - 73x78  8 - 84x85
C       - 58x51    - 65x62    - 76x73    - 87x84
C       - 56x58    - 67x65    - 78x76    - 85x87
C=======================================================================
      SUBROUTINE DONRM3 (NODES, ISBND, COSIN, NUMEL, LINK, X, Y, Z,
     $     NLINK, NUMNP)
C=======================================================================
      integer NODES(*)
      logical isbnd(*)
      REAL X(*), Y(*), Z(*), COSIN(3,*)
      INTEGER LINK(NLINK, *)
      integer map(3,8,2)

      parameter (tol = 1.0e-6)

      data map /2,4,5,  6,3,1,  7,4,2,  8,1,3,
     $          1,8,6,  2,5,7,  3,6,8,  4,7,5,
     $          5,2,4,  1,6,3,  2,7,4,  3,8,1,
     $          6,1,8,  7,2,5,  8,3,6,  5,4,7/

      do 10 inod = 1, numnp
         nodes(inod) = 0
 10   continue

      do 30 ilnk = 1, nlink
         do 20 iel = 1, numel
            nodes(link(ilnk, iel)) = 1
 20      continue
 30   continue

      do 60 iel = 1, numel
         do 50 iseg=1,nlink
            XI = x( link(ISEG,iel) )
            YI = y( link(ISEG,iel) )
            ZI = z( link(ISEG,iel) )

            do 40 ino = 1, 3
               XJ = x( link(map(ino,iseg,1),iel) ) - xi
               YJ = y( link(map(ino,iseg,1),iel) ) - yi
               ZJ = z( link(map(ino,iseg,1),iel) ) - zi

               XK = x( link(map(ino,iseg,2),iel) ) - xi
               YK = y( link(map(ino,iseg,2),iel) ) - yi
               ZK = z( link(map(ino,iseg,2),iel) ) - zi

               ax = yj * zk - zj * yk
               ay = zj * xk - xj * zk
               az = xj * yk - yj * xk

               rmag = sqrt(ax**2 + ay**2 + az**2)

               cosin(1, link(iseg, iel)) = cosin(1, link(iseg, iel)) +
     $              ax / rmag
               cosin(2, link(iseg, iel)) = cosin(2, link(iseg, iel)) +
     $              ay / rmag
               cosin(3, link(iseg, iel)) = cosin(3, link(iseg, iel)) +
     $              az / rmag
 40         continue
 50      continue
 60   continue

      do 70 inod = 1, numnp
         if (nodes(inod) .ne. 0) then
            if ( abs(cosin(1, inod)) .gt. tol .or.
     $           abs(cosin(2, inod)) .gt. tol .or.
     $           abs(cosin(3, inod)) .gt. tol) then
               nodes(inod) = -1
               isbnd(inod) = .TRUE.
            end if
         end if
 70   continue

      RETURN
      END
