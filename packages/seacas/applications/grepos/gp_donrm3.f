C Copyright(C) 2011-2017 National Technology & Engineering Solutions of
C Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C * Redistributions of source code must retain the above copyright
C    notice, this list of conditions and the following disclaimer.
C           
C * Redistributions in binary form must reproduce the above
C   copyright notice, this list of conditions and the following
C   disclaimer in the documentation and/or other materials provided
C   with the distribution.
C                         
C * Neither the name of NTESS nor the names of its
C   contributors may be used to endorse or promote products derived
C   from this software without specific prior written permission.
C                                                 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

C     1 - 12x15  2 - 26x21  3 - 37x32  4 - 48x43
C       - 14x12    - 23x26    - 34x37    - 41x48
C       - 15x14    - 21x23    - 32x34    - 43x41
C        
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
C
            do 40 ino = 1, 3
               XJ = x( link(map(ino,iseg,1),iel) ) - xi
               YJ = y( link(map(ino,iseg,1),iel) ) - yi
               ZJ = z( link(map(ino,iseg,1),iel) ) - zi
C
               XK = x( link(map(ino,iseg,2),iel) ) - xi
               YK = y( link(map(ino,iseg,2),iel) ) - yi
               ZK = z( link(map(ino,iseg,2),iel) ) - zi
C
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
      
C
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
