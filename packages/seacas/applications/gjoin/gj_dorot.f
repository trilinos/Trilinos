C Copyright (c) 2008 Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C 
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.
C 
C     * Neither the name of Sandia Corporation nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
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
C 

C $Id: dorot.f,v 1.1 1999/01/18 19:21:21 gdsjaar Exp $
C $Log: dorot.f,v $
C Revision 1.1  1999/01/18 19:21:21  gdsjaar
C ExodusII version of gjoin, needs testing and syncing with exodus 1 version, but is being committed to permit easier testing and modifications.  This was created by Dave Fry at Goodyear
C
c Revision 1.1.1.1  1998/11/05  16:23:25  a294617
c Initial import == gjoin 1.36
c
C Revision 1.1  1992/11/11 22:00:09  gdsjaar
C Added revolve and revcen transformation commands.
C
      
      subroutine dorot(ndim, numnp, xn, yn, zn, rotmat, rotcen)
      real xn(*), yn(*), zn(*)
      real rotmat(3,3), rotcen(3)

      if (NDIM .EQ. 3) THEN
         DO 10 JNP = 1, NUMNP
            X = XN(JNP) - ROTCEN(1)
            Y = YN(JNP) - ROTCEN(2)
            Z = ZN(JNP) - ROTCEN(3)
            XN(JNP) = X*ROTMAT(1,1) + Y*ROTMAT(2,1) + Z*ROTMAT(3,1)
     $           + ROTCEN(1)
            YN(JNP) = X*ROTMAT(1,2) + Y*ROTMAT(2,2) + Z*ROTMAT(3,2)
     $           + ROTCEN(2)
            ZN(JNP) = X*ROTMAT(1,3) + Y*ROTMAT(2,3) + Z*ROTMAT(3,3)
     $           + ROTCEN(3)
 10         continue
         ELSE IF (NDIM .EQ. 2) THEN
            DO 20 JNP = 1, NUMNP
               X = XN(JNP) - ROTCEN(1)
               Y = YN(JNP) - ROTCEN(2)
               XN(JNP) = X*ROTMAT(1,1) + Y*ROTMAT(2,1) + ROTCEN(1)
               YN(JNP) = X*ROTMAT(1,2) + Y*ROTMAT(2,2) + ROTCEN(2)
 20         continue
         END IF
         return
         end
