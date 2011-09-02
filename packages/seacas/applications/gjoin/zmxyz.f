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

C=======================================================================
      SUBROUTINE ZMXYZ (NDIM, NUMNP, IXNP, XN, YN, ZN)
C=======================================================================
C $Id: zmxyz.f,v 1.1 1999/01/18 19:21:28 gdsjaar Exp $
C $Log: zmxyz.f,v $
C Revision 1.1  1999/01/18 19:21:28  gdsjaar
C ExodusII version of gjoin, needs testing and syncing with exodus 1 version, but is being committed to permit easier testing and modifications.  This was created by Dave Fry at Goodyear
C
c Revision 1.1.1.1  1998/11/05  16:23:28  a294617
c Initial import == gjoin 1.36
c
C Revision 1.1.1.1  1990/11/12 14:36:31  gdsjaar
C GJOIN - X1.00.40 - 7/17/90
C
c Revision 1.1  90/11/12  14:36:30  gdsjaar
c Initial revision
c 

C   --*** ZMXYZ *** (GJOIN) Compress coordinates
C   --   Written by Amy Gilkey - revised 01/20/88
C   --
C   --ZMXYZ compresses the coordinates by removing deleted nodes.
C   --
C   --Parameters:
C   --   NDIM - IN - the number of coordinates per node
C   --   NUMNP - IN/OUT - the number of nodes; returned compressed
C   --   IXNP - IN - the index of the compressed node; 0 if deleted
C   --   XN, YN, ZN - IN/OUT - the coordinates, returned compressed

      INTEGER IXNP(*)
      REAL XN(*), YN(*), ZN(*)

      JNP = 0
      DO 100 INP = 1, NUMNP
         IF (IXNP(INP) .GT. 0) THEN
            JNP = JNP + 1
            XN(JNP) = XN(INP)
            YN(JNP) = YN(INP)
            IF (NDIM .GE. 3) ZN(JNP) = ZN(INP)
         END IF
  100 CONTINUE

      NUMNP = JNP

      END
