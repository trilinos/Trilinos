C Copyright (C) 2009 Sandia Corporation.  Under the terms of Contract
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

C $Id: mxvect.f,v 1.3 1993/07/16 22:56:19 gdsjaar Exp $ 
C $Log: mxvect.f,v $
C Revision 1.3  1993/07/16 22:56:19  gdsjaar
C Unrolled loops for faster execution
C
c Revision 1.2  1993/07/16  19:30:48  gdsjaar
c Restructured to optimize faster
c
c Revision 1.1  1993/07/16  16:47:36  gdsjaar
c Changed plt to library rather than single source file.
c 
C=======================================================================
      SUBROUTINE MXVECT(N,VEC,MAT,RES)
      REAL VEC(*),MAT(N,*),RES(*)

      IF (N .EQ. 4) THEN
        RES(1) = MAT(1,1)*VEC(1) + MAT(2,1)*VEC(2) + MAT(3,1)*VEC(3) +
     *    MAT(4,1)*VEC(4)
        
        RES(2) = MAT(1,2)*VEC(1) + MAT(2,2)*VEC(2) + MAT(3,2)*VEC(3) +
     *    MAT(4,2)*VEC(4)
        
        RES(3) = MAT(1,3)*VEC(1) + MAT(2,3)*VEC(2) + MAT(3,3)*VEC(3) +
     *    MAT(4,3)*VEC(4)
        
        RES(4) = MAT(1,4)*VEC(1) + MAT(2,4)*VEC(2) + MAT(3,4)*VEC(3) +
     *    MAT(4,4)*VEC(4)
        
      ELSE
        DO 2980 J = 1,N
          RES(J) = 0.0
 2980   CONTINUE
        DO 3010 I = 1,N
          DO 2990 J = 1,N
            RES(J) = RES(J) + MAT(I,J)*VEC(I)
 2990     CONTINUE
 3010   CONTINUE
      END IF
      RETURN

      END
