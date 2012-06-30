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

C $Id: mxmult.f,v 1.3 1993/07/16 22:50:49 gdsjaar Exp $ 
C $Log: mxmult.f,v $
C Revision 1.3  1993/07/16 22:50:49  gdsjaar
C Unrolled loops for faster execution
C
c Revision 1.2  1993/07/16  19:14:19  gdsjaar
c Restructured to optimize faster
c
c Revision 1.1  1993/07/16  16:47:35  gdsjaar
c Changed plt to library rather than single source file.
c 
C=======================================================================
      SUBROUTINE MXMULT(N,MAT1,MAT2,MATR)
      REAL MAT1(N,*),MAT2(N,*),MATR(N,*)

      IF (N .EQ. 4) THEN
          MATR(1,1) = MAT1(1,1)*MAT2(1,1) + MAT1(1,2)*MAT2(2,1) +
     *                MAT1(1,3)*MAT2(3,1) + MAT1(1,4)*MAT2(4,1)
          MATR(1,2) = MAT1(1,1)*MAT2(1,2) + MAT1(1,2)*MAT2(2,2) +
     *                MAT1(1,3)*MAT2(3,2) + MAT1(1,4)*MAT2(4,2)
          MATR(1,3) = MAT1(1,1)*MAT2(1,3) + MAT1(1,2)*MAT2(2,3) +
     *                MAT1(1,3)*MAT2(3,3) + MAT1(1,4)*MAT2(4,3)
          MATR(1,4) = MAT1(1,1)*MAT2(1,4) + MAT1(1,2)*MAT2(2,4) +
     *                MAT1(1,3)*MAT2(3,4) + MAT1(1,4)*MAT2(4,4)

          MATR(2,1) = MAT1(2,1)*MAT2(1,1) + MAT1(2,2)*MAT2(2,1) +
     *                MAT1(2,3)*MAT2(3,1) + MAT1(2,4)*MAT2(4,1)
          MATR(2,2) = MAT1(2,1)*MAT2(1,2) + MAT1(2,2)*MAT2(2,2) +
     *                MAT1(2,3)*MAT2(3,2) + MAT1(2,4)*MAT2(4,2)
          MATR(2,3) = MAT1(2,1)*MAT2(1,3) + MAT1(2,2)*MAT2(2,3) +
     *                MAT1(2,3)*MAT2(3,3) + MAT1(2,4)*MAT2(4,3)
          MATR(2,4) = MAT1(2,1)*MAT2(1,4) + MAT1(2,2)*MAT2(2,4) +
     *                MAT1(2,3)*MAT2(3,4) + MAT1(2,4)*MAT2(4,4)

          MATR(3,1) = MAT1(3,1)*MAT2(1,1) + MAT1(3,2)*MAT2(2,1) +
     *                MAT1(3,3)*MAT2(3,1) + MAT1(3,4)*MAT2(4,1)
          MATR(3,2) = MAT1(3,1)*MAT2(1,2) + MAT1(3,2)*MAT2(2,2) +
     *                MAT1(3,3)*MAT2(3,2) + MAT1(3,4)*MAT2(4,2)
          MATR(3,3) = MAT1(3,1)*MAT2(1,3) + MAT1(3,2)*MAT2(2,3) +
     *                MAT1(3,3)*MAT2(3,3) + MAT1(3,4)*MAT2(4,3)
          MATR(3,4) = MAT1(3,1)*MAT2(1,4) + MAT1(3,2)*MAT2(2,4) +
     *                MAT1(3,3)*MAT2(3,4) + MAT1(3,4)*MAT2(4,4)

          MATR(4,1) = MAT1(4,1)*MAT2(1,1) + MAT1(4,2)*MAT2(2,1) +
     *                MAT1(4,3)*MAT2(3,1) + MAT1(4,4)*MAT2(4,1)
          MATR(4,2) = MAT1(4,1)*MAT2(1,2) + MAT1(4,2)*MAT2(2,2) +
     *                MAT1(4,3)*MAT2(3,2) + MAT1(4,4)*MAT2(4,2)
          MATR(4,3) = MAT1(4,1)*MAT2(1,3) + MAT1(4,2)*MAT2(2,3) +
     *                MAT1(4,3)*MAT2(3,3) + MAT1(4,4)*MAT2(4,3)
          MATR(4,4) = MAT1(4,1)*MAT2(1,4) + MAT1(4,2)*MAT2(2,4) +
     *                MAT1(4,3)*MAT2(3,4) + MAT1(4,4)*MAT2(4,4)

      ELSE
        
        DO 230 K = 1,N
          DO 200 J = 1,N
            MATR(K,J) = 0.0
 200      CONTINUE
          DO 220 I = 1,N
            DO 210 J = 1,N
              MATR(K,J) = MATR(K,J) + MAT1(K,I)*MAT2(I,J)
 210        CONTINUE
 220      CONTINUE
 230    CONTINUE
      END IF
      RETURN
      END
