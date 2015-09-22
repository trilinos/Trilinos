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

C $Id: mxcopy.f,v 1.3 1993/07/19 18:08:44 gdsjaar Exp $ 
C $Log: mxcopy.f,v $
C Revision 1.3  1993/07/19 18:08:44  gdsjaar
C Added special case for n=4 since that is how plt calls it primarily
C
c Revision 1.2  1993/07/16  19:35:45  gdsjaar
c Restructured to optimize faster
c
c Revision 1.1  1993/07/16  16:47:33  gdsjaar
c Changed plt to library rather than single source file.
c 
C=======================================================================
      SUBROUTINE MXCOPY(N,MAT1,MAT2)
      REAL MAT1(N,*),MAT2(N,*)
      
      IF (N .EQ. 4) THEN
        DO 100 I=1, 4
          MAT2(I,1) = MAT1(I,1)
          MAT2(I,2) = MAT1(I,2)
          MAT2(I,3) = MAT1(I,3)
          MAT2(I,4) = MAT1(I,4)
 100    CONTINUE
      ELSE
        DO 2910 I = 1,N
          DO 2890 J = 1,N
            MAT2(I,J) = MAT1(I,J)
 2890     CONTINUE
 2910   CONTINUE
      end if
      RETURN
      
      END
