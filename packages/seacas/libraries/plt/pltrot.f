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

C $Id: pltrot.f,v 1.2 1993/07/16 22:27:49 gdsjaar Exp $ 
C $Log: pltrot.f,v $
C Revision 1.2  1993/07/16 22:27:49  gdsjaar
C Unrolled do loops to speed up execution.
C
c Revision 1.1  1993/07/16  16:49:18  gdsjaar
c Changed plt to library rather than single source file.
c 
C=======================================================================
      SUBROUTINE PLTROT(R,A,B)
      REAL R(3,3),A(3,3),B(3,3)

      B(1,1) = R(1,1)*A(1,1) + R(2,1)*A(2,1) + R(3,1)*A(3,1)
      B(2,1) = R(1,2)*A(1,1) + R(2,2)*A(2,1) + R(3,2)*A(3,1)
      B(3,1) = R(1,3)*A(1,1) + R(2,3)*A(2,1) + R(3,3)*A(3,1)
      
      B(1,2) = R(1,1)*A(1,2) + R(2,1)*A(2,2) + R(3,1)*A(3,2)
      B(2,2) = R(1,2)*A(1,2) + R(2,2)*A(2,2) + R(3,2)*A(3,2)
      B(3,2) = R(1,3)*A(1,2) + R(2,3)*A(2,2) + R(3,3)*A(3,2)
      
      B(1,3) = R(1,1)*A(1,3) + R(2,1)*A(2,3) + R(3,1)*A(3,3)
      B(2,3) = R(1,2)*A(1,3) + R(2,2)*A(2,3) + R(3,2)*A(3,3)
      B(3,3) = R(1,3)*A(1,3) + R(2,3)*A(2,3) + R(3,3)*A(3,3)

      RETURN

      END
