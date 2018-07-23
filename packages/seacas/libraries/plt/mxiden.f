C Copyright (C) 2009-2017 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
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
C     * Neither the name of NTESS nor the names of its
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

C $Id: mxiden.f,v 1.3 1993/07/19 18:08:45 gdsjaar Exp $ 
C $Log: mxiden.f,v $
C Revision 1.3  1993/07/19 18:08:45  gdsjaar
C Added special case for n=4 since that is how plt calls it primarily
C
c Revision 1.2  1993/07/16  19:35:47  gdsjaar
c Restructured to optimize faster
c
c Revision 1.1  1993/07/16  16:47:34  gdsjaar
c Changed plt to library rather than single source file.
c 
C=======================================================================
      SUBROUTINE MXIDEN(N,MAT)
      REAL MAT(N,*)

      IF (N .EQ. 4) THEN
        MAT(1,1) = 1.0
        MAT(1,2) = 0.0
        MAT(1,3) = 0.0
        MAT(1,4) = 0.0
        
        MAT(2,1) = 0.0
        MAT(2,2) = 1.0
        MAT(2,3) = 0.0
        MAT(2,4) = 0.0
        
        MAT(3,1) = 0.0
        MAT(3,2) = 0.0
        MAT(3,3) = 1.0
        MAT(3,4) = 0.0
        
        MAT(4,1) = 0.0
        MAT(4,2) = 0.0
        MAT(4,3) = 0.0
        MAT(4,4) = 1.0
      ELSE
        DO 2930 I = 1,N
          DO 2950 J = 1,N
            MAT(I,J) = 0.
 2950     CONTINUE
 2930   CONTINUE
        DO 2970 I = 1,N
          MAT(I,I) = 1.
 2970   CONTINUE
      end if
      RETURN
      
      END
