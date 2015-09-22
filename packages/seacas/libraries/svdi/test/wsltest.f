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

      PROGRAM WTEST05
      REAL X(102), Y(102)
      CHARACTER COLOR(7)*4
C
      DATA COLOR/'RED','GREE','YELL','BLUE','MAGE','CYAN','WHIT'/
C
      CALL WSTART(0.,0)
      CALL WTTYPE('SOFT')
      DO 10 I=1,102
         X(I) = I/100.*2.*3.141592
         Y(I) = SIN(X(I))
   10 CONTINUE
      CALL WXAXIS(.2,.2,.7,0.,13.,1.,'Sines of the Times')
      CALL WYAXIS(.2,.2,.5,-1.,1.,.5,' ')
      CALL WLNSTY('DOTT')
      CALL WGRID(0,1)
      CALL WLNSTY('SOLI')
      CALL WMARK('DIAM')
      DO 20 I=1,7
         CALL WSPACE('NDC')
         CALL WMSIZE(I/500.)
         CALL WSPACE('USER')
         CALL WCOLOR(COLOR(I))
         CALL WDRAW(X,Y,102,10)
         DO 15 J=1,102
            X(J) = X(J) + 1.
   15    CONTINUE
   20 CONTINUE
      CALL WCOLOR(COLOR(3))
      CALL WEND
      STOP
      END
