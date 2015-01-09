C    Copyright (c) 2014, Sandia Corporation.
C    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C    the U.S. Government retains certain rights in this software.
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C        * Redistributions of source code must retain the above copyright
C          notice, this list of conditions and the following disclaimer.
C    
C        * Redistributions in binary form must reproduce the above
C          copyright notice, this list of conditions and the following
C          disclaimer in the documentation and/or other materials provided
C          with the distribution.
C    
C        * Neither the name of Sandia Corporation nor the names of its
C          contributors may be used to endorse or promote products derived
C          from this software without specific prior written permission.
C    
C    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C    

C $Id: rotate.f,v 1.1 1990/11/30 11:15:08 gdsjaar Exp $
C $Log: rotate.f,v $
C Revision 1.1  1990/11/30 11:15:08  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]ROTATE.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE FQ_ROTATE (N, X, Y, NID, NEWF)
C***********************************************************************
C
C  SUBROUTINE ROTATE = CIRCULARLY SHIFTS THE DATA IN X,  Y,  AND NID
C
C***********************************************************************
C
      DIMENSION X (N), Y (N), NID (N)
C
      IF ((NEWF .LE. 1) .OR. (NEWF .GT. N)) RETURN
C
C  BUBBLE UP THROUGH THE ARRAYS AS MANY TIMES AS NEEDED
C
      DO 110 I = 1, NEWF - 1
         XLAST = X (1)
         YLAST = Y (1)
         NLAST = NID (1)
         DO 100 J = 1, N - 1
            X(J) = X (J + 1)
            Y(J) = Y (J + 1)
            NID(J) = NID (J + 1)
  100    CONTINUE
         X(N)   = XLAST
         Y(N)   = YLAST
         NID(N) = NLAST
  110 CONTINUE
C
      RETURN
C
      END
