C    Copyright (c) 2014, Sandia Corporation.
C    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C    the U.S. Governement retains certain rights in this software.
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

C $Id: sortia.f,v 1.1 1990/11/30 11:16:11 gdsjaar Exp $
C $Log: sortia.f,v $
C Revision 1.1  1990/11/30 11:16:11  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]SORTIA.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE SORTIA(N,IAR,IRANGE,I)
C***********************************************************************
C
C  SUBROUTINE SORTIA = THIS SUBROUTINE SORTS AN INTEGER ARRAY IN
C                      ASCENDING ORDER
C
C***********************************************************************
C
C VARIABLES   IN : N ...   NUMBER OF ELEMENTS IN THE ARRAY
C                  IAR ... INTEGER ARRAY WITH DATA TO BE SORTED
C             OUT: I ...   INDEX ARRAY WITH ITS SORTED VALUES IN
C                          ASCENDING ORDER.
C
C WRITTEN BY:  HORACIO RECALDE                 DATE: FEB 25,  1988
C***********************************************************************
C
      INTEGER I(N),IAR(N)
C
C-- COPY ELEMENTS IN THE I ARRAY
C
      DO 100 J = 1,IRANGE
         I(J) = IAR(J)
  100 CONTINUE
C
C---  PERFORM AN EXCHANGE SORT ON THE FIRST IRANGE-1
C
      DO 120 K = 1,IRANGE - 1
         MIN = I(K)
C
C---  EXCHANGE THE K-TH ELEMENTS WITH THE MINIMUM ELEMENT REMAIN
C
         DO 110 J = K+1,IRANGE
            IF (I(J) .LT. MIN) THEN
               L = I(J)
               I(J) = I(K)
               I(K) = L
               MIN = I(K)
            ENDIF
  110    CONTINUE
  120 CONTINUE
C
      RETURN
      END
