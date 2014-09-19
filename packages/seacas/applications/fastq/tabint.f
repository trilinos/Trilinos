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

C $Id: tabint.f,v 1.1 1990/11/30 11:17:04 gdsjaar Exp $
C $Log: tabint.f,v $
C Revision 1.1  1990/11/30 11:17:04  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]TABINT.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE TABINT (X1, X2, Y1, Y2, CT, ST, SCALE, XX1, YY1, XX2,
     &   YY2, DRWTAB)
C***********************************************************************
C
C     SUBROUTINE TABINT = INITIALIZES THE TABLET TO THE PLOT LIMITS
C
C***********************************************************************
C
      LOGICAL DRWTAB
C
      IF (DRWTAB) THEN
         THETA = ATAN2 (YY2 - YY1, XX2 - XX1) -
     &      ATAN2 (Y2 - Y1, X2 - X1)
         CT = COS (THETA)
         ST = SIN (THETA)
         SCALE = SQRT (((X2 - X1) ** 2 + (Y2 - Y1) ** 2 ) /
     &      ((XX2 - XX1) ** 2 + (YY2 - YY1) ** 2 ))
      ELSE
         CT = 1.
         ST = 0.
         XX1 = 2000
         XX2 = 15000
         YY1 = 2000
         YY2 = 10000
         SCALEX =  (X2 - X1) / (XX2 - XX1)
         SCALEY =  (Y2 - Y1) / (YY2 - YY1)
         IF (SCALEX .GT. SCALEY) THEN
            SCALE = SCALEX
            YY1 =  (YY2 - YY1) -  ( (Y2 - Y1) / SCALE)
         ELSE
            SCALE = SCALEY
            XX1 =  (XX2 - XX1) -  ( (X2 - X1) / SCALE)
         ENDIF
      ENDIF
C
      RETURN
C
      END
