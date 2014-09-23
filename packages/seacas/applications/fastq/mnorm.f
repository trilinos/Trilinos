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

C $Id: mnorm.f,v 1.1 1990/11/30 11:12:25 gdsjaar Exp $
C $Log: mnorm.f,v $
C Revision 1.1  1990/11/30 11:12:25  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]MNORM.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE MNORM (MXND, XN, YN, NXL, LLL, STDLEN)
C***********************************************************************
C
C  SUBROUTINE MNORM = FINDS THE AVERAGE LENGTH OF THOSE LINES NOT MUCH
C                     LONGER THAN THE AVERAGE
C
C***********************************************************************
C
      DIMENSION NXL (2, 3 * MXND), XN (MXND), YN (MXND)
C
      STDLEN = 0.
      NUML = 0
      S = 0.0
      DO 100 L = 1, LLL
         N1 = NXL (1, L)
         N2 = NXL (2, L)
         IF (N1 .GT. 0) THEN
            D = SQRT ((XN (N1) - XN (N2)) **2 + (YN (N1) - YN (N2)) **2)
            S = S + D
            NUML = NUML + 1
         ENDIF
  100 CONTINUE
C
      IF (NUML .LE. 0) RETURN
      TOL = 1.25 * S / FLOAT (NUML)
      NUML = 0
      S = 0.0
      DO 110 L = 1, LLL
         N1 = NXL (1, L)
         N2 = NXL (2, L)
         IF (N1 .GT. 0) THEN
            D = SQRT ((XN (N1) - XN (N2)) **2 + (YN (N1) - YN (N2)) **2)
            IF (D .LT. TOL) THEN
               S = S + D
               NUML = NUML + 1
            ENDIF
         ENDIF
  110 CONTINUE
      STDLEN = S / FLOAT (NUML)
C
      RETURN
C
      END
