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

C $Id: spaced.f,v 1.1 1990/11/30 11:16:18 gdsjaar Exp $
C $Log: spaced.f,v $
C Revision 1.1  1990/11/30 11:16:18  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]SPACED.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE SPACED (MXND, MXCORN, MLN, ILEN, NCORN, LCORN, LNODES,
     &   ICOMB, ITEST, LTEST, ERR)
C***********************************************************************
C
C  SUBROUTINE SPACED = COUNTS THE INTERVAL SPACINGS FOR A COMBINATION
C
C***********************************************************************
C
      DIMENSION ICOMB (MXCORN), LCORN (MXCORN)
      DIMENSION ITEST (ILEN), LTEST (ILEN)
      DIMENSION LNODES (MLN, MXND)
C
      LOGICAL ERR
C
      ERR = .TRUE.
      KLEN = 0
      KOUNTC = 0
C
      DO 100 I = 1, NCORN
C
         IF (ICOMB (I) .EQ. 1) THEN
            KLEN = KLEN + 1
            IF (KLEN .GT. ILEN) THEN
               CALL MESAGE ('PROBLEMS IN SPACED - COUNTERS DON''T '//
     &            'MATCH DATA')
               RETURN
            ENDIF
C
            ITEST (KLEN) = LCORN(I)
            LTEST (KLEN) = KOUNTC
            KOUNTC = LNODES (7, LCORN(I))
         ELSE
            KOUNTC = KOUNTC + LNODES (7, LCORN (I))
         ENDIF
  100 CONTINUE
C
C  NOW ADD THE REMAINING KOUNTC ONTO THE FRONT
C
      LTEST (1) = LTEST (1) + KOUNTC
C
C  NOW SWITCH THE COUNTS TO BE FOLLOWING THE CORNERS INSTEAD OF
C  BEFORE THE CORNERS
C
      IHOLD = LTEST (1)
      DO 110 I = 2, KLEN
         LTEST (I - 1) = LTEST (I)
  110 CONTINUE
      LTEST (KLEN) = IHOLD
C
      ERR = .FALSE.
      RETURN
C
      END
