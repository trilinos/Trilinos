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

C $Id: chkkxl.f,v 1.1 1990/11/30 11:04:34 gdsjaar Exp $
C $Log: chkkxl.f,v $
C Revision 1.1  1990/11/30 11:04:34  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]CHKKXL.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE CHKKXL (MXND, LXK, KXL, LLL, ERR)
C***********************************************************************
C
C  SUBROUTINE CHKKXL = CHECKS TO SEE IF THE KXL COMPARES CORRECTLY TO
C                      THE LXK ARRAY
C
C***********************************************************************
C
      DIMENSION LXK (4, MXND), KXL (2, 3 * MXND)
C
      LOGICAL ERR
C
      ERR = .TRUE.
C
      DO 130 L = 1, LLL
         DO 120 IK = 1, 2
            K = KXL (IK, L)
            IF (K .NE. 0) THEN
               DO 100 I = 1, 4
                  IF (LXK (I, K) .EQ. L)GOTO 110
  100          CONTINUE
               WRITE ( * , 10000)IK, L, K
               RETURN
            ENDIF
  110       CONTINUE
  120    CONTINUE
  130 CONTINUE
      ERR = .FALSE.
C
      RETURN
C
10000 FORMAT ('KXL(', I4, ',', I4,') = ', I4,
     &   ' IS NOT IN LXK ARRAY  -  CHKKXL')
C
      END
