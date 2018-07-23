C    Copyright(C) 2014-2017 National Technology & Engineering Solutions of
C    Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C    * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C    
C    * Redistributions in binary form must reproduce the above
C      copyright notice, this list of conditions and the following
C      disclaimer in the documentation and/or other materials provided
C      with the distribution.
C    
C    * Neither the name of NTESS nor the names of its
C      contributors may be used to endorse or promote products derived
C      from this software without specific prior written permission.
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

C $Id: nodord.f,v 1.1 1990/11/30 11:12:46 gdsjaar Exp $
C $Log: nodord.f,v $
C Revision 1.1  1990/11/30 11:12:46  gdsjaar
C Initial revision
C
C
CC* FILE: [.RENUM]NODORD.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE NODORD (NPNODE, XN, YN, LISTN, NUID, NNN)
C***********************************************************************
C
C  SUBROUTINE NODORD = ORDER THE NODE TABLE INTO INCREASING VALUES OF
C                      THE VARIABLE LISTN
C
C***********************************************************************
C
      DIMENSION LISTN (NPNODE)
      DIMENSION XN (NPNODE), YN (NPNODE), NUID (NPNODE)
C
      NN = NNN
      M = NN
  100 CONTINUE
      M =  (9 * M) / 16
      IF (M .LE. 0) RETURN
      M1 = M + 1
      DO 120 J = M1, NN
         L = J
         I = J - M
  110    CONTINUE
         IF (LISTN (L) .LT. LISTN (I)) THEN
            KLISTN = LISTN (I)
            KNUID = NUID (I)
            TXN = XN (I)
            TYN = YN (I)
            LISTN (I) = LISTN (L)
            NUID (I) = NUID (L)
            XN (I) = XN (L)
            YN (I) = YN (L)
            LISTN (L) = KLISTN
            NUID (L) = KNUID
            XN (L) = TXN
            YN (L) = TYN
            L = I
            I = I - M
            IF (I .GE. 1)GOTO 110
         ENDIF
  120 CONTINUE
      GOTO 100
C
      END
