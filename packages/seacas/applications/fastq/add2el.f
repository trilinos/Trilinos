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

C $Id: add2el.f,v 1.1 1990/11/30 11:02:44 gdsjaar Exp $
C $Log: add2el.f,v $
C Revision 1.1  1990/11/30 11:02:44  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]ADD2EL.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE ADD2EL (MXND, MLN, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &   ANGLE, LNODES, NNN, KKK, LLL, NNNOLD, LLLOLD, NLOOP, I1, I4,
     &   IAVAIL, NAVAIL, GRAPH, VIDEO, NOROOM, ERR)
C***********************************************************************
C
C  SUBROUTINE ADD2EL = CLOSES A SIX SIDED REGION BY FORMING 2 ELEMENTS
C
C***********************************************************************
C
      DIMENSION XN (MXND), YN (MXND), NUID (MXND)
      DIMENSION LXK (4, MXND), KXL (2, 3*MXND)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)
      DIMENSION ANGLE (MXND), LNODES (MLN, MXND)
C
      LOGICAL GRAPH, VIDEO, ERR, NOROOM
C
      ERR = .FALSE.
C
      I2 = LNODES (3, I1)
      I3 = LNODES (3, I2)
      I5 = LNODES (3, I4)
      I6 = LNODES (3, I5)
C
      CALL CONNOD (MXND, MLN, XN, YN, NUID, LXK, KXL, NXL, LXN, ANGLE,
     &   LNODES, NNN, KKK, LLL, NNNOLD, LLLOLD, I1, I2, I3, I4, IDUM,
     &   NLOOP, IAVAIL, NAVAIL, GRAPH, VIDEO, NOROOM, ERR)
      IF ((NOROOM) .OR. (ERR)) GOTO 100
      CALL CLOSE4 (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   I1, I4, I5, I6, KKK, ERR)
C
  100 CONTINUE
C
      RETURN
C
      END
