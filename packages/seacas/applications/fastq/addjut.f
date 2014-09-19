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

C $Id: addjut.f,v 1.1 1990/11/30 11:02:52 gdsjaar Exp $
C $Log: addjut.f,v $
C Revision 1.1  1990/11/30 11:02:52  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]ADDJUT.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE ADDJUT (MXND, XN, YN, LXK, KXL, NXL, LXN,
     +   ANGLE, LNODES, XNEW, YNEW, NNN, LLL, NOLD, NLOOP, JUTTED)
C***********************************************************************
C
C  SUBROUTINE ADDJUT = ADDS A NEW NODE JUTTING OUT FROM AN EXISTING
C                      NODE
C
C***********************************************************************
C
      DIMENSION XN (MXND), YN (MXND)
      DIMENSION LXK (4, MXND), KXL (2, 3*MXND)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)
      DIMENSION ANGLE (MXND), LNODES (7, MXND)
C
      LOGICAL JUTTED
C
      NNN = NNN+1
      XN (NNN) = XNEW
      YN (NNN) = YNEW
C
C  MAKE LXN AND NXL ARRAYS
C
C  ADD THE NEW NODE'S LINES
C
      LLL = LLL+1
      NXL (1, LLL) = NNN
      NXL (2, LLL) = NOLD
C
      DO 100 I = 1, 4
         LXN (I, NNN) = 0
  100 CONTINUE
C
      KXL (1, LLL) = 0
      KXL (2, LLL) = 0
C
C  REDO THE LNODES ARRAY
C
      LNODES (1, NNN) = 0
      LNODES (2, NNN) = NOLD
      LNODES (3, NNN) = NOLD
      LNODES (4, NNN) = - 1
      LNODES (5, NNN) = LLL
C
      LNODES (1, NOLD) = 0
      LNODES (3, NOLD) = NNN
      LNODES (5, NOLD) = LLL
C
      NLOOP = NLOOP + 2
      JUTTED = .TRUE.
C
      RETURN
C
      END
