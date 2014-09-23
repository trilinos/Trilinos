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

C $Id: lowfnd.f,v 1.1 1990/11/30 11:11:28 gdsjaar Exp $
C $Log: lowfnd.f,v $
C Revision 1.1  1990/11/30 11:11:28  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]LOWFND.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE LOWFND (MXND, NUID, N, INDX, I, IOLD)
C***********************************************************************
C
C  SUBROUTINE LOWFND  =  LINEAR INDEXED SEARCH FOR MATCHING NUID VALUES
C
C***********************************************************************
C
      DIMENSION NUID(MXND), INDX(N)
C
      IBOT = 1
      ITOP = N
C
  100 CONTINUE
      II = (IBOT + ITOP)/2
      IF (NUID(INDX(II)) .EQ. NUID(I)) THEN
         IOLD = INDX(II)
         RETURN
      ELSE IF (NUID(INDX(II)) .GT. NUID(I)) THEN
         ITOP = II - 1
      ELSE
         IBOT = II + 1
      ENDIF
      IF (IBOT .LE. ITOP) GO TO 100
C
      IOLD = 0
      RETURN
C
      END
