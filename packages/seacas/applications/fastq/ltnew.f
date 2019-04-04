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

C $Id: ltnew.f,v 1.1 1990/11/30 11:11:40 gdsjaar Exp $
C $Log: ltnew.f,v $
C Revision 1.1  1990/11/30 11:11:40  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]LTNEW.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE LTNEW (MDIM, LINK)
C***********************************************************************
C
C  SUBROUTINE LTNEW = LOOKUP TABLE CLEARING FOR DATA POINTER ARRAYS
C
C***********************************************************************
C
C  VARIABLES USED:
C     MDIM   = DIMENSION OF LINK ARRAY, AND BASE FOR LOOKUP START
C     LINK   = LOOKUP TABLE ARRAY OF ID'S AND POINTERS
C              LINK(1,I) = ID VALUE STORED IN I'TH ROW (0 IF EMPTY)
C              LINK(2,I) = DATA POINTER ASSOCIATED W/THIS I'TH ID VALUE
C
C***********************************************************************
C
      DIMENSION LINK(2,MDIM)
C
C  ZERO OUT THE ID ROW ONLY
C
      DO 100 I = 1, MDIM
         LINK(1,I) = 0
  100 CONTINUE
C
      RETURN
C
      END
