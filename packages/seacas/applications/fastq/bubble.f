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

C $Id: bubble.f,v 1.1 1990/11/30 11:04:13 gdsjaar Exp $
C $Log: bubble.f,v $
C Revision 1.1  1990/11/30 11:04:13  gdsjaar
C Initial revision
C
CC* FILE: [.QMESH]BUBBLE.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE BUBBLE (X, KARRY, NORD, N)
C***********************************************************************
C
C  SUBROUTINE BUBBLE=SORTS ALL VALUES X(I), KARRY(I) INTO DECREASING
C                      ORDER, ASSUMING THAT VALUES 1 TO NORD ARE SORTED
C
C***********************************************************************
C
      DIMENSION X (N), KARRY (N)
C
      IF (N .LE. 1) RETURN
C
      ISTART = MAX0 (NORD + 1, 2)
      IF (ISTART .GT. N) RETURN
      DO 120 J = ISTART, N
         XVAL = X (J)
         KVAL = KARRY (J)
         JM1 = J - 1
         I = J
         DO 100 II = 1, JM1
            IF  (XVAL .LE. X (I - 1)) GO TO 110
            X (I) = X (I - 1)
            KARRY (I) = KARRY (I - 1)
            I = I - 1
  100    CONTINUE
  110    CONTINUE
         X (I) = XVAL
         KARRY (I) = KVAL
  120 CONTINUE
C
      RETURN
C
      END
