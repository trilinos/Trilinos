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

C $Id: indexi.f,v 1.1 1990/11/30 11:09:26 gdsjaar Exp $
C $Log: indexi.f,v $
C Revision 1.1  1990/11/30 11:09:26  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]INDEXI.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE INDEXI_FQ (NMAX, IARRAY, N, INDX)
C***********************************************************************
C
C     INDEXI - SORT THE N ELEMENTS OF IARRAY WHOSE POSITION IS STORED
C              IN INDX.  ONLY THE ORDER OF THE INDEX ARRAY, INDX, IS
C              MODIFIED
C
C***********************************************************************
C
      DIMENSION IARRAY(NMAX), INDX(N)
C
      L = N/2 + 1
      IR = N
  100 CONTINUE
      IF (L .GT. 1) THEN
         L = L - 1
         INDXT = INDX(L)
         JQ = IARRAY(INDXT)
      ELSE
         INDXT = INDX(IR)
         JQ = IARRAY(INDXT)
         INDX(IR) = INDX(1)
         IR = IR - 1
         IF (IR .EQ. 1) THEN
            INDX(1) = INDXT
            RETURN
         END IF
      END IF
C
      I = L
      J = L + L
  110 CONTINUE
      IF (J .LE. IR) THEN
         IF (J .LT. IR) THEN
            IF (IARRAY(INDX(J)) .LT. IARRAY(INDX(J + 1))) J = J + 1
         END IF
         IF (JQ .LT. IARRAY(INDX(J))) THEN
            INDX(I) = INDX(J)
            I = J
            J = J + J
         ELSE
            J = IR + 1
         END IF
         GO TO 110
      END IF
      INDX(I) = INDXT
      GO TO 100
C
      END
