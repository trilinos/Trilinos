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

C $Id: esolve.f,v 1.1 1990/11/30 11:06:51 gdsjaar Exp $
C $Log: esolve.f,v $
C Revision 1.1  1990/11/30 11:06:51  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]ESOLVE.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      REAL FUNCTION ESOLVE (A7, A8, A2, ANG1, ANG2)
C***********************************************************************
C
C  FUNCTION ESOLVE = FINDS A SOLUTION TO THE ELIPSE EQUATION
C                    GIVEN AN INTERVAL THAT CONTAINS THE SOLUTION
C
C***********************************************************************
C
      EPS = 1.E-6
C
      F1 = ELIPSE (A7, A8, A2, ANG1)
      IF (ABS(F1) .LT. EPS) THEN
         ESOLVE = ANG1
         GO TO 110
      END IF
      F2 = ELIPSE (A7, A8, A2, ANG2)
      IF (ABS(F2) .LT. EPS) THEN
         ESOLVE = ANG2
         GO TO 110
      END IF
C
  100 CONTINUE
      IF (ABS(ANG1 - ANG2) .LT. EPS) THEN
         ESOLVE = (ANG1 + ANG2)/2.0
      ELSE
         ANG3 = (ANG1 + ANG2)/2.0
         F3 = ELIPSE (A7, A8, A2, ANG3)
C
         IF (ABS(F3) .LT. EPS) THEN
            ESOLVE = ANG3
            GO TO 110
         ELSE IF (F1/F3 .LT. 0.0) THEN
            ANG2 = ANG3
            F2 = F3
         ELSE
            ANG1 = ANG3
            F1 = F3
         END IF
         GO TO 100
      END IF
C
  110 CONTINUE
C
      RETURN
C
      END
