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

C $Id: lpntin.f,v 1.1 1990/11/30 11:11:30 gdsjaar Exp $
C $Log: lpntin.f,v $
C Revision 1.1  1990/11/30 11:11:30  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]LPNTIN.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      LOGICAL FUNCTION LPNTIN (MAXNP, CX, CY, NPNT, X, Y)
C***********************************************************************
C
C  FUNCTION LPNTIN = .TRUE. IF THE POINT IS WITHIN THE PERIMETER
C
C***********************************************************************
C
      LOGICAL IN
      REAL CX(MAXNP), CY(MAXNP), X, Y
C
      I1 = 1
      IN = .TRUE.
  100 CONTINUE
      IF (I1 .LE. NPNT .AND. IN) THEN
         IF (I1 .EQ. NPNT) THEN
            I2 = 1
         ELSE
            I2 = I1 + 1
         END IF
C
C  CHECK FOR VERTICAL LINE
         IF (ABS(CX(I1) - CX(I2)) .LT. 1.0E-06) THEN
            IF (CY(I1) .LT. CY(I2)) THEN
               IN = X .LT. CX(I2)
            ELSE
               IN = X .GT. CX(I2)
            END IF
C
C  CHECK FOR HORIZONTAL LINE
         ELSE IF (ABS(CY(I1) - CY(I2)) .LT. 1.0E-06) THEN
            IF (CX(I1) .LT. CX(I2)) THEN
               IN = Y .GT. CY(I2)
            ELSE
               IN = Y .LT. CY(I2)
            END IF
C
C  MUST BE INCLINED LINE
         ELSE
            U = (X - CX(I1))/(CX(I2) - CX(I1))
            U = MIN(U, 1.0)
            U = MAX(0.0, U)
            V = (Y - CY(I1))/(CY(I2) - CY(I1))
            V = MIN(V, 1.0)
            V = MAX(0.0, V)
            IF (CX(I1) .LT. CX(I2)) THEN
               IF (CY(I1) .LT. CY(I2)) THEN
                  IN = U .LE. V
               ELSE
                  IN = U .GE. V
               END IF
            ELSE
               IF (CY(I1) .GT. CY(I2)) THEN
                  IN = U .LE. V
               ELSE
                  IN = U .GE. V
               END IF
            END IF
         END IF
C
         I1 = I1 + 1
         GO TO 100
      END IF
      LPNTIN = IN
C
      RETURN
      END
