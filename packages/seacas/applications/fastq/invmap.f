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

C $Id: invmap.f,v 1.1 1990/11/30 11:10:26 gdsjaar Exp $
C $Log: invmap.f,v $
C Revision 1.1  1990/11/30 11:10:26  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]INVMAP.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE INVMAP (X0, Y0, X1, Y1, X2, Y2, X3, Y3, X4, Y4, SXI,
     &   SETA, INSIDE)
C***********************************************************************
C
C  THIS IS A TEST OF THE INVERTED MAPPING OF AN ELEMENT
C
C***********************************************************************
C
      DOUBLE PRECISION AX, BX, CX, DX, AY, BY, CY, DY
      DOUBLE PRECISION ALPHA, BETA, GAMMA, RAD
      DOUBLE PRECISION XI, ETA, XI1, ETA1, XI2, ETA2
C
      LOGICAL INSIDE
C
      EPS = 1.E-3
      EPS2 = 1.E-10
C
C  GET THE A, B, C, AND D VALUES FOR X AND Y.
C
      AX = X1 - X0
      BX = X2 - X1
      CX = X1 - X2 + X3 -X4
      DX = X4 - X1
      AY = Y1 - Y0
      BY = Y2 - Y1
      CY = Y1 - Y2 + Y3 -Y4
      DY = Y4 - Y1
C
C  CALCULATE THE ALPHA, BETA, AND GAMMA VALUES.
C
      ALPHA = (CY * DX) - (CX * DY)
      BETA  = (AX * CY) - (AY * CX) + (BY * DX) - (BX * DY)
      GAMMA = (AX * BY) - (AY * BX)
C
C  CALCULATE THE XI AND ETA VALUES.
C
      IF (ALPHA .EQ. 0.) THEN
         ETA = -GAMMA / BETA
         IF ((ETA .EQ. 0) .AND. (BX .EQ. 0)) THEN
            XI = (Y0 - Y1) / (Y2 - Y1)
         ELSE IF ((BX .EQ. -CX) .AND. (ETA .EQ. 1.)) THEN
            XI = (Y0 - Y3)/(Y4 - Y3)
         ELSE IF (BX .EQ. (-CX * ETA)) THEN
            XI = -1000.
         ELSE
            XI = (- AX - (DX * ETA)) / (BX + (CX * ETA))
         ENDIF
      ELSE
         RAD = BETA**2 - (4. * ALPHA * GAMMA)
         IF (RAD .LT. 0.) THEN
C
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/19/90
CC* MODIFICATION: COMMENTED OUT THE ERROR MESSAGE FOR THE
C**               NEGATIVE RADICAL PROBLEM AS IT APPEARS THAT
C**               THIS MAY OCCUR - IT JUST MEANS THAT THE POINT
C**               TRULLY IS NOT IN THE ELEMENT.
C
C            CALL MESAGE ('** ERROR - NEGATIVE RADICAL IN INVMAP **')
            INSIDE = .FALSE.
            GOTO 100
         ENDIF
         RAD = DSQRT (RAD)
         ETA1 = (- BETA + RAD) / (2. * ALPHA)
         ETA2 = (- BETA - RAD) / (2. * ALPHA)
C
         IF ((ABS(ETA1) .LT. EPS2) .AND. (ABS(BX) .LT. EPS2)) THEN
            XI1 = (Y0 - Y1) / (Y2 - Y1)
         ELSE IF ((BX .EQ. -CX) .AND. (ETA1 .EQ. 1.)) THEN
            XI1 = (Y0 - Y3)/(Y4 - Y3)
         ELSE IF (BX .EQ. (-CX * ETA1)) THEN
            XI1 = -1000.
         ELSE
            XI1 = (- AX - (DX * ETA1)) / (BX + (CX * ETA1))
         ENDIF
C
         IF ((ABS(ETA2) .LT. EPS2) .AND. (ABS(BX) .LT. EPS2)) THEN
            XI2 = (Y0 - Y1) / (Y2 - Y1)
         ELSE IF ((BX .EQ. -CX) .AND. (ETA2 .EQ. 1.)) THEN
            XI2 = (Y0 - Y3)/(Y4 - Y3)
         ELSE IF (BX .EQ. (-CX * ETA2)) THEN
            XI2 = -1000.
         ELSE
            XI2 = (- AX - (DX * ETA2)) / (BX + (CX * ETA2))
         ENDIF
C
         D1 = DSQRT (ETA1*ETA1 + XI1*XI1)
         D2 = DSQRT (ETA2*ETA2 + XI2*XI2)
         IF (D1 .LT. D2) THEN
            ETA = ETA1
            XI = XI1
         ELSE
            ETA = ETA2
            XI = XI2
         ENDIF
      END IF
C
C  CHECK TO SEE IF ETA AND XI ARE WITHIN THE ELEMENT
C
      IF (.NOT. ((ETA .LE. 1.0 + EPS) .AND.
     &   (ETA .GE. 0.0 - EPS)) ) THEN
         INSIDE = .FALSE.
         GOTO 100
      ELSE IF (.NOT. ((XI .LE. 1.0 + EPS) .AND.
     &   (XI .GE. 0.0 - EPS)) ) THEN
         INSIDE = .FALSE.
         GOTO 100
      ELSE
         INSIDE = .TRUE.
      ENDIF
      SXI = XI
      SETA = ETA
C
  100 CONTINUE
      RETURN
C
      END
