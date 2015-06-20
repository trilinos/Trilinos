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

C $Id: l2smo.f,v 1.1 1990/11/30 11:10:53 gdsjaar Exp $
C $Log: l2smo.f,v $
C Revision 1.1  1990/11/30 11:10:53  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]L2SMO.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE L2SMO (MXND, XN, YN, NXL, LXN, NNN, NNNOLD, NIT, EPS,
     &   RO)
C***********************************************************************
C
C  SUBROUTINE L2SMO = LAPLACIAN SQUARED METHOD OF MESH SMOOTHING
C
C***********************************************************************
C
C  NOTE:
C     THIS METHOD IS LIKE LAPLACIAN SMOOTHING EXCEPT EACH VECTOR
C     POINTING TO A NEIGHBOR NODE HAS A LENGTH OF
C      (DISTANCE TO THAT NODE)**2.
C
C***********************************************************************
C
C  VARIABLES USED:
C     NIT = MAX NUMBER OF ITERATIONS TO DO
C     EPS = NODE MOVEMENT TOLERANCE FOR CONVERGENCE
C     RO  = UNDER OR OVER-RELAXATION FACTOR.
C
C***********************************************************************
C
      DIMENSION LINES (20)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND), XN (MXND), YN (MXND)
C
      LOGICAL BIG, ERR
C
      EPS2 = (EPS * RO) ** 2
C
C  ITERATION LOOP
C
      DO 120 IT = 1, NIT
         BIG = .FALSE.
C
C  NODE LOOP
C
         DO 110 NODE = NNNOLD + 1, NNN
            IF ((LXN (1, NODE).GT.0) .AND. (LXN (2, NODE).GT.0))THEN
               FX = 0.
               FY = 0.
               SL = 0.
C
C  LOOP THRU ALL LINES CONNECTED TO NODE
C
               CALL GETLXN (MXND, LXN, NODE, LINES, KOUNT, ERR)
               IF (ERR) RETURN
               DO 100 IL = 1, KOUNT
                  L = LINES (IL)
                  NEND = NXL (1, L) + NXL (2, L) - NODE
                  DX = XN (NEND) - XN (NODE)
                  DY = YN (NEND) - YN (NODE)
                  AL2 = DX * DX + DY * DY
                  AL = SQRT (AL2)
                  FX = FX + DX * AL
                  FY = FY + DY * AL
                  SL = SL + AL
  100          CONTINUE
C
C  MOVE THE NODE
C
               DELX = RO * FX/SL
               DELY = RO * FY/SL
               XN (NODE) = XN (NODE) + DELX
               YN (NODE) = YN (NODE) + DELY
               IF (DELX ** 2 + DELY ** 2 .GT. EPS2) BIG = .TRUE.
            ENDIF
  110    CONTINUE
         IF (.NOT.BIG)RETURN
  120 CONTINUE
      RETURN
      END
